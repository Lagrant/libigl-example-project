#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_face_normals.h>
#include <iostream>
#include <string>
#include "constraint.hpp"
#include "halfedge.h"

#define random(x) (rand()%x)

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd N_faces;
Eigen::Vector3d* d;
set<int>* labels;
int* result;
vertex* ver;
edge* e;
face* triFace;
int NfaceRows;
int faceNum;
int total;

enum class suffixs{
    OFF,
    OBJ
};

static std::map<std::string, suffixs> s_mapStringValues =
{
    {"off",suffixs::OFF},
    {"obj",suffixs::OBJ},
    {"OFF",suffixs::OFF},
    {"OBJ",suffixs::OBJ}
};

using namespace std;

void loadFile(Eigen::MatrixXd& V, Eigen::MatrixXi& F, string filename, string path = "/Users/Lagrant/libigl/tutorial/shared/"){
    string suffix = "";
    int dot = 0;
    while (dot < filename.size()) {
        if(filename[dot] == '.'){
            if(dot == filename.size()-1){
                cerr<<"incorrect file name";
                exit(-1);
            }
            while(dot < filename.size()-1){
                dot++;
                suffix += filename[dot];
            }
            break;
        } else dot++;
    }
    
    switch(s_mapStringValues[suffix])
    {
        case suffixs::OFF:
            igl::readOFF(path+filename,V,F);
            break;
        case suffixs::OBJ:
            igl::readOBJ(path+filename,V,F);
            break;
        default:
            cerr<<"invalid file name";
            exit(-1);
    }
}

int main(int argc, char *argv[])
{
    string path = "/Users/Lagrant/libigl/tutorial/shared/";
    string filename = "";
    cin>>filename;
    
    // Load a mesh in OFF format
    loadFile(V, F, filename);
//    igl::readOFF("/Users/Lagrant/libigl/tutorial/shared/bunny.off",V,F);
    
    int labelCost, totalVer, htotal;
    total       = 600;
    labelCost   = 0.5*10;
    if(total % 2)
        total--;
    htotal      = total/2;
    
    // Compute per-face normals
    igl::per_face_normals(V,F,N_faces);
    NfaceRows = N_faces.rows();
    faceNum   = F.rows();
    totalVer  = faceNum*3;
    
    Eigen::RowVector3d M        = V.colwise().maxCoeff().transpose();
    Eigen::RowVector3d* color   = new Eigen::RowVector3d[total];
    Eigen::MatrixXd C(faceNum,3);
    igl::viewer::Viewer viewer;
    queue<int> tracked;
    queue<int> untracked;
    vector<int> directionTracking;
    
    ver          = (vertex*) malloc(sizeof(vertex)*totalVer);
    e            = (edge*) malloc(sizeof(edge)*totalVer);
    triFace      = (face*) malloc(sizeof(face)*faceNum);
    result       = (int*) malloc(sizeof(int)*faceNum);
    int* mark    = (int*) malloc(sizeof(int)*faceNum);
    int* un_mark = (int*) malloc(sizeof(int)*faceNum);
    d            = new Eigen::Vector3d[total];
    memset(result, -1, sizeof(int)*faceNum);
    memset(mark, 0 , sizeof(int)*faceNum);
    memset(un_mark, 0 , sizeof(int)*faceNum);
    
    for(int i = 0;i < htotal;i++){
        d[i](0) = Guass();
        d[i](1) = Guass();
        d[i](2) = Guass();
        d[i].normalize();
        d[i+htotal] = -d[i];
    }
    
    for(int i = 0;i < total;i++){
        color[i](0) = (double) rand() / RAND_MAX;
        color[i](1) = (double) rand() / RAND_MAX;
        color[i](2) = (double) rand() / RAND_MAX;
    }
    
    //initialize the set class
    
    void* rawMemory = operator new(NfaceRows*sizeof(set<int>));
    labels = reinterpret_cast<set<int>*>(rawMemory);
    for(int i = 0;i < faceNum;i++)
        new (&labels[i])set<int>(total);
    
    printf("face number = %d\n", faceNum);
    
    
    // To apply the halfedge data structure
    getHalfEdge(&ver, &e, &triFace);
    
    build(&ver, totalVer);
    
    // To compute the connected components
    heightField();
    
    GeneralGraph_DArraySArraySpatVarying(labelCost);
    
    integrate();
    
    integrate();
    
//    directionTracking.push_back(triFace[0].label);
    append(viewer, directionTracking, triFace[0].label, M, color[triFace[0].label]);
    tracked.push(0);
    int preLabel = -1;
    do {
        if(tracked.empty()){
            while(!untracked.empty() && mark[untracked.front()])
                untracked.pop();
            if(untracked.empty())
                break;
            else {
                tracked.push(untracked.front());
                untracked.pop();
            }
        }
        
        int currentFaceNumber = tracked.front();
        tracked.pop();
        
        int currentFaceLabel  = triFace[currentFaceNumber].label;
        mark[currentFaceNumber] = 1;
        
        if(preLabel != currentFaceLabel){
            preLabel = currentFaceLabel;
            //change color[reuslt[i]] here
//            color[currentFaceLabel](0) *= 0.9;
//            color[currentFaceLabel](1) = (double)rand()/RAND_MAX;
//            color[currentFaceLabel](2) *= 0.7;
            color[currentFaceLabel].row(0) *= 0.7;
        }
        C.row(currentFaceNumber) = color[currentFaceLabel];
        
        edge* currentEdge = triFace[currentFaceNumber].adjacentEdge;
        edge* runEdge = currentEdge;
        
        do{
            int label = runEdge->pair->face->label;
            int numbering = runEdge->pair->face->numbering;
            
            if(!mark[numbering]){
                if(label == currentFaceLabel){
                    C.row(numbering) = color[label];
                    tracked.push(numbering);
                    mark[numbering] = 1;
                    runEdge = runEdge->next;
                    
                } else if(!un_mark[numbering]){
                    append(viewer, directionTracking, label, M, color[label]);
                    untracked.push(numbering);
                    un_mark[numbering] = 1;
                    runEdge = runEdge->next;
                    
                } else runEdge = runEdge->next;
                
            } else runEdge = runEdge->next;
        } while(currentEdge != runEdge);
        
    } while(!tracked.empty() || !untracked.empty());
    
    
    
    
//    for(int i = 0 ;i < faceNum;i++)
//        C.row(i) = color[result[i]];
    
    // Plot the mesh
    viewer.core.show_lines = false;
    viewer.data.set_mesh(V, F);
    viewer.data.set_normals(N_faces);
    viewer.data.set_colors(C);
//    viewer.data.add_points(M,Eigen::RowVector3d(0,0,0));
    
//    int componentNumber = directionTracking.size();
//    for(int i = 0; i < componentNumber; i++){
//        viewer.data.add_edges
//        (
//         M,
//         d[directionTracking.at(i)].transpose() + M,
//         color[directionTracking.at(i)]*0.7
//         );
//    }
    
    viewer.data.add_edges
    (
     M,
     Eigen::RowVector3d(M(0)+1, M(1), M(2)),
     Eigen::RowVector3d(0,0,0)
     );
    viewer.data.add_edges
    (
     M,
     Eigen::RowVector3d(M(0), M(1)+1, M(2)),
     Eigen::RowVector3d(0,0,0)
     );
    viewer.data.add_edges
    (
     M,
     Eigen::RowVector3d(M(0), M(1), M(2)+1),
     Eigen::RowVector3d(0,0,0)
     );
    viewer.data.add_edges
    (
     M,
     Eigen::RowVector3d(M(0)-1, M(1), M(2)),
     Eigen::RowVector3d(0,0,0)
     );
    viewer.data.add_edges
    (
     M,
     Eigen::RowVector3d(M(0), M(1)-1, M(2)),
     Eigen::RowVector3d(0,0,0)
     );
    viewer.data.add_edges
    (
     M,
     Eigen::RowVector3d(M(0), M(1), M(2)-1),
     Eigen::RowVector3d(0,0,0)
     );
    
    viewer.launch();
    
    for(int i = 0;i < faceNum;i++)
        labels[i].~set();
    
    operator delete(rawMemory);
    delete[] color;
    delete[] d;
    free(result);
    free(ver);
    free(e);
    free(triFace);
    free(mark);
    free(un_mark);
    
    return 0;
}

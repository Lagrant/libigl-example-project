#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_face_normals.h>
#include <iostream>
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
//connectedComponents* ccp;
int NfaceRows;
int faceNum;
int total;

using namespace std;

int main(int argc, char *argv[])
{
    // Load a mesh in OFF format
    igl::readOFF( "/Users/Lagrant/libigl/tutorial/shared/bunny.off", V, F);
    
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
    
    Eigen::Vector3d M           = V.colwise().maxCoeff();
    Eigen::RowVector3d* color   = new Eigen::RowVector3d[total];
    Eigen::MatrixXd C(faceNum,3);
    
    queue<int> tracked;
    queue<int> untracked;
    vector<int> directionTracking;
    
    ver          = (vertex*) malloc(sizeof(vertex)*totalVer);
    e            = (edge*) malloc(sizeof(edge)*totalVer);
    triFace      = (face*) malloc(sizeof(face)*faceNum);
    //  ccp          = (connectedComponents*) malloc(sizeof(connectedComponents)*faceNum*total);
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
        color[i].normalize();
    }
    
    //initialize the set class
    
    void* rawMemory = operator new(NfaceRows*sizeof(set<int>));
    labels = reinterpret_cast<set<int>*>(rawMemory);
    for(int i = 0;i < faceNum;i++)                     //NfaceRows = faceNum, NfaceRows means the total numebr of norms of faces while faceNum means the total  number of faces
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
    
    directionTracking.push_back(triFace[0].label);
    tracked.push(0);
    
    do {
        if(tracked.empty()){
            while(mark[untracked.front()] && !untracked.empty())
                untracked.pop();
            if(untracked.empty())
                break;
            else {
                tracked.push(untracked.front());
                untracked.pop();
            }
        }
        
        int currentFaceLabel = tracked.front();
        tracked.pop();
        mark[currentFaceLabel] = 1;
        
        edge* currentEdge = triFace[currentFaceLabel].adjacentEdge;
        edge* runEdge = currentEdge;
        
        do{
            int label = runEdge->pair->face->label;
            int numbering = runEdge->pair->face->numbering;
            if(!mark[numbering]){
                if(label == triFace[currentFaceLabel].label){
                    tracked.push(runEdge->pair->face->numbering);
                    mark[numbering] = 1;
                    runEdge = runEdge->next;
                } else if(!un_mark[numbering]){
                    insert(directionTracking, label);
                    untracked.push(runEdge->pair->face->numbering);
                    un_mark[numbering] = 1;
                    runEdge = runEdge->next;
                }
            }
        } while(currentEdge != runEdge);
        
    } while(!tracked.empty() && !untracked.empty());
    
    int componentNumber = directionTracking.size();
    
    Eigen::MatrixXd VIndicator(componentNumber+1,3);
    
    
    //change color[reuslt[i]] here
    
    for(int i = 0 ;i < faceNum;i++)
        C.row(i) = color[result[i]];
    
    
    
    // Plot the mesh
    igl::viewer::Viewer viewer;
    viewer.core.show_lines = false;
    viewer.data.set_mesh(V, F);
    viewer.data.set_normals(N_faces);
    viewer.data.set_colors(C);
    
    viewer.launch();
    
    for(int i = 0;i < faceNum;i++){
        
        labels[i].~set();
    }
    
    operator delete(rawMemory);
    delete[] color;
    delete[] d;
    free(result);
    free(ver);
    free(e);
    free(triFace);
    //    free(ccp);
    
    return 0;
}

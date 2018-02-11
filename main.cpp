#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_face_normals.h>
#include <Eigen/LU>
#include <iostream>
#include <vector>
#include "set.h"
#include "constraints.h"
#include "halfedge.h"

#define random(x) (rand()%x)

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd N_faces;
Eigen::Vector3d* d;
set<int>* components;
set<vector<int>>* edges;
set<int>* labels;
int* result;
vertex* ver;
edge* e;
face* triFace;
connectedComponents* ccp;

using namespace std;

int main(int argc, char *argv[])
{
    int total;
    total = 120;
    d = new Eigen::Vector3d[total];

    for(int i = 0;i < total;i++){
        d[i](0) = Guass();
        d[i](1) = Guass();
        d[i](2) = Guass();
        d[i].normalize();
    }
    
    // Load a mesh in OFF format
    igl::readOFF( "/Users/Lagrant/libigl/tutorial/shared/bunny1.off", V, F);
    
    // Compute per-face normals
    igl::per_face_normals(V,F,N_faces);
    int NfaceRows = N_faces.rows();
    int faceNum = F.rows();
    
    // To apply the halfedge data structure
    int totalVer = faceNum*3;
    ver = (vertex*) malloc(sizeof(vertex)*totalVer);
    e = (edge*) malloc(sizeof(edge)*totalVer);
    triFace = (face*) malloc(sizeof(face)*faceNum);
    ccp = (connectedComponents*) malloc(sizeof(connectedComponents)*faceNum);
    
    Eigen::RowVector3d* color = new Eigen::RowVector3d[total];
    result = (int*) malloc(sizeof(int)*faceNum);
    memset(result, -1, sizeof(int)*faceNum);
    
    //initialize the set class
    void* rawMemory = operator new(faceNum*sizeof(set<int>));
    void* rawMemory1 = operator new(faceNum*sizeof(set<vector<int>>));
    void* rawMemory3 = operator new(NfaceRows*sizeof(set<int>));
    
    components = reinterpret_cast<set<int>*>(rawMemory);
    edges = reinterpret_cast<set<vector<int>>*>(rawMemory1);
    labels = reinterpret_cast<set<int>*>(rawMemory3);
    
    for(int i = 0;i < faceNum;i++){
        new (&components[i])set<int>(faceNum); //NfaceRows = faceNum, NfaceRows means the total numebr of norms of faces while faceNum means the total  number of faces
        new (&edges[i])set<vector<int>>(3*faceNum);
        new (&labels[i])set<int>(faceNum);
    }
    printf("face number = %d\n", faceNum);

    for(int i = 0;i < total;i++){
        color[i](0) = (double) rand() / RAND_MAX;
        color[i](1) = (double) rand() / RAND_MAX;
        color[i](2) = (double) rand() / RAND_MAX;
        color[i].normalize();
    }
    Eigen::MatrixXd C(faceNum,3);
    
    //compute the components class, the components are not connected yet
    heightField(total);
    
    GeneralGraph_DArraySArraySpatVarying(faceNum,total);
        
//    getHalfEdge(V.rows(), faceNum, result, &ver, &e, &triFace);
    
//    build(&ver, totalVer);
    
//    integrate(&triFace, result, faceNum);
    
//    integrate(&triFace, result, faceNum);
    
//    integrate(&triFace, result, faceNum);
    
    
    for(int i = 0 ;i < faceNum;i++){
//        C.row(i) = color[labels[i].visitItem(0)];
        C.row(i) = color[result[i]];
    }
    
    
    
    // Plot the mesh
    igl::viewer::Viewer viewer;
    viewer.core.show_lines = false;
    viewer.data.set_mesh(V, F);
    viewer.data.set_normals(N_faces);
    viewer.data.set_colors(C);
    
    viewer.launch();
    
    for(int i = 0;i < faceNum;i++){
        components[i].~set();
        edges[i].~set();
        labels[i].~set();
    }
    operator delete(rawMemory);
    operator delete(rawMemory1);
    operator delete(rawMemory3);
    delete[] color;
    delete[] d;
    free(result);
    free(ver);
    free(e);
    free(triFace);
    free(ccp);
    
    return 0;
}

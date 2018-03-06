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
connectedComponents* ccp;
int NfaceRows;
int faceNum;

using namespace std;

int main(int argc, char *argv[])
{
    int total, labelCost;
    total = 600; labelCost = 0.5*10;
    d = new Eigen::Vector3d[total];

    for(int i = 0;i < total/2;i++){
        d[i](0) = Guass();
        d[i](1) = Guass();
        d[i](2) = Guass();
        d[i].normalize();
        d[i+total/2] = -d[i];
    }
    
    // Load a mesh in OFF format
    igl::readOFF( "/Users/Lagrant/libigl/tutorial/shared/bunny.off", V, F);
    
    // Compute per-face normals
    igl::per_face_normals(V,F,N_faces);
    NfaceRows = N_faces.rows();
    faceNum = F.rows();
    
    // To apply the halfedge data structure
    int totalVer = faceNum*3;
    ver = (vertex*) malloc(sizeof(vertex)*totalVer);
    e = (edge*) malloc(sizeof(edge)*totalVer);
    triFace = (face*) malloc(sizeof(face)*faceNum);
    ccp = (connectedComponents*) malloc(sizeof(connectedComponents)*faceNum*total);
    
    Eigen::RowVector3d* color = new Eigen::RowVector3d[total];
    result = (int*) malloc(sizeof(int)*faceNum);
    memset(result, -1, sizeof(int)*faceNum);
    
    //initialize the set class

    void* rawMemory = operator new(NfaceRows*sizeof(set<int>));
    
    labels = reinterpret_cast<set<int>*>(rawMemory);
    
    for(int i = 0;i < faceNum;i++){
//NfaceRows = faceNum, NfaceRows means the total numebr of norms of faces while faceNum means the total  number of faces
        new (&labels[i])set<int>(total);
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
    getHalfEdge(&ver, &e, &triFace);
    
    build(&ver, totalVer);
    
    heightField(total);
    
    GeneralGraph_DArraySArraySpatVarying(faceNum, total, labelCost);
    
    
    integrate();
    
    integrate();
    
//    integrate();
    
    
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

        labels[i].~set();
    }

    operator delete(rawMemory);
    delete[] color;
    delete[] d;
    free(result);
    free(ver);
    free(e);
    free(triFace);
    free(ccp);
    
    return 0;
}

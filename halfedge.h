//
//  halfedge.h
//  example
//
//  Created by 张鑫禹 on 2018/1/29.
//
//

#ifndef halfedge_h
#define halfedge_h
#include <stdio.h>
#include <Eigen/LU>
#include "set.h"
struct edge;

struct vertex {
    Eigen::Matrix<double,3,1> vert;
    
    edge* startEdge;
};

struct face {
    edge* adjacentEdge;
    int numbering;
    int label = -1;
};

struct edge {
    vertex* endVertex;
    edge* pair = NULL;
    edge* next;
    face* face;
};

extern Eigen::MatrixXd V;
extern Eigen::MatrixXi F;
extern Eigen::MatrixXd N_faces;
extern Eigen::Vector3d* d;
extern set<int>* components;
extern set<vector<int>>* edges;
extern set<int>* labels;
extern int* result;



inline bool equal(vertex* vi, vertex* vj){
    if(vi->vert[0] == vj->vert[0] && vi->vert[2] == vj->vert[2] && vi->vert[1] == vj->vert[1])
        return true;
    else return false;
}

inline void assignEdge(vertex** ver, int end, edge** e, int i){
    
    (*e)[3*i].endVertex = &(*ver)[end-1];
    (*e)[3*i].next = &(*e)[3*i+1];
    (*ver)[end-2].startEdge = &(*e)[3*i];
    
    (*e)[3*i+1].endVertex = &(*ver)[end];
    (*e)[3*i+1].next = &(*e)[3*i+2];
    (*ver)[end-1].startEdge = &(*e)[3*i+1];
    
    (*e)[3*i+2].endVertex = &(*ver)[end-2];
    (*e)[3*i+2].next = &(*e)[3*i];
    (*ver)[end].startEdge = &(*e)[3*i+2];
}


void getHalfEdge(int verNum, int faceNum, vertex** ver,  edge** e, face** triFace);

void build(vertex** ver, int totalVer);

#endif /* halfedge_h */

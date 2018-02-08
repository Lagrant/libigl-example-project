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
struct edge;

struct vertex {
    double x;
    double y;
    double z;
    
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
bool equal(vertex* vi, vertex* vj);

void assignEdge(vertex** ver, int end, edge** e, int i);

void getHalfEdge(int verNum, int faceNum, int* result, vertex** ver,  edge** e, face** triFace);

void build(vertex** ver, int totalVer);

//void adjacentFace(face* f);

#endif /* halfedge_h */

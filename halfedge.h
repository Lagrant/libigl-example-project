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

void getHalfEdge(int verNum, int faceNum, vertex** ver,  edge** e, face** triFace);

void build(vertex** ver, int totalVer);

//void adjacentFace(face* f);

#endif /* halfedge_h */

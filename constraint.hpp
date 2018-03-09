//
//  constraint.hpp
//  example
//
//  Created by 张鑫禹 on 2018/3/5.
//
//

#ifndef constraint_hpp
#define constraint_hpp

#include <stdio.h>
#include <vector>
#include <Eigen/LU>
#include <math.h>
#include <stdlib.h>
#include <queue>
#include <igl/viewer/Viewer.h>
#include "GCoptimization.h"
#include "halfedge.h"

#define Random(x) rand()%x

typedef struct {
    Eigen::Vector3d direction;
    vector<int> conFaces;
} connectedComponents;

extern Eigen::MatrixXd V;
extern Eigen::MatrixXi F;
extern Eigen::MatrixXd N_faces;
extern Eigen::Vector3d* d;
extern set<int>* labels;
extern int* result;
extern face* triFace;
extern int total;
extern int faceNum;

double Guass();

bool interEmpty(Eigen::Matrix3d P, Eigen::Matrix<double,1,3> j0, Eigen::Matrix<double,1,3> j1, Eigen::Matrix<double,1,3> k0, Eigen::Matrix<double,1,3> k1);

bool overlap(Eigen::Matrix3d P, connectedComponents L, int label);

int seekTheSameLabel(int label, int face);

bool interLock(edge* e, Eigen::Vector3d d);

int heightField();

void GeneralGraph_DArraySArraySpatVarying(int label_cost);

void getNeighborLabels(int dLabel[], const face* f);

inline int compare(int dLabel[]){
    if(dLabel[1] == dLabel[2] || dLabel[1] == dLabel[0]){
        return dLabel[1];
    } else if(dLabel[2] == dLabel[0])
        return dLabel[2];
    else return -1;
}

void integrate();

void merge();

int append(igl::viewer::Viewer& viewer, vector<int>& vec, int k, Eigen::RowVector3d M, Eigen::RowVector3d color);

#endif /* constraint_hpp */

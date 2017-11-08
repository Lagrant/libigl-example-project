//
//  constraints.h
//  example
//
//  Created by 张鑫禹 on 2017/11/3.
//
//

#ifndef constraints_h
#define constraints_h
#include <vector>
#include "set.h"
#include <Eigen/LU>

void sort(int* q,int len){
    int temp;
    for(int i = 0;i < len;i++){
        for(int j = i+1;j < len;j++){
            if(q[i] < q[j]){
                temp = q[i];
                q[i] = q[j];
                q[j] = temp;
            }
        }
    }
}

bool interEmpty(Eigen::Matrix<double,1,3> j0, Eigen::Matrix<double,1,3> j1, Eigen::Matrix<double,1,3> k0, Eigen::Matrix<double,1,3> k1){
    Eigen::Matrix<double,1,3> vec1 = j0 - j1, vec2 = k0 - k1, vec3 = k0-j0, vec4 = k1 - j0;
    vec1.normalize(); vec2.normalize(); vec3.normalize(); vec4.normalize();
    if(vec1 == vec2)
        return true;
    Eigen::Matrix3d space;
    space.row(0) = vec1;
    space.row(1) = vec3;
    space.row(2) = vec4;
    Eigen::Matrix3d inverseSpace;
    bool invertible;
    double determinat;
    space.computeInverseAndDetWithCheck(inverseSpace,determinat,invertible);
    return invertible ? true : false;
 }

void heightField(Eigen::MatrixXd N_faces, Eigen::Vector3d* d, set<int>* components, set<vector<int>>* edges, Eigen::MatrixXi F, int total){
    int NfaceRows = N_faces.rows();
    for(int j = 0;j < total;j++){
        for(int i = 0;i < NfaceRows;i++){
            if(N_faces.row(i).dot(d[j]) >= 0){
                components[j].addItem(i);
                Eigen::Matrix<int,1,3> e = F.row(i);
                int q[] = {e(0),e(1),e(2)};
                sort(q,3);
                vector<int> k(2);
                k.at(0) = q[0]; k.at(1) = q[1];
                edges[j].addItem(k);
                k.at(0) = q[0]; k.at(1) = q[2];
                edges[j].addItem(k);
                k.at(0) = q[1]; k.at(1) = q[2];
                edges[j].addItem(k);
            }
        }
    }
}

void overlap(Eigen::Vector3d d, set<int> componets, set<vector<int>> edges, Eigen::MatrixXd V, Eigen::MatrixXi F, int total){
    Eigen::Matrix<double,1,3> dT = d.transpose();
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();  //the unit matrix
    Eigen::Matrix3d P = I - d*dT;
    for(int i = 0;i < total;i++){
        for(int j = 0;j < total;j++){
            vector<int> e1 = edges.visitItem(i);
            vector<int> e2 = edges.visitItem(j);
            if(e1.empty() || e2.empty())
                cerr<<"the index is out of boundary in overlap";
            if(e1 == e2)
                continue;
            Eigen::Matrix<double,1,3> j0 = V.row(e1.at(0));
            Eigen::Matrix<double,1,3> j1 = V.row(e1.at(1));
            Eigen::Matrix<double,1,3> k0 = V.row(e2.at(0));
            Eigen::Matrix<double,1,3> k1 = V.row(e2.at(1));
            bool isEmpty = interEmpty(j0,j1,k0,k1);
        }
    }
}

#endif /* constraints_h */

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

bool interEmpty(Eigen::Matrix3d P, Eigen::Matrix<double,1,3> j0, Eigen::Matrix<double,1,3> j1, Eigen::Matrix<double,1,3> k0, Eigen::Matrix<double,1,3> k1){
    Eigen::Matrix<double,3,1> j0T = j0.transpose();
    Eigen::Matrix<double,3,1> j1T =j1.transpose();
    Eigen::Matrix<double,3,1> k0T =k0.transpose();
    Eigen::Matrix<double,3,1> k1T = k1.transpose();
    j0T = P*j0T;
    j1T = P*j1T;
    k0T = P*k0T;
    k1T = P*k1T;
    Eigen::Matrix<double,3,1> vec1 = j0T-j1T, vec2 = k0T-k1T, vec3 = k0T-j0T, vec4 = k1T-j0T;
    vec1.normalize(); vec2.normalize(); vec3.normalize(); vec4.normalize();
    if(vec1 == vec2)
        return true;
    Eigen::Matrix3d space;
    space.col(0) = vec1;
    space.col(1) = vec3;
    space.col(2) = vec4;
    Eigen::Matrix3d inverseSpace;
    bool invertible;
    double determinat;
    space.computeInverseAndDetWithCheck(inverseSpace,determinat,invertible);
    return invertible ? true : false;
 }

void appendSeams(set<vector<int>>* seams, const int pos,vector<int> k){
    int scale = seams[pos].scale();
    for(int i = 0;i < scale;i++){
        vector<int> item = seams[pos].visitItem(i);
        if(k.at(0)==item.at(0) && k.at(1)==item.at(1)){
            seams[pos].removeItem(k);
            return;
        }
    }
    seams[pos].addItem(k);
}

void removeSeams(Eigen::Matrix<int,1,3> e, set<vector<int>>* seams, const int current, const int pos, const int faceIndex){
    int q[] = {e(0),e(1),e(2)};
    sort(q,3);
    vector<int> k(3);
    k.at(0) = q[0]; k.at(1) = q[1]; k.at(2) = faceIndex;
    appendSeams(seams, current, k);
    appendSeams(seams, pos, k);
    k.at(0) = q[0]; k.at(1) = q[2]; k.at(2) = faceIndex;
    appendSeams(seams, current, k);
    appendSeams(seams, pos, k);
    k.at(0) = q[1]; k.at(1) = q[2]; k.at(2) = faceIndex;
    appendSeams(seams, current, k);
    appendSeams(seams, pos, k);
}

void heightField(Eigen::MatrixXd N_faces, Eigen::Vector3d* d, set<int>* components, set<vector<int>>* edges, set<vector<int>>* seams, Eigen::MatrixXi F, int total){
    int NfaceRows = N_faces.rows();
    int* cover = new int[NfaceRows];
    memset(cover,1,sizeof(cover));
    for(int j = 0;j < total;j++){
        for(int i = 0;i < NfaceRows;i++){
            if(N_faces.row(i).dot(d[j])>=0 && cover[i]){
                cover[i] = 0;
                components[j].addItem(i);
                Eigen::Matrix<int,1,3> e = F.row(i);
                int q[] = {e(0),e(1),e(2)};
                sort(q,3);
                vector<int> k(3);
                k.at(0) = q[0]; k.at(1) = q[1]; k.at(2) = i;
                edges[j].addItem(k);
                appendSeams(seams, j, k);
                k.at(0) = q[0]; k.at(1) = q[2]; k.at(2) = i;
                edges[j].addItem(k);
                appendSeams(seams, j, k);
                k.at(0) = q[1]; k.at(1) = q[2]; k.at(2) = i;
                edges[j].addItem(k);
                appendSeams(seams, j, k);
            }
        }
    }
    delete[] cover;
}

void overlap(const Eigen::Vector3d d, set<vector<int>>* edges, set<int>* components, set<vector<int>>* seams, const Eigen::MatrixXd V, Eigen::MatrixXi F, int total, int& pos, const int current){
    Eigen::Matrix<double,1,3> dT = d.transpose();
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();  //the unit matrix
    Eigen::Matrix3d P = I - d*dT;
    for(int i = 0;i < total;i++){
        vector<int> e1 = edges[current].visitItem(i);
        Eigen::Matrix<double,1,3> j0 = V.row(e1.at(0));
        Eigen::Matrix<double,1,3> j1 = V.row(e1.at(1));
        for(int j = i+1;j < total;j++){
            vector<int> e2 = edges[current].visitItem(j);
            if(e1.empty() || e2.empty()){
                cerr<<"the index is out of boundary in overlap";
                return;
            }
            if(e1 == e2)
                continue;
            
            Eigen::Matrix<double,1,3> k0 = V.row(e2.at(0));
            Eigen::Matrix<double,1,3> k1 = V.row(e2.at(1));
            bool isEmpty = interEmpty(P,j0,j1,k0,k1);
            if(isEmpty)
                continue;
            
            int faceIndex = e2.at(2);
            Eigen::Matrix<int,1,3> e = F.row(faceIndex);
            removeSeams(e, seams, current, pos, faceIndex);
            edges[current].removeItem(e2);
            edges[pos].addItem(e2);
            components[current].removeItem(faceIndex);
            components[pos].addItem(faceIndex);
        }
    }
    if(!components[pos].isEmpty())
        pos++;
}

/*
 for every direction d and the components related to d, calculate the interlcok individually.
 */
void interlock(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N_faces, const Eigen::Vector3d d, set<vector<int>>* edges, set<int>* components, set<vector<int>>* seams, int total, int& pos, const int current){
    for(int i = 0;i < total;i++){
        vector<int> sEdges = seams[current].visitItem(i);
        int s = sEdges.at(2);
        Eigen::Matrix<double,1,3> v1 = V.row(sEdges.at(0)),
        v2 = V.row(sEdges.at(1)),
        v0 = v1-v2,
        n = N_faces.row(s),
        nt = n.cross(v0);
        
        nt.normalize();
        if(nt*d >= 0)
            continue;
        
        components[current].removeItem(s);
        components[pos].addItem(s);
        edges[current].removeItem(sEdges);
        edges[pos].addItem(sEdges);
        
        Eigen::Matrix<int,1,3> e = F.row(s);
        removeSeams(e, seams, current, pos, s);
    }
}
/*
 calculate the value total in the later coding
 construct an array to store the directions
 */
#endif /* constraints_h */

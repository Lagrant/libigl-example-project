//
//  halfedge.cpp
//  example
//
//  Created by 张鑫禹 on 2018/1/29.
//
//

#include "halfedge.h"
#include <Eigen/LU>
#include <vector>
#include <queue>
#include <iostream>
#include <assert.h>
#include "set.h"

extern Eigen::MatrixXd V;
extern Eigen::MatrixXi F;
extern Eigen::MatrixXd N_faces;
extern Eigen::Vector3d* d;
extern set<int>* components;
extern set<vector<int>>* edges;
extern set<int>* labels;
extern int* result;



inline bool equal(vertex* vi, vertex* vj){
    if(vi->x == vj->x && vi->y == vj->y && vi->z == vj->z)
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

void getHalfEdge(int verNum, int faceNum, int* result, vertex** ver, edge** e, face** triFace){
    int f[3], end = 0;
    for(int i = 0; i < faceNum; i++){
        (*triFace)[i].adjacentEdge = &(*e)[3*i];
        (*triFace)[i].label = result[i];
        (*triFace)[i].numbering = i;
        (*e)[3*i+1].face = &(*triFace)[i];
        (*e)[3*i+2].face = &(*triFace)[i];
        (*e)[3*i].face = &(*triFace)[i];
        for(int j = 0; j< 3; j++){
            f[j] = F(i,j);
            assert(f[j] <= verNum);
            (*ver)[end].x = V(f[j],0);
            (*ver)[end].y = V(f[j],1);
            (*ver)[end].z = V(f[j],2);
            end++;
        }
        
        assignEdge(ver, end-1, e, i);
    }
}

void build(vertex** ver, int totalVer) {
    int* isPaired = (int*)malloc(sizeof(int)*totalVer); // new int[totalVer];
    memset(isPaired, 0, sizeof(int)*totalVer);
    
    for(int i = 0; i < totalVer; i++){
        if(isPaired[i])
            continue;
        
        for(int j = 0; j < totalVer; j++){
            if(isPaired[j])
                continue;
            
            if(equal((*ver)[i].startEdge->endVertex, &(*ver)[j]) && equal((*ver)[j].startEdge->endVertex, &(*ver)[i])){
                (*ver)[i].startEdge->pair = (*ver)[j].startEdge;
                (*ver)[j].startEdge->pair = (*ver)[i].startEdge;
                isPaired[j] = 1;
                isPaired[i] = 1;
                break;
            }
        }
        
    }
    //    delete [] isPaired;
    free(isPaired);
}

//void adjacentFace(face* f) {
//    std::queue<face*> adjFaces;
//    adjFaces.push(f);
//    do{
//        f = adjFaces.front();
//        adjFaces.pop();
//        f->visited = 1;
//        std::cout<<f<<"\n";
//        vertex v = (*f->adjacentEdge->endVertex);
//        std::cout<<"x = "<<v.x<<" y = "<<v.y<<" z = "<<v.z<<std::endl;
//        edge* e = f->adjacentEdge;
//        edge* eRun = f->adjacentEdge;
//        do{
//            if(eRun->pair != NULL && !(eRun->pair->face->visited))
//                adjFaces.push(eRun->pair->face);
//            eRun = eRun->next;
//        } while (e != eRun);
//    } while (!adjFaces.empty());
//}

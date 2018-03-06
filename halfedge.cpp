//
//  halfedge.cpp
//  example
//
//  Created by 张鑫禹 on 2018/1/29.
//
//

#include "halfedge.h"


void getHalfEdge(vertex** ver, edge** e, face** triFace){
    int f[3], end = 0;
    for(int i = 0; i < faceNum; i++){
        (*triFace)[i].adjacentEdge = &(*e)[3*i];
        (*triFace)[i].numbering = i;
        (*e)[3*i+1].face = &(*triFace)[i];
        (*e)[3*i+2].face = &(*triFace)[i];
        (*e)[3*i].face = &(*triFace)[i];
        for(int j = 0; j< 3; j++){
            f[j] = F(i,j);
            (*ver)[end].vert = V.row(f[j]);
            end++;
        }
        
        assignEdge(ver, end-1, e, i);
    }
}

void build(vertex** ver, int totalVer) {
    int* isPaired = (int*)malloc(sizeof(int)*totalVer);
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
    
    free(isPaired);
}

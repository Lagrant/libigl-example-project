//
//  constraints.h
//  example
//
//  Created by 张鑫禹 on 2017/11/3.
//
//

#ifndef constraints_h
#define constraints_h
#include "set.h"

void heightField(Eigen::MatrixXd N_faces, Eigen::Vector3d* d, set* components, int total){
    int NfaceRows = N_faces.rows();
    for(int j = 0;j < total;j++){
        for(int i = 0;i < NfaceRows;i++){
            if(N_faces.row(i).dot(d[j]) >= 0){
                components[j].addItem(i);
            }
        }
    }
}

void overlap(){
    
}

#endif /* constraints_h */

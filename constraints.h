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
#include <Eigen/LU>
#include <math.h>
#include <stdlib.h>
#include <queue>
#include "GCoptimization.h"
#include "set.h"
#include "halfedge.h"

#define Random(x) rand()%x
#define MAX 1000

typedef struct {
    Eigen::Vector3d direction;
    vector<int> conFaces;
} connectedComponents;

extern Eigen::MatrixXd V;
extern Eigen::MatrixXi F;
extern Eigen::MatrixXd N_faces;
extern Eigen::Vector3d* d;
extern set<int>* components;
extern set<vector<int>>* edges;
extern set<int>* labels;
extern int* result;
extern face* triFace;
extern connectedComponents* ccp;

double Guass(){
    double u;
    double v;
    double r;
    do{
        v = ((double) rand() / (RAND_MAX)) * 2 - 1;
        u = ((double) rand() / (RAND_MAX)) * 2 - 1;
        r = u * u + v * v;
    }while(r >= 1 || r == 0);
    
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}

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
    if(vec1 == vec2){
        return true;
    }
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

double match(int k[], int l[], int total){
    int c[3] = {-1,-1,-1}, cur = 0;
    for(int i = 0;i < total;i++){
        for(int j = 0;j < total;j++){
            if(k[i] == l[j])
                c[cur++] = k[i];
        }
    }
    if(c[1] == -1)
        return -1;
    Eigen::Matrix<double,3,1> temp1 = V.row(c[0]);
    Eigen::Matrix<double,3,1> temp2 = V.row(c[1]);
    temp1 -= temp2;
    return temp1.norm();
}

int heightField(const int total){
    int NfaceRows = N_faces.rows(), end = 0;
    int *mark = (int*) malloc(sizeof(int)*NfaceRows);
    
    for(int j = 0;j < total;j++){
        memset(mark,0,sizeof(int)*NfaceRows);

        for(int i = 0;i < NfaceRows;i++){
            if(mark[i])
                continue;
            if(N_faces.row(i)*d[j] >= 0){
                mark[i] = 1;
                ccp[end].direction = d[j];
                components[j].addItem(end);
                queue<face*> q;
                face* f;
                q.push(&triFace[i]);
                do{
                    f = q.front();
                    q.pop();
                    ccp[end].conFaces.push_back(f->numbering);
                    mark[f->numbering] = 1;
                    labels[f->numbering].addItem(j);
                    
                    edge* e = f->adjacentEdge;
                    f = e->pair->face;
                    if(N_faces.row(f->numbering)*d[j] >= 0 && !mark[f->numbering])
                        q.push(f);
                    
                    e = e->next;
                    f = e->pair->face;
                    if(N_faces.row(f->numbering)*d[j] >= 0 && !mark[f->numbering])
                        q.push(f);
                    
                    e = e->next;
                    f = e->pair->face;
                    if(N_faces.row(f->numbering)*d[j] >= 0 && !mark[f->numbering])
                        q.push(f);
                    
                } while(!q.empty());
                end++;
            }
        }
    }
    free(mark);
    return end;
}

bool overlap(Eigen::Matrix3d P, Eigen::Vector3d d, connectedComponents L, int label){
    
    int l0 = F(label,0), l1 = F(label,1), l2 = F(label,2);
    Eigen::Matrix<double,1,3> v0 = V.row(l0), v1 = V.row(l1), v2 = V.row(l2);
    int total = L.conFaces.size();
    for(int j = 0;j < total;j++){
        int item = L.conFaces.at(j);
        int k0 = F(item,0), k1 = F(item,1), k2 = F(item,2);
        Eigen::Matrix<double,1,3> w0 = V.row(k0), w1 = V.row(k1), w2 = V.row(k2);
        //cout<<"direction  = ("<<d.row(0)<<", "<<d.row(1)<<", "<<d.row(2)<<")\n";
        if(!interEmpty(P,v2,v1,w0,w2))
            return false;
        if(!interEmpty(P,v0,v1,w0,w1) || !interEmpty(P,v0,v1,w0,w2) || !interEmpty(P,v0,v1,w2,w1) || !interEmpty(P,v0,v2,w0,w1) || !interEmpty(P,v0,v2,w0,w2) || !interEmpty(P,v0,v2,w2,w1) || !interEmpty(P,v2,v1,w0,w1) || !interEmpty(P,v2,v1,w2,w1))
            return true;
       
    }
    return false;
}

void GeneralGraph_DArraySArraySpatVarying(int num_pixels,int num_labels)
{
//    int *result = new int[num_pixels];   // stores result of optimization
    
    // first set up the array for data costs
    int *data = (int*) malloc(sizeof(int)*num_pixels*num_labels);
    int *smooth = (int*) malloc(sizeof(int)*num_labels*num_labels);
    
    for ( int i = 0; i < num_pixels; i++ ){
        for (int l = 0; l < num_labels; l++ ){
            if(labels[i].isExist(l) != -1){
                data[i*num_labels+l] = 0;
            }
            else data[i*num_labels+l] = MAX;
        }
    }
    
    // next set up the array for smooth costs
    for ( int l1 = 0; l1 < num_labels; l1++ )
        for (int l2 = 0; l2 < num_labels; l2++ ){
            (l1 == l2)? smooth[l1+l2*num_labels] = 0: smooth[l1+l2*num_labels] = 1;
        }
    
    try{
        GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_pixels,num_labels);
        gc->setDataCost(data);
        gc->setSmoothCost(smooth);
        
        for(int i = 0;i < num_pixels;i++){
            int scale = labels[i].scale();
            int li = labels[i].visitItem(Random(scale));
            Eigen::Matrix<double,1,3> dT = d[li].transpose();
            Eigen::Matrix3d I = Eigen::Matrix3d::Identity();  //the unit matrix
            Eigen::Matrix3d P0 = I - d[li]*dT;
            
            gc->setNeighbors(i,i,0);
            for(int j = i+1;j < num_pixels;j++){
                int scale = labels[j].scale();
                int lj = labels[j].visitItem(Random(scale));   //fix a bug that changing label[k] to label[j]
                Eigen::Matrix<double,1,3> dT = d[lj].transpose();
                Eigen::Matrix3d I = Eigen::Matrix3d::Identity();  //the unit matrix
                Eigen::Matrix3d P = I - d[lj]*dT;
                
                if(li == lj){
                    //cout<<"label = "<<li<<"\n";
                    gc->setNeighbors(i,j,0);
                    gc->setNeighbors(j,i,0);
                }
                else if(N_faces.row(i)*d[li]<0 || N_faces.row(j)*d[lj]<0){
                    gc->setNeighbors(i,j,MAX);
                    gc->setNeighbors(j,i,MAX);
                }
                else if(overlap(P0,d[li],components[li],li) || overlap(P,d[lj],components[lj],lj)){
                    gc->setNeighbors(i,j,MAX);
                    gc->setNeighbors(j,i,MAX);
                }
                else{
                    int p[3],q[3];
                    for(int i = 0;i < 3;i++){
                        p[i] = F(li,i);
                        q[i] = F(lj,i);
                    }
                    double value = match(p,q,3);
                    if(value < 0){
                        gc->setNeighbors(i,j,0);
                        gc->setNeighbors(j,i,0);
                    }
                    else{
                        gc->setNeighbors(i,j,value);
                        gc->setNeighbors(j,i,value);
                    }
                }
            }
        }
        
        printf("\nBefore optimization energy is %lld",gc->compute_energy());
        gc->expansion();// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
        printf("\nAfter optimization energy is %lld",gc->compute_energy());
        
        for ( int  i = 0; i < num_pixels; i++ )
            result[i] = gc->whatLabel(i);
        
        delete gc;
    } catch (GCException e){
        e.Report();
    }
    
//    delete [] result;
    free(data);
    free(smooth);
    
}

void getNeighborLabels(int dLabel[], face* f) {
    int count = 0;
    edge* e = f -> adjacentEdge;
    edge* eRun = f ->adjacentEdge;
    
    do {
        if(eRun -> pair != NULL){
            face* pairFace = eRun -> pair -> face;
            dLabel[count] = pairFace ->label;
            count++;
        }
        eRun = eRun -> next;
    } while(e != eRun);
}

inline int compare(int dLabel[]){
    if(dLabel[1] == dLabel[2] || dLabel[1] == dLabel[0]){
        return dLabel[1];
    } else if(dLabel[2] == dLabel[0])
        return dLabel[2];
    else return -1;
}

void integrate(face** triFace, int* result, int faceNum) {
//    std::queue<face*> q;
//    q.push((*triFace));
    for(int i = 0; i < faceNum; i++){
        
        int dLabel[3];
        memset(dLabel, -1, sizeof(int)*3);
        
        getNeighborLabels(dLabel, &(*triFace)[i]);
        
        int dl = compare(dLabel);
        if(dl != -1 && dl != result[i] /*&& N_faces.row(i)*d[dl] >= 0*/){
//            std::cout<<"before result["<<i<<"] = "<<result[i]<<"\t";
            result[i] = dl;
            (*triFace)[i].label = dl;
//            std::cout<<"after result["<<i<<"] = "<<result[i]<<"\n";
        }
    }
}

/*
 for every direction d and the components related to d, calculate the interlock individually.
 */
/*void interlock(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N_faces, const Eigen::Vector3d d, set<vector<int>>* edges, set<int>* components, set<vector<int>>* seams, int total, int& pos, const int current){
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
    if(!components[pos].isEmpty())
       pos++;
}*/
/*
 calculate the value total in the later coding
 construct an array to store the directions
 */
#endif /* constraints_h */

//
//  constraint.cpp
//  example
//
//  Created by 张鑫禹 on 2018/3/5.
//
//

#include "constraint.hpp"

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

bool overlap(Eigen::Matrix3d P, connectedComponents L, int label){

    int l0 = F(label,0), l1 = F(label,1), l2 = F(label,2);
    Eigen::Matrix<double,1,3> v0 = V.row(l0), v1 = V.row(l1), v2 = V.row(l2);
    int total = L.conFaces.size();
    for(int j = 0;j < total;j++){
        int item = L.conFaces.at(j);
        int k0 = F(item,0), k1 = F(item,1), k2 = F(item,2);
        Eigen::Matrix<double,1,3> w0 = V.row(k0), w1 = V.row(k1), w2 = V.row(k2);

        if(!interEmpty(P,v2,v1,w0,w2))
            return false;
        if(!interEmpty(P,v0,v1,w0,w1) || !interEmpty(P,v0,v1,w0,w2) || !interEmpty(P,v0,v1,w2,w1) || !interEmpty(P,v0,v2,w0,w1) || !interEmpty(P,v0,v2,w0,w2) || !interEmpty(P,v0,v2,w2,w1) || !interEmpty(P,v2,v1,w0,w1) || !interEmpty(P,v2,v1,w2,w1))
            return true;

    }
    return false;
}

int seekTheSameLabel(int label, int face){
    int scale = labels[face].scale();
    for(int i = 0; i < scale; i++){
        if(labels[face].visitItem(i) == label)
            return 1;
    }
    return 0;
}


int heightField(const int total){
    int NfaceRows = N_faces.rows(), end = 0;
    int *mark = (int*) malloc(sizeof(int)*NfaceRows);

    for(int j = 0;j < total;j++){
        memset(mark,0,sizeof(int)*NfaceRows);

        for(int i = 0; i < NfaceRows; i++){
            if(mark[i])
                continue;
            if(N_faces.row(i)*d[j] >= 0){
                mark[i] = 1;
                ccp[end].direction = d[j];
                //                components[j].addItem(end);

                queue<face*> q;
                face* f;
                edge* e;
                q.push(&triFace[i]);

                do{
                    f = q.front();
                    q.pop();
                    ccp[end].conFaces.push_back(f->numbering);
                    labels[f->numbering].addItem(j);

                    e = f->adjacentEdge;
                    f = e->pair->face;
                    if(N_faces.row(f->numbering)*d[j] >= 0 && !mark[f->numbering]){
                        q.push(f);
                        mark[f->numbering] = 1;
                    }

                    e = e->next;
                    f = e->pair->face;
                    if(N_faces.row(f->numbering)*d[j] >= 0 && !mark[f->numbering]){
                        q.push(f);
                        mark[f->numbering] = 1;
                    }

                    e = e->next;
                    f = e->pair->face;
                    if(N_faces.row(f->numbering)*d[j] >= 0 && !mark[f->numbering]){
                        q.push(f);
                        mark[f->numbering] = 1;
                    }

                } while(!q.empty());
                end++;
            }
        }
    }
    free(mark);
    return end;
}

void GeneralGraph_DArraySArraySpatVarying(int num_pixels,int num_labels, int label_cost){

    // first set up the array for data costs
    int *data = (int*) malloc(sizeof(int)*num_pixels*num_labels);
    int *smooth = (int*) malloc(sizeof(int)*num_labels*num_labels);

    for ( int i = 0; i < num_pixels; i++ ){
        for (int l = 0; l < num_labels; l++ ){
            if(labels[i].isExist(l) != -1){
                data[i*num_labels+l] = 0;
            }
            else data[i*num_labels+l] = GCO_MAX_ENERGYTERM;
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
        gc->setLabelCost(label_cost);

        int* mark = (int*) malloc(sizeof(int)*num_pixels);
        memset(mark,0,sizeof(int)*num_pixels);
        face* f, *adjF;
        edge* e;
        queue<face*> q;
        q.push(&triFace[0]);
        mark[0] = 1;
        result[0] = labels[0].visitItem(Random(labels[0].scale()));

        do {
            f = q.front();
            q.pop();
            e = f->adjacentEdge;
            int i = f->numbering, j, li;
            edge *eRun = e;
            if(result[i] == -1)
                result[i] = labels[i].visitItem(Random(labels[i].scale()));

            //select a random initial label to be the seed. Traverse the adjacent faces. If the same label is found, set the label on the adjacent face. Traverse and seek the same label if the label is set, or set a random inital label and then seek the same label for each face on the front of the queue. Remove overlap, consider seam conflicts only.
            do {
                adjF = eRun->pair->face;
                j = adjF->numbering;

                if(!mark[j]){
                    mark[j] = 1;
                    q.push(&triFace[j]);
                }

                if(result[j] == -1){
                    li = seekTheSameLabel(result[i], j);

                    if(!li){
                        //                        if(interLock(eRun, d[result[i]])){
                        //                            Eigen::Matrix<double,3,1> value = eRun->endVertex->vert - eRun->pair->endVertex->vert;
                        //                            gc->setNeighbors(i,j,value.norm());
                        //                        } else gc->setNeighbors(i, j, GCO_MAX_ENERGYTERM);
                        Eigen::Matrix<double,3,1> value = eRun->endVertex->vert - eRun->pair->endVertex->vert;
                        gc->setNeighbors(i,j,value.norm());
                    } else {
                        result[j] = result[i];
                        gc->setNeighbors(i, j, 0);
                    }
                } else if(result[i] == result[j])
                    gc->setNeighbors(i, j, 0);
                else {
                    //                    if(interLock(eRun, d[result[i]])){
                    //                        Eigen::Matrix<double,3,1> value = eRun->endVertex->vert - eRun->pair->endVertex->vert;
                    //                        gc->setNeighbors(i,j,value.norm());
                    //                    } else gc->setNeighbors(i, j, GCO_MAX_ENERGYTERM);
                    Eigen::Matrix<double,3,1> value = eRun->endVertex->vert - eRun->pair->endVertex->vert;
                    gc->setNeighbors(i,j,value.norm());
                }

                eRun = eRun->next;

            } while(eRun != e);

        } while(!q.empty());


        printf("\nBefore optimization energy is %lld",gc->compute_energy());
        gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
        printf("\nAfter optimization energy is %lld",gc->compute_energy());

        for ( int  i = 0; i < num_pixels; i++ ){
            result[i] = gc->whatLabel(i);
            triFace[i].label = result[i];
        }

        free(mark);
        delete gc;
    } catch (GCException e){
        e.Report();
    }

    free(data);
    free(smooth);

}

void getNeighborLabels(int dLabel[], const face* f){
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

void integrate(){
    
    for(int i = 0; i < faceNum; i++){

        int dLabel[3];
        memset(dLabel, -1, sizeof(int)*3);

        getNeighborLabels(dLabel, &triFace[i]);

        int dl = compare(dLabel);
        if(dl != -1 && dl != result[i]){
            result[i] = dl;
            triFace[i].label = dl;
        }
    }
}

void merge(){
    
}

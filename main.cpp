#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_face_normals.h>
#include <iostream>
#include <vector>
#include "set.h"
#include "constraints.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd N_faces;
using namespace std;

int main(int argc, char *argv[])
{
    int total;
    double x,y,z;
    cin>>total;
    Eigen::Vector3d* d =new Eigen::Vector3d[total];
    for(int i = 0;i < total;i++){
        cin>>x>>y>>z;
        d[i]<<x,y,z;
    }
    
    // Load a mesh in OFF format
    igl::readOFF( "/Users/Lagrant/libigl/tutorial/shared/fandisk.off", V, F);
    
    // Compute per-face normals
    igl::per_face_normals(V,F,N_faces);
    int NfaceRows = N_faces.rows();
    int faceNum = F.rows();
    
    //initialize the set class
    void* rawMemory = operator new(total*sizeof(set<int>));
    void* rawMemory1 = operator new(3*total*sizeof(set<vector<int>>));
    set<int>* components = reinterpret_cast<set<int>*>(rawMemory);
    set<vector<int>>* edges = reinterpret_cast<set<vector<int>>*>(rawMemory1);
    for(int i = 0;i < total;i++){
        new (&components[i])set<int>(faceNum);
        new (&edges[i])set<vector<int>>(faceNum);
    }
    
    
    //compute the components class, the components are not connected yet
//    heightField(N_faces, d, components, edges, F, total);
    
   
    
    
    // Plot the mesh
    igl::viewer::Viewer viewer;
    viewer.core.show_lines = false;
    viewer.data.set_mesh(V, F);
    viewer.data.set_normals(N_faces);
    
    
    
    viewer.launch();
    
    for(int i = 0;i < total;i++){
        components[i].~set();
    }
    for(int i = 0;i < total;i++){
        edges[i].~set();
    }
    operator delete(rawMemory);
    operator delete(rawMemory1);
}

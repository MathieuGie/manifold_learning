#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include <igl/jet.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <stack>

#include "spiral.cpp"
#include "eps_isomap.cpp"
#include "knn.cpp"
#include "simple_MDS.cpp"

using namespace Eigen;
using namespace std;


/*
Manifold learning

in the data use the swiss roll

1. Implement simple isomap method
2. Implement the recent variant 
*/


bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;
  if (key == 'D')
  {std::cout<<"it's useless to press 1";}
  return false;}

void set_pc(igl::opengl::glfw::Viewer &viewer) {
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  //viewer.data().add_points(V1,Eigen::RowVector3d(1,0.1,0.1)); // Eigen::RowVector3d(0.1,0.1,0.1) allows to put a color
  //std::cout<<give_color(REDDER(0), V1, 0);
  
  for (int i=0; i<V1.rows();i++){
    double lev_red;
    lev_red=give_color(REDDER(0), V1, i);
    double lev_green;
    lev_green=give_color(GREENER(0), V1, i);
    double lev_blue;
    lev_blue=give_color(BLUER(0), V1, i);
    
    viewer.data().add_points(V1.row(i),Eigen::RowVector3d(lev_red,lev_green,lev_blue));
  }
}

void displayer(){
  //igl::readOFF("../data/sphere.off", V1, F1); //can be replaced 
  
  //REGULAR SPIRAL: (comment otherwise)
  //V1=extend_on_z(make_spiral_equilibrated(5,15,1),2,5);

  //SPIRAL WITH HOLE: (comment otherwise)
  //V1=extend_on_z_hole(make_spiral_equilibrated(5,7,0.5),5,15);

  //WITH MDS: (comment otherwise)
  //V1=MDS(extend_on_z(make_spiral_equilibrated(5,10,0.5),10,10));

  //FIRST ISOMAP - with epsilon method: (comment otherwise)
  /*
  V1=little_MDS(Gramm(extend_on_z(make_spiral_equilibrated(5,7,1),10,10),1.2));
  
  MatrixXd V2d (V1.rows(), 3);
  std::cout<<"V2d"<<"\n";
  for (int i=0;i<V1.rows();i++){
      V2d(i,0)=V1(i,0);
      V2d(i,1)=V1(i,1);
      V2d(i,2)=0;
  }
  V1=V2d;
  */

  //SECOND ISOMAP - with knn method: (comment otherwise)
  
  V1=little_mMDS(Prepare_metric(extend_on_z_hole(make_spiral_equilibrated(5,8,0.5),5,15),20, 50));
  
  MatrixXd V2d (V1.rows(), 3);
  std::cout<<"V2d"<<"\n";
  for (int i=0;i<V1.rows();i++){
      V2d(i,0)=V1(i,0);
      V2d(i,1)=V1(i,1);
      V2d(i,2)=0;
  }
  V1=V2d;
  
  

  igl::opengl::glfw::Viewer viewer;
  set_pc(viewer);
  viewer.launch();
}

MatrixXd attempt(){
MatrixXd attempt(7,3);
attempt(0,0)=0;
attempt(0,1)=0;
attempt(0,2)=0;
attempt(1,0)=0;
attempt(1,1)=1;
attempt(1,2)=0;
attempt(2,0)=1;
attempt(2,1)=1;
attempt(2,2)=0;
attempt(3,0)=2;
attempt(3,1)=0;
attempt(3,2)=0;
attempt(4,0)=6;
attempt(4,1)=2;
attempt(4,2)=0;
attempt(5,0)=7;
attempt(5,1)=2;
attempt(5,2)=0;
attempt(6,0)=8;
attempt(6,1)=2;
attempt(6,2)=0;
return attempt;}

int main(int argc, char *argv[]) {
  displayer();
  //std::cout<<Dijstra(KNN(V1,20,50),1);
  //std::cout<<Dijstra(KNN(attempt(), 2, attempt().rows()), 4);
  }
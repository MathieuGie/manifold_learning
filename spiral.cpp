#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

/* Make the cheese roll to display it in main */


MatrixXd V1; 
MatrixXi F1;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

MatrixXd make_spiral(int start, int stop, int precision){
  MatrixXd Spiral(precision*(stop-start), 2);
  double radius=start;
  double d_radius;
  d_radius=double(1/double(precision));
  for(int i=0; i<precision*(stop-start); i++){
    double cosine;
    double sine;
    cosine=cos(double(double(i+start)/double(precision)));
    sine=sin(double(double(i+start)/double(precision)));
    std::cout<<radius<<i<<"\n";
    Spiral(i,0)=radius*cosine;
    Spiral(i,1)=radius*sine;
    radius+=d_radius;
  }
  return Spiral;
}

MatrixXd make_spiral_equilibrated(int start, int stop, double length){//the precision parameter becomes the length
  int count=1;
  double alpha;
  alpha=double(start);
  while(alpha<=stop){
    alpha+=double(length)/alpha;
    count++;
  }
   MatrixXd Spiral(count, 2);
   alpha=double(start);
   for(int i=0; i<count; i++){
    double cosine;
    double sine;
    cosine=cos(alpha);
    sine=sin(alpha);
    Spiral(i,0)=alpha*cosine;
    Spiral(i,1)=alpha*sine;
    alpha+=double(length)/alpha;
   }
  return Spiral;
}

MatrixXd extend_on_z(MatrixXd spiral, int length, int quantity){
  MatrixXd Roll(quantity*spiral.rows(),3);
  for (int i=0; i<spiral.rows();i++){
    for (int j=0; j<quantity; j++){
      Roll(i*quantity+j,0)=spiral(i,0);
      Roll(i*quantity+j,1)=spiral(i,1);
      Roll(i*quantity+j,2)=fRand(0,length);
    }
  }
  return Roll;
}

MatrixXd extend_on_z_hole(MatrixXd spiral, int length, int quantity){
  
  double total_length=spiral.rows();
  double start1=total_length*2/5;
  double stop1=total_length*3/5;
  int start=int (start1);
  int stop=int (stop1);
  //std::cout<<"start"<<"\t"<<total_length<<"\t"<<start1<<"\t"<<start;
  int small_quant=quantity-int(double(quantity)/2);
  int len=(stop-start)*small_quant+(total_length-(stop-start))*quantity;
  //std::cout<<len<<"\t"<<start<<"\t"<<stop<<"\n";

  MatrixXd Roll(len,3);

  for (int i=0; i<spiral.rows();i++){
    if (i>=start){
      if (i<stop){
          for (int j=0; j<small_quant;j++){
          //std::cout<<start*quantity+(i-start)*small_quant+j<<"\t"<<i<<"\n";
          Roll(start*quantity+(i-start)*small_quant+j,0)=spiral(i,0);
          Roll(start*quantity+(i-start)*small_quant+j,1)=spiral(i,1);
          int pick_side=fRand(0,10);
          //std::cout<<pick_side<<"\t";
          if (pick_side<5){
            Roll(start*quantity+(i-start)*small_quant+j,2)=fRand(0,double(length)/double(5));
          } else {
            Roll(start*quantity+(i-start)*small_quant+j,2)=fRand(4*double(length)/double(5),length);
          }
          }
        } else{
          for (int j=0; j<quantity; j++){
          //std::cout<<"normal"<<"\t"<<start*quantity+(stop-start)*small_quant+(i-stop)*quantity+j<<"\t"<<i<<"\n";
          Roll(start*quantity+(stop-start)*small_quant+(i-stop)*quantity+j,0)=spiral(i,0);
          Roll(start*quantity+(stop-start)*small_quant+(i-stop)*quantity+j,1)=spiral(i,1);
          Roll(start*quantity+(stop-start)*small_quant+(i-stop)*quantity+j,2)=fRand(0,length);
        }}
      } else{
          for (int j=0; j<quantity; j++){
          //std::cout<<"normal"<<"\t"<<i*quantity+j<<"\t"<<i<<"\n";
          Roll(i*quantity+j,0)=spiral(i,0);
          Roll(i*quantity+j,1)=spiral(i,1);
          Roll(i*quantity+j,2)=fRand(0,length);
          }
        }
  }

  return Roll;
}

double give_color(MatrixXd colors, MatrixXd V, int index){
  MatrixXd Interp(colors.rows()-1,2);
  for (int i=0;i<Interp.rows();i++){
    //std::cout<<"hellooooo"<<i<<"\n"<<Interp.rows()<<"\t"<<Interp<<"\n";
    Interp(i,0)=double((colors(i+1,0)-colors(i,0)));
    Interp(i,1)=double(colors(i,0))-double(Interp(i,0))*double(i);
  }
  //std::cout<<Interp<<"\n";
  double ratio=double(index)*double(Interp.rows())/double(V.rows());
  int index_interp=floor(ratio);
  //std::cout<<ratio<<"\t"<<double(Interp(index_interp,0))*double(ratio)+double(Interp(index_interp,1))<<"\n";
  return double(Interp(index_interp,0))*double(ratio)+double(Interp(index_interp,1));
}

MatrixXd REDDER(int a){
MatrixXd COLOR(7,1);
COLOR(0,0)=double(148)/double(255);
COLOR(1,0)=double(75)/double(255);
COLOR(2,0)=double(0)/double(255);
COLOR(3,0)=double(0)/double(255);
COLOR(4,0)=double(1);
COLOR(5,0)=double(1);
COLOR(6,0)=double(1);
return COLOR;
}

MatrixXd GREENER(int a){
MatrixXd COLOR(7,1);
COLOR(0,0)=double(0)/double(255);
COLOR(1,0)=double(0)/double(255);
COLOR(2,0)=double(0)/double(255);
COLOR(3,0)=double(1);
COLOR(4,0)=double(1);
COLOR(5,0)=double(127)/double(255);
COLOR(6,0)=double(0)/double(255);
return COLOR;
}

MatrixXd BLUER(int a){
MatrixXd COLOR(7,1);
COLOR(0,0)=double(211)/double(255);
COLOR(1,0)=double(130)/double(255);
COLOR(2,0)=double(255)/double(255);
COLOR(3,0)=double(0)/double(255);
COLOR(4,0)=double(0)/double(255);
COLOR(5,0)=double(0)/double(255);
COLOR(6,0)=double(0)/double(255);
return COLOR;
}

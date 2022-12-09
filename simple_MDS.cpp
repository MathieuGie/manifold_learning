#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>

#include <list>

using namespace Eigen;
using namespace std;

/* Here implement the MDS algorithm */

MatrixXd MDS(MatrixXd P){
    Vector3d bar = P.colwise().mean();
    std::cout<<"bar \n"<<bar;
    MatrixXd P_bar (P.rows(),P.cols());
    for (int i=0;i<P.rows();i++){
        P_bar(i,0)=P(i,0)-bar(0);
        P_bar(i,1)=P(i,1)-bar(1);
        P_bar(i,2)=P(i,2)-bar(2);
    }
    //std::cout<<"P_bar \n"<<P_bar.row(0);
    MatrixXd G(P.rows(),P.rows());
    G=P_bar*P_bar.transpose();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(G);

    MatrixXd Lambda(2,1);
    Lambda(0,0)=eigensolver.eigenvalues()(P.rows()-1);
    Lambda(1,0)=eigensolver.eigenvalues()(P.rows()-2);

    std::cout<<"Lambda"<<Lambda;

    MatrixXd Q (P.rows(),P.rows());
    Q=eigensolver.eigenvectors();

    MatrixXd result(P.rows(),2);
    for(int i=0;i<P.rows();i++){
        for (int j=0;j<P.rows();j++){
            if (j==0){
                result(i,j)=Q(i,0)*Lambda(0);
            } else if (j==1){
                result(i,j)=Q(i,1)*Lambda(1);
            }
        }
    }
    return result;
}
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include <stack>

#include <list>

using namespace Eigen;
using namespace std;

/* Here implement the first isomap method using the k-nearest neigbor */

//L2 dist and little MDS are imported from eps_isomap so
//#include "eps_isomap.cpp"

MatrixXd KNN(MatrixXd V, int k, int lookout){
    MatrixXd knn(V.rows(),V.rows());
    for (int i=0; i<V.rows(); i++){
      for (int j=0; j<V.rows(); j++){
        knn(i,j)=0;
      }
    }
    std::cout<<V.rows()<<"\n";
    for(int i=0; i<V.rows(); i++){
        //std::cout<<i<<"\n";
        int half=lookout/2;
        int before=half;
        int after=half;
        if (i<half){
            before=i;
            after=lookout-i;
        } else if (V.rows()-i<half){
            before=lookout-V.rows()+i;
            after=V.rows()-i;
        }
        //std::cout<<before<<"\t"<<after<<"\n";

        MatrixXd max_corresp(k, 2); //col 0, the index of point in V, col1, distance with the point considered

        //first initialise the max_corresp to 1000 everywhere (like infinity)
        for (int line=0; line<k; line++){
            max_corresp(line, 0)=0;
            max_corresp(line, 1)=1000;
        }
        //std::cout<<max_corresp<<"\n";

        for (int neighb=i-before; neighb<i+after; neighb++){
            //std::cout<<neighb<<"\n";
            if (neighb!=i){
                double distance=L2_dist(V.row(neighb), V.row(i));
                //std::cout<<max_corresp.col(1).maxCoeff()<<"\n";

                //Find index at which one should insert the new neighbor
                if (distance<max_corresp.col(1).maxCoeff()){
                    int index=0;
                    while (max_corresp(index,1)>distance){
                        index+=1;
                        if (index==k){
                            break;
                        }
                    }
                    index-=1;

                    //Update max_corresp which tracks the closest neighbors
                    MatrixXd new_max_corresp(k, 2);
                    for (int line=0; line<k; line++){
                        //std::cout<<line<<"\n";
                        if (line>index){
                            new_max_corresp(line, 0)=max_corresp(line, 0);
                            new_max_corresp(line, 1)=max_corresp(line, 1);
                        } else if (line==index){
                            new_max_corresp(line, 0)=neighb;
                            new_max_corresp(line, 1)=distance;
                        } else {
                            new_max_corresp(line, 0)=max_corresp(line+1, 0);
                            new_max_corresp(line, 1)=max_corresp(line+1, 1);
                        }
                    } max_corresp=new_max_corresp;
                }
       
            } //std::cout<<max_corresp<<"\n";
        }
        for (int l=0; l<k; l++){
            knn(i,max_corresp(l,0))=1; //instead of l+1
            knn(max_corresp(l,0),i)=1;
        }
    }
    return knn;
}

MatrixXd Gramm2(MatrixXd V, int k, int lookout){
    MatrixXd Edges(V.rows(), V.rows());
    Edges=KNN(V, k, lookout);
    MatrixXd gram(V.rows(),V.rows());

    for (int i=0; i<V.rows(); i++){
      int l=V.rows()-i-1;
      MatrixXd DIJ(V.rows(), 3);
      DIJ=Dijstra(Edges, l);
      std::cout<<"done Dij"<<i<<"\n";
      for (int j=0; j<V.rows(); j++){
          gram(l,j)=DIJ(j,1);
          gram(j,l)=DIJ(j,1);
        }
    }
    //std::cout<<gram<<"\n";
    return gram;
}

MatrixXd Prepare_metric(MatrixXd V, int k, int lookout){
  MatrixXd Edges(V.rows(),V.rows());
  Edges=Gramm2(V, k, lookout);
  //std::cout<<Edges;
  MatrixXd e(V.rows(), 1);
  MatrixXd I(V.rows(), V.rows());
  for (int i=0; i<V.rows(); i++){
    e(i)=1;
    for (int j=0; j<V.rows(); j++){
      Edges(i,j)=pow(Edges(i,j),2);

      if(i==j){
        I(i,j)=1;
      }else{
        I(i,j)=0;
      }
    }
  }
  MatrixXd H(V.rows(), V.rows());
  H=I - (double(1)/double(V.rows())) * e * e.transpose();
  //std::cout<<"helooooooo"<<Edges.row(0)<<"\n";
  //std::cout<<"helooooooo"<<Edges.row(V.rows()-1)<<"\n";
  //std::cout<<-0.5*H*Edges*H;
  return -0.5*H*Edges*H;
}

MatrixXd little_mMDS(MatrixXd gram){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(gram);
    //std::cout<<gram;

    MatrixXd Lambda(2,1);
    Lambda(0,0)=eigensolver.eigenvalues()(gram.rows()-1);
    Lambda(1,0)=eigensolver.eigenvalues()(gram.rows()-2);
    std::cout<<eigensolver.eigenvalues()<<"\n";
    //std::cout<<Lambda;
    //std::cout<<eigensolver.eigenvalues();

    MatrixXd Q (gram.rows(),gram.rows());
    Q=eigensolver.eigenvectors();
    //std::cout<<"Q"<<"\t"<<Q.row(0)<<"\n";
    //std::cout<<"eigen"<<"\t"<<eigensolver.eigenvectors()<<"\n";

    MatrixXd result(gram.rows(),2);
    for(int i=0;i<gram.rows();i++){
        for (int j=0;j<gram.rows();j++){
            if (j==0){
                result(i,j)=Q(i,gram.rows()-1)*pow(Lambda(0,0),0.5);
            } else if (j==1){
                result(i,j)=Q(i,gram.rows()-2)*pow(Lambda(1,0),0.5);
            }
        }
    }
    return result;
}











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

/* Here implement the first isomap method using the epsilon neigbor */

//Some basic functions:

double L2_dist(MatrixXd A, MatrixXd B){
    double distance=0;
    for(int i=0;i<A.cols();i++){
        distance+=pow((A(0,i)-B(0,i)),2);
    }
    return pow(distance, 0.5);
}

MatrixXd make_fake(MatrixXd A){
    int fake=100000000;
    for (int i=0; i<A.rows(); i++){
        for (int j=0; j<A.cols(); j++){
            A(i,j)=fake;
        }
    }
    return A;
}

int check_ele_mtx(MatrixXd A, int element){
    for (int i=0; i<A.rows(); i++){
        for (int j=0; j<A.cols(); j++){
            if (A(i,j)==element){
                return 1;
            }
        }
    }
    return 0;
}

//Build the edges matrix according to epsilon from V

MatrixXd eps_edges(MatrixXd V, double epsilon){
    MatrixXd Edges(V.rows(),V.rows());
    for(int i=0;i<V.rows();i++){
        for(int j=0; j<V.rows(); j++){
            if (L2_dist(V.row(i),V.row(j))<epsilon && i!=j){
                Edges(i,j)=1;
                Edges(j,i)=1;
            } else{
               Edges(i,j)=0;
               Edges(j,i)=0; 
            }
        }
    }
    std::cout<<Edges.rows();
    for (int i=0; i<Edges.rows();i++){
        Eigen::Index maxRow, maxCol;
        if ((Edges.row(i).maxCoeff)(&maxRow, &maxCol)==0){
            std::cout<<"alert"<<"\t"<<i<<"\n";
        }
    }
    return Edges;
}

//Implement Dijstra's algorithm

MatrixXd Dijstra(MatrixXd Edges, int index_start){
    int fake=100000000;
    MatrixXd tracker(Edges.rows(),3);
    for (int j=0;j<Edges.rows(); j++){
        tracker(j,0)=j;
        tracker(j,1)=fake;
        tracker(j,2)=fake; //this number means there is no one previously
        }

    tracker(index_start, 1)=0;
    tracker(index_start, 2)=fake;

    //so we will start with index_start
    MatrixXd see_later(Edges.rows(),1);
    see_later=make_fake(see_later);
    see_later(index_start)=1;

    MatrixXd explored(Edges.rows(), 1);
    explored=make_fake(explored);

    MatrixXd final (Edges.rows(), 1);
    for (int i=0; i<Edges.rows(); i++){
        final(i,0)=i;
    }
    int repetition=0;
    int studied_before;
    studied_before=100000000;

    while(explored!=final){
        
        //Find a suitable study point
        Eigen::Index minRow, minCol;
        float min = see_later.minCoeff(&minRow,&minCol);
        int study_point=minRow;
        see_later=make_fake(see_later);
        
        //std::cout<<study_point<<"\t"<<studied_before<<"\n";

        //To prevent looping around and missing a point: trick to jump in the graph to a non-studied point
        if (repetition>2){
            Eigen::Index maxRow, maxCol;
            float max = explored.maxCoeff(&maxRow,&maxCol);
            study_point=maxRow;
        }

        //Updating the tracker
        for(int index=0; index<Edges.rows(); index++){
            if (Edges(index, study_point)>0 && check_ele_mtx(explored, index)==0){
                see_later(index,0)=Edges(index, study_point);
                if (tracker(index,2)==fake){
                    tracker(index,1)=tracker(study_point, 1)+Edges(index, study_point);
                    tracker(index,2)=study_point;
                } else if (tracker(index, 1)>(tracker(study_point, 1)+Edges(index, study_point))){
                    tracker(index,1)=tracker(study_point, 1)+Edges(index, study_point);
                    tracker(index,2)=study_point;
                }
            }
        }
        explored(study_point,0)=study_point;
        
        if (studied_before==study_point){
            repetition+=1;
        }
        studied_before=study_point;

        //std::cout<<"explored"<<"\t"<<explored<<"\n";
        //std::cout<<"see_later"<<"\t"<<see_later<<"\n";
    }
    //std::cout<<tracker<<"\n";
    //std::cout<<explored;
    return tracker;
}

//Implement Gram Matrix and from it do MDS


MatrixXd Gramm(MatrixXd V, double eps){
    MatrixXd Edges(V.rows(), V.rows());
    Edges=eps_edges(V, 1.5);
    MatrixXd gram(V.rows(),V.rows());

    for (int i=0; i<V.rows(); i++){

        MatrixXd DIJ(V.rows(), 3);
        DIJ=Dijstra(Edges, i);
        std::cout<<"done Dij"<<i<<"\n";
        for (int j=0; j<V.rows(); j++){
            gram(i,j)=DIJ(j,1);
            gram(j,i)=DIJ(j,1);
        }
    }
    //std::cout<<"done gram"<<"\n";
    return gram;
}

MatrixXd little_MDS(MatrixXd gram){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(gram);
    //std::cout<<gram;

    MatrixXd Lambda(2,1);
    Lambda(0,0)=eigensolver.eigenvalues()(0);
    Lambda(1,0)=eigensolver.eigenvalues()(1);
    //std::cout<<eigensolver.eigenvalues()<<"\n";
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
                result(i,j)=Q(i,0)*pow(-Lambda(0,0),0.5);
            } else if (j==1){
                result(i,j)=Q(i,1)*pow(-Lambda(1,0),0.5);
            }
        }
    }
    return result;
}



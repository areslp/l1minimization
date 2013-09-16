#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <omp.h>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <fstream>
#include <vector>

#include <Eigen/Sparse>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>

using namespace Eigen;
using namespace std;
using namespace pcl;

void mxArrayToEigen(const mxArray *m, MatrixXd& A)
{
    const int nRows = mxGetM(m); // number of rows
    const int nCols = mxGetN(m); // number of columns
    A.resize(nRows, nCols);
    const double *p = mxGetPr(m);
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            A(i, j) = p[j * nRows + i];
        }
    }
}

void mat2pcl(const mxArray* array, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud){
    MatrixXd mat;    
    mxArrayToEigen(array, mat);
    int n=mat.rows();
    assert(mat.cols()==3);
    for(int i=0;i<n;i++){
        cloud->points.push_back(pcl::PointXYZ(mat(i,0),mat(i,1),mat(i,2)));
    }
}

// 返回mat'的三元组，行列是反的
void EigenSparseTomxArray(Eigen::SparseMatrix<double>& mat,mxArray* ia, mxArray* ja, mxArray* va){
    double *pi=mxGetPr(ia);
    double *pj=mxGetPr(ja);
    double *pv=mxGetPr(va);
    int idx=0;
    for (int k=0; k<mat.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
            double v=it.value();
            int row=it.row();   // row index
            int col=it.col();   // col index (here it is equal to k)
            pi[idx]=col;
            pj[idx]=row;
            pv[idx]=v;
            ++idx;
        }
}

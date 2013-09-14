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

// void mxArrayToEigen(const mxArray *m, MatrixXd& A)
// {
    // const int nRows = mxGetM(m); // number of rows
    // const int nCols = mxGetN(m); // number of columns
    // A.resize(nRows, nCols);
    // const double *p = mxGetPr(m);
    // for (int i = 0; i < nRows; i++) {
        // for (int j = 0; j < nCols; j++) {
            // A(i, j) = p[j * nRows + i];
        // }
    // }
// }


void EigenSparseTomxArray(Eigen::SparseMatrix<double>& mat,mxArray* array){
    double *pa=mxGetPr(array);
    int m=mxGetM(array);
    int n=mxGetN(array);
    for (int k=0; k<mat.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
            double v=it.value();
            int row=it.row();   // row index
            int col=it.col();   // col index (here it is equal to k)
            pa[col*m+row]=v; 
        }
}

#include "util.h"
#include <pcl/common/time.h>

// parameters points,k 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// declare variables
	const mxArray *points;
    int k;
	//associate inputs
    points= (prhs[0]);
	k= mxGetScalar(prhs[1]);
    int dim=3;
    int n=mxGetM(points);
    int Asize=n*dim;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    mat2pcl(points,cloud);

    //kdtree
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud);
    //计算平均边长
    // int kk=3;
    // double av_len=0.0;
    // for (int i=0;i<n;i++){
        // double av_len_p=0.0;
        // // get k-nearest neighbors
        // std::vector<int> pointIdxNKNSearch(kk);
        // std::vector<float> pointNKNSquaredDistance(kk);
        // kdtree.nearestKSearch (cloud->points[i], kk+1, pointIdxNKNSearch, pointNKNSquaredDistance); //排除当前点自己
        // for (int j=0;j<(int)pointNKNSquaredDistance.size();j++)
        // {
            // if (pointIdxNKNSearch[j]==i)
            // {
                // continue;
            // }
            // av_len_p+=pointNKNSquaredDistance[j];
        // }
        // av_len_p/=kk;
        // av_len+=av_len_p;
    // }
    // av_len/=n;
    // printf("av_len:%f\n",av_len);

    int nnz=(k+1)*Asize;
    SparseMatrix<double> A(Asize,Asize);
    A.reserve(VectorXi::Constant(Asize,k+1)); //预分配空间

    // double sigma=2*av_len;
    StopWatch timer;
    double t;
    omp_set_num_threads(4);
    #pragma omp parallel for
    for(int i=0;i<n;i++){
        // timer.reset();
        std::vector<int> pointIdxNKNSearch(k);
        std::vector<float> pointNKNSquaredDistance(k);
        kdtree.nearestKSearch (cloud->points[i], k+1, pointIdxNKNSearch, pointNKNSquaredDistance); 
        // double sum=0;
        int kn=pointNKNSquaredDistance.size();
        assert(kn==k+1);
        //exp(-d.^2/sigma^2) and sum for normalize
        // for(int j=0;j<kn;j++){
            // float v=pointNKNSquaredDistance[j];
            // v=exp(-v*v/(sigma*sigma));
            // pointNKNSquaredDistance[j]=v;
            // sum+=v*v;
        // }
        // sum=sqrt(sum);
        // t=timer.getTime();
        // printf("kdtree takes:%f\n",t);
        // timer.reset();
        float v=1.0/(float)kn;
        for(int j=0;j<kn;j++){
            int jidx=pointIdxNKNSearch[j];
            // float v=pointNKNSquaredDistance[j];
            // v=v/sum;
            for(int dd=0;dd<dim;++dd){
                // insert rows for each dim
                A.insert(jidx*dim+dd,i*dim+dd)=v;
            }
        }
        // t=timer.getTime();
        // printf("SparseMatrix insert takes:%f\n",t);
    }

    //output
    mxArray* ia=mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* ja=mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* va=mxCreateDoubleMatrix(nnz,1,mxREAL);
    EigenSparseTomxArray(A,ia,ja,va);
    nlhs=3;
    plhs[0] = ia;
    plhs[1] = ja;
    plhs[2] = va;
}

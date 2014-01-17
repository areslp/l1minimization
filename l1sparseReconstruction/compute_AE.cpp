#include "util.h"
#include <pcl/common/time.h>

void construct_E(mxArray* E,SparseMatrix<double>& A);

void construct_E(mxArray* E,SparseMatrix<double>& mat){
    int m=mxGetM(E);
    int n=mxGetM(E);
    printf("row of E: %d\n",m);
    double* pe=mxGetPr(E);
    int rowe=1;
    for (int k=0; k<mat.outerSize(); ++k){
        vector<int> idxs;
        for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
            int row=it.row();   // row index
            int col=it.col();   // col index (here it is equal to k)
            idxs.push_back(row);
        }
        assert(idxs.size()==2);
        pe[rowe-1]=idxs[0];
        pe[m+rowe-1]=idxs[1];
        ++rowe;
    }
}

// parameters points,normals,k,dim
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// declare variables
	const mxArray *points, *normals;
    int k,dim;
	//associate inputs
    points= (prhs[0]);
    normals= (prhs[1]);
	k= mxGetScalar(prhs[2]);
	dim= mxGetScalar(prhs[3]);
    int n=mxGetM(points);
    int total_pair=n*k*dim;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudN(new pcl::PointCloud<pcl::PointXYZ>);
    mat2pcl(points,cloud);
    mat2pcl(normals,cloudN);

    //kdtree
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud);

    int nnz=total_pair*2;
    SparseMatrix<double> A(n*dim,total_pair);
    A.reserve(VectorXi::Constant(total_pair,k+1)); //预分配空间

    //计算平均边长
    int kk=6; //curve实验里用两个点计算平均边长，point cloud里用4-6个点
    double av_len=0.0;
    for (int i=0;i<n;i++){
        double av_len_p=0.0;
        // get k-nearest neighbors
        std::vector<int> pointIdxNKNSearch(kk);
        std::vector<float> pointNKNSquaredDistance(kk);
        kdtree.nearestKSearch (cloud->points[i], kk+1, pointIdxNKNSearch, pointNKNSquaredDistance); //排除当前点自己
        for(int j=0; j < pointNKNSquaredDistance.size(); j++)
        {
            if(pointIdxNKNSearch[j]==i){
                continue;
            }
            av_len_p+=sqrt(pointNKNSquaredDistance[j]);
        }
        av_len_p/=kk;
        av_len+=av_len_p;
    }
    av_len/=n;
    printf("av_len:%f\n",av_len);

    double sigma=2*av_len;

    StopWatch timer;
    double t;
    omp_set_num_threads(4);
    #pragma omp parallel for
    for(int i=0;i<n;i++){
        // timer.reset();
        std::vector<int> pointIdxNKNSearch(k);
        std::vector<float> pointNKNSquaredDistance(k);
        pcl::PointXYZ query_point=cloud->points[i];
        kdtree.nearestKSearch (query_point, k+1, pointIdxNKNSearch, pointNKNSquaredDistance); 

        int kn=pointNKNSquaredDistance.size();
        assert(kn==k+1);

        // for ith point, var should be normalized
        double tv=0;
        for(int j=0;j<kn;j++){
            int jidx=pointIdxNKNSearch[j];
            if ( jidx==i )
            {
                continue;
            }
            pcl::PointXYZ tp(query_point.x-cloud->points[jidx].x,query_point.y-cloud->points[jidx].y,query_point.z-cloud->points[jidx].z); 
            double dis=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
            double var=exp(-dis*dis/(sigma*sigma));
            tv+=var;
        }

        int jc=0;
        for(int j=0;j<kn;j++){
            int jidx=pointIdxNKNSearch[j];
            // 跳过当前点
            if ( jidx==i )
            {
                continue;
            }
            pcl::PointXYZ tp(query_point.x-cloud->points[jidx].x,query_point.y-cloud->points[jidx].y,query_point.z-cloud->points[jidx].z); 
            double dis=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
            double var=exp(-dis*dis/(sigma*sigma));
            // var=var/tv; // 对第i点进行归一化，这个貌似不能加，加了法矢、位置都不对了
            for(int dd=0;dd<dim;++dd){
                int col=(i*k+jc)*dim+dd; //这里不能用j，因为有一个continue的，j还是会++
                A.insert(jidx*dim+dd,col)=-var;
                A.insert(i*dim+dd,col)=var;
                // printf("inserting to A %d col\n",col);
            }
            jc++;
        }
        // t=timer.getTime();
        // printf("SparseMatrix insert takes:%f\n",t);
    }
    printf("caculate A complete\n");
    //output
    //
    //Sparse matrix A
    mxArray* ia=mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* ja=mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* va=mxCreateDoubleMatrix(nnz,1,mxREAL);
    EigenSparseTomxArray(A,ia,ja,va);
    //E
    mxArray* E=mxCreateDoubleMatrix(total_pair,2,mxREAL);
    construct_E(E,A);
    nlhs=4;
    plhs[0] = ia;
    plhs[1] = ja;
    plhs[2] = va;
    plhs[3] = E;
}

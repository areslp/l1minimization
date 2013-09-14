#include "util.h"


void mat2pcl(const mxArray* array, pcl::PointCloud<pcl::PointXYZ>& cloud);

void mat2pcl(const mxArray* array, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud){
    MatrixXd mat;    
    // mxArrayToEigen(array, mat);
    int n=mat.rows();
    assert(mat.cols()==3);
    for(int i=0;i<n;i++){
        cloud->points.push_back(pcl::PointXYZ(mat(i,0),mat(i,1),mat(i,2)));
    }
}


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
    // mat2pcl(points,cloud);

    //kdtree
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud);
    //计算平均边长
    int kk=3;
    double av_len=0.0;
    for (int i=0;i<n;i++){
        double av_len_p=0.0;
        // get k-nearest neighbors
        std::vector<int> pointIdxNKNSearch(kk);
        std::vector<float> pointNKNSquaredDistance(kk);
        kdtree.nearestKSearch (cloud->points[i], kk+1, pointIdxNKNSearch, pointNKNSquaredDistance); //排除当前点自己
        for (int j=0;j<(int)pointNKNSquaredDistance.size();j++)
        {
            if (pointIdxNKNSearch[j]==i)
            {
                continue;
            }
            av_len_p+=pointNKNSquaredDistance[j];
        }
        av_len_p/=kk;
        av_len+=av_len_p;
    }
    av_len/=n;

    int nnz=(k+1)*Asize;
    SparseMatrix<double> A(Asize,Asize);
    A.reserve(VectorXi::Constant(Asize,k+1)); //预分配空间

    double sigma=2*av_len;
    omp_set_num_threads(4);
    #pragma omp parallel for
    for(int i=0;i<n;i++){
        std::vector<int> pointIdxNKNSearch(k);
        std::vector<float> pointNKNSquaredDistance(k);
        kdtree.nearestKSearch (cloud->points[i], k+1, pointIdxNKNSearch, pointNKNSquaredDistance); 
        double sum=0;

        int kn=pointNKNSquaredDistance.size();
        //exp(-d.^2/sigma^2) and sum for normalize
        for(int j=0;j<kn;j++){
            float v=pointNKNSquaredDistance[j];
            v=exp(-v*v/(sigma*sigma));
            pointNKNSquaredDistance[j]=v;
            sum+=v*v;
        }
        sum=sqrt(sum);
        for(int j=0;j<kn;j++){
            float v=pointNKNSquaredDistance[j];
            pointNKNSquaredDistance[j]=v/sum;
        }
        for(int j=0;j<kn;j++){
            A.insert(i,pointIdxNKNSearch[j])=pointNKNSquaredDistance[j];
        }
    }

    //output
    // mxArray* out=mxCreateSparse(Asize,Asize,nnz,mxREAL);
    // EigenSparseTomxArray(A,out);
    // plhs[0] = out;
}

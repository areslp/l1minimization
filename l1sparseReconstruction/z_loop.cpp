#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <omp.h>
#include <iostream>
// #include <windows.h>
#include <time.h>
#include <assert.h>
#include <fstream>
#include <vector>
using namespace std;


double norm(double* x,int size);
mxArray* compute_x(mxArray* V, mxArray* q, mxArray* D, double t);
mxArray* get_column(mxArray* z, int index);
void swap_column(mxArray* z,int index,mxArray* c);
mxArray* Atb(mxArray* A,mxArray* b);
mxArray* AB(mxArray* A,mxArray* B);
mxArray* update_x(mxArray* A, mxArray* b, double kappa, mxArray* V, mxArray* Dig);
void vst(double k, double* v, double* r,int size);
void get_sub_vec(int start, double* v, int dim, double* r); //得到向量v的子向量，从start开始取dim个
void set_sub_vec(double* src, double* dst, int start, int size); //从start开始，把src的内容拷贝到dst里

void set_sub_vec(double* src, double* dst, int start, int size){
    for(int i=0;i<size;i++){
        dst[start+i]=src[i];
    }
}

void get_sub_vec(int start, double* v, int dim, double* r){
    for(int i=0;i<dim;i++){
        r[i]=v[start+i];
    }
}

void vst(double k, double* v, double* r, int size){
    double normv=norm(v, size);
    if(normv==0){
        for(int i=0;i<size;i++){
            r[i]=0;
        }
    } else {
        double co=1-k/normv;
        if(co<0){
            co=0;
        }
        for(int i=0;i<size;i++){
            r[i]=v[i]*co;
        }
    }
}

double norm(double* x,int size){
    double sum=0;
    for(int i=0;i<size;i++){
        sum+=x[i]*x[i];
    }
    sum=sqrt(sum);
    return sum;
}

// parameters B=Dx-u,lambda/rho,dim
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// declare variables
	mxArray *z, *B, *V, *Dig;
    double kappa;
	int N, dim;

	//associate inputs
    B= const_cast<mxArray*>(prhs[0]);
	kappa= mxGetScalar(prhs[1]);
	dim= mxGetScalar(prhs[2]);

	// output, z is const, cannot modify
	int mz = (int)mxGetM(B);
	int nz = (int)mxGetN(B);
    N=mz/dim;
    mxArray* zz = mxCreateDoubleMatrix(mz, nz, mxREAL);
    plhs[0] = zz;
    // get pointer
    double *pzz=mxGetPr(zz);
    double *pB=mxGetPr(B);
    
    double b[3];
    double xx[3];
    omp_set_num_threads(4);
    #pragma omp parallel for private(b,xx)
	for (int i=0;i<N;i++)
	{
        int idx=i*dim;
        get_sub_vec(idx,pB,dim,b); //取B的第i列
        vst(kappa,b,xx,dim);
        set_sub_vec(xx,pzz,idx,dim);
    }
}


// x = V*((V'*q)./(D + t)); nx1
// V nxn
// D nx1
// q nx1
// mxArray* compute_x(mxArray* V, mxArray* q, mxArray* D, double t){
    // int md=(int)mxGetM(D);
    // int nd=(int)mxGetN(D);
    // mxArray* tmp=mxCreateDoubleMatrix(md,nd,mxREAL);
    // mxArray* Vtq=Atb(V,q);
    // double* ptmp=mxGetPr(tmp);
    // double* pvtq=mxGetPr(Vtq);
    // double* pd=mxGetPr(D);
    // for(int i=0;i<md;i++){
        // for(int j=0;j<nd;j++){
            // ptmp[md*j+i]=pvtq[md*j+i]/(pd[md*j+i]+t);
        // }
    // }
    // mxArray* r=AB(V,tmp);
    // mxDestroyArray(tmp);
    // mxDestroyArray(Vtq);
    // return r;
// }

// mxArray* get_column(mxArray* z, int index){
    // int k=(int)mxGetM(z);
    // int n=(int)mxGetN(z);
    // mxArray *c=mxCreateDoubleMatrix(k,1,mxREAL);
    // double *p = mxGetPr(c);
    // double *q = mxGetPr(z);
    // for(int j=0;j<k;j++){
        // p[j]=q[k*index+j];                    
    // }
    // return c;
// }

// void swap_column(mxArray* z,int index,mxArray* c){
    // int k=(int)mxGetM(z);
    // int n=(int)mxGetN(z);
    // double *p = mxGetPr(c);
    // double *q = mxGetPr(z);
    // for(int j=0;j<k;j++){
        // q[k*index+j]=p[j];                    
    // }
// }

// A'*b
// mxArray* Atb(mxArray* A,mxArray* b){
    // // A mxn
    // // b mx1
    // int m=(int)mxGetM(A); 
    // int n=(int)mxGetN(A); 
    // mxArray* r=mxCreateDoubleMatrix(n,1,mxREAL);
    // double *pa=mxGetPr(A);
    // double *pb=mxGetPr(b);
    // double *pr=mxGetPr(r);
    // for(int i=0;i<n;i++){
        // //column of A dot b
        // for(int j=0;j<m;j++){
            // pr[i]+=pa[m*i+j]*pb[j];
        // }
    // }
    // return r;
// }
//  AB
// mxArray* AB(mxArray* A,mxArray* B){
    // int m=(int)mxGetM(A); 
    // int k=(int)mxGetN(A); 
    // int n=(int)mxGetN(B);
    // mxArray* r=mxCreateDoubleMatrix(m,n,mxREAL);
    // double* pa=mxGetPr(A);
    // double* pb=mxGetPr(B);
    // double* pr=mxGetPr(r);
    // for(int i=0;i<m;i++){
        // for(int j=0;j<n;j++){
            // for(int l=0;l<k;l++){
                // pr[m*j+i]+=pa[m*l+i]*pb[j*k+l];
            // }
        // }
    // }
    // return r;
// }

// mxArray* update_x(mxArray* A, mxArray* b, double kappa, mxArray* V, mxArray* Dig){
    // // LARGE_INTEGER  large_interger;  
    // // LONGLONG dff;  
    // // __int64  c1, c2;  

    // // QueryPerformanceFrequency(&large_interger);  
    // // dff = large_interger.QuadPart;  
    // // QueryPerformanceCounter(&large_interger);  
    // // c1 = large_interger.QuadPart;  

    // // [m,n] = size(A);
    // int m=(int)mxGetM(A); 
    // int n=(int)mxGetN(A); 

    // mxArray *input_array[10], *output_array[10]; //mexCallMATLAB 参数

    // mxArray* xx = mxCreateDoubleMatrix(n, 1, mxREAL);
    // // q = A'*b;
    // mxArray* q = Atb(A,b);

    // double tmp=norm(q);

    // if (tmp<=kappa) {
        // //置0
        // double *p = mxGetPr(xx); 
        // for(int i=0;i<n;i++){
            // *p=0;
            // p++;
        // }
    // } else {
        // double lower=0, upper=1e10;
        // for(int i=0;i<100;i++){
            // double t=(upper+lower)/2.0;
            // // x = V*((V'*q)./(D + t)); nx1
            // mxDestroyArray(xx);
            // xx=compute_x(V,q,Dig,t);
            // if(t>kappa/norm(xx))
                // upper=t;
            // else
                // lower=t;
            // if(upper-lower<1e-6)
                // break;
        // }
    // }
    // mxDestroyArray(q);

    // // QueryPerformanceCounter(&large_interger);  
    // // c2 = large_interger.QuadPart;  
    // // printf("times:%lf ms\n", (c2 - c1) * 1000 / dff);  
    // return xx;
// }

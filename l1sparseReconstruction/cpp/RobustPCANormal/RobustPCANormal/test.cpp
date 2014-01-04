#include <iostream>
#include <fstream>
#include <string>

#include <pcl/filters/extract_indices.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>

using namespace std;

int main(int argc, char** argv)
{
	if (argc<=1)
	{
		cout<<"Plz select input mesh!"<<endl;
		return 0;
	}
	//read points
	string filename=argv[1];
	int KN=20;
	if (argc>2)
	{
		KN=boost::lexical_cast<int>(argv[2]);
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	ifstream ifh(filename.c_str());
	float x,y,z;
	while (ifh>>x>>y>>z)
	{
		cloud->push_back(pcl::PointXYZ(x,y,z));
	}
	ifh.close();
	int pnum=cloud->points.size();
	cout<<"point size:"<<pnum<<endl;
	//kdtree
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud(cloud);
	int K=3;
	//计算平均边长
	float av_len=0.0;
	for (int i=0;i<pnum;i++){
		float av_len_p=0.0;
		// get k-nearest neighbors
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		kdtree.nearestKSearch (cloud->points[i], K+1, pointIdxNKNSearch, pointNKNSquaredDistance); //排除当前点自己
		for (int j=0;j<(int)pointNKNSquaredDistance.size();j++)
		{
			if (pointIdxNKNSearch[j]==i)
			{
				continue;
			}
			av_len_p+=pointNKNSquaredDistance[j];
		}
		av_len_p/=K;
		av_len+=av_len_p;
	}
	av_len/=pnum;
	cout<<"average edge len:"<<av_len<<endl;
	//compute normal
	std::vector<pcl::Normal> normals(pnum); 
	K = 50;
	for (int i=0;i<pnum;i++)
	{
		std::vector<int> inliers;
		// get k-nearest neighbors
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		kdtree.nearestKSearch (cloud->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance);
		pcl::PointCloud<pcl::PointXYZ>::Ptr neighbors(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::copyPointCloud(*cloud,pointIdxNKNSearch,*neighbors);
		// created RandomSampleConsensus object and compute the appropriated model
		pcl::SampleConsensusModelPlane<pcl::PointXYZ>::Ptr
			model_p (new pcl::SampleConsensusModelPlane<pcl::PointXYZ> (neighbors));
		pcl::RandomSampleConsensus<pcl::PointXYZ> ransac (model_p);
		ransac.setDistanceThreshold (av_len);
		ransac.computeModel();
		ransac.getInliers(inliers);
		//compute normal for point i
		pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
		ne.setInputCloud(neighbors);
		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
		ne.setSearchMethod (tree);
		float nx,ny,nz,curvature;
		ne.computePointNormal(*neighbors,inliers,nx,ny,nz,curvature);
		normals[i]=pcl::Normal(nx,ny,nz);
	}
	
	std::string name=filename.substr(0,filename.find_last_of("."));
	cout<<"filename:"<<filename<<",name:"<<name<<endl;

	std::string out=name+".xyzn";
	ofstream ofh(out.c_str());
	for (int i=0;i<(int)cloud->points.size();i++)
	{
		pcl::PointXYZ& p=cloud->points[i];
		pcl::Normal& n=normals[i];
		ofh<<p.x<<" "<<p.y<<" "<<p.z<<" "<<n.normal_x<<" "<<n.normal_y<<" "<<n.normal_z<<endl;
	}
	ofh.close();

	return 0;
}

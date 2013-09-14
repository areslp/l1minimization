% mex z_loop.cpp COMPFLAGS="/openmp $COMPFLAGS"
mex -I"C:\Program Files (x86)\PCL 1.6.0\include\pcl-1.6" -I"D:\lib\boost_1_53_0" -I"C:\Program Files (x86)\Eigen\include" -I"C:\Program Files (x86)\flann\include" -L"C:\Program Files (x86)\PCL 1.6.0\lib" -lpcl_common_release -lpcl_kdtree_release -lpcl_search_release -L"C:\Program Files (x86)\flann\lib" -lflann_cpp_s -lflann compute_WL.cpp COMPFLAGS="/openmp $COMPFLAGS"


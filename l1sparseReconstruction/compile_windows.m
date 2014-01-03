clear mex;

mex z_loop.cpp COMPFLAGS="/openmp $COMPFLAGS"

mex -I"C:\Program Files\PCL 1.6.0\include\pcl-1.6" -I"C:\Program Files\PCL 1.6.0\3rdParty\Boost\include" -I"C:\Program Files\PCL 1.6.0\3rdParty\Eigen\include" -I"C:\Program Files\PCL 1.6.0\3rdParty\FLANN\include" -L"C:\Program Files\PCL 1.6.0\3rdParty\Boost\lib" -L"C:\Program Files\PCL 1.6.0\lib" -lpcl_common_release -lpcl_kdtree_release -lpcl_search_release -L"C:\Program Files\PCL 1.6.0\3rdParty\FLANN\lib" -lflann_cpp_s -lflann compute_WL.cpp COMPFLAGS="/openmp $COMPFLAGS"

mex -I"C:\Program Files\PCL 1.6.0\include\pcl-1.6" -I"C:\Program Files\PCL 1.6.0\3rdParty\Boost\include" -I"C:\Program Files\PCL 1.6.0\3rdParty\Eigen\include" -I"C:\Program Files\PCL 1.6.0\3rdParty\FLANN\include" -L"C:\Program Files\PCL 1.6.0\3rdParty\Boost\lib" -L"C:\Program Files\PCL 1.6.0\lib" -lpcl_common_release -lpcl_kdtree_release -lpcl_search_release -L"C:\Program Files\PCL 1.6.0\3rdParty\FLANN\lib" -lflann_cpp_s -lflann compute_AE.cpp COMPFLAGS="/openmp $COMPFLAGS"

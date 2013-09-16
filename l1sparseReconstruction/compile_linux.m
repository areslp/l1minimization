mex z_loop.cpp CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" COMPFLAGS="/openmp $COMPFLAGS"

mex -I"/usr/local/include/pcl-1.7" -I"/usr/local/include/boost" -I"/usr/local/include/eigen3" -I"/usr/include/flann" -L"/usr/local/lib" -lpcl_common -lpcl_kdtree -lpcl_search -lflann_cpp_s -lflann compute_WL.cpp CFLAGS="\$CFLAGS -fopenmp -O3" LDFLAGS="\$LDFLAGS -fopenmp" COMPFLAGS="/openmp $COMPFLAGS"

mex -I"/usr/local/include/pcl-1.7" -I"/usr/local/include/boost" -I"/usr/local/include/eigen3" -I"/usr/include/flann" -L"/usr/local/lib" -lpcl_common -lpcl_kdtree -lpcl_search -lflann_cpp_s -lflann compute_AE.cpp CFLAGS="\$CFLAGS -fopenmp -O3" LDFLAGS="\$LDFLAGS -fopenmp" COMPFLAGS="/openmp $COMPFLAGS"


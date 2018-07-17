README

Default Algorithm: Primal Dual

Compilation:
	1) cd build
	2) cmake ..
	3) make

Run:
	cd build
	./StereoMatching

Result: The depth will be stored in build folder as 


To use brute force algorithm:
	1) Please uncomment/enable line number 38 and 39
	2) comment/disable line number 42-56
	3) Compile and run as explained above


Libraries used:
	1) openexr
	2) ilmbase
	3) FreeImage
	4) Eigen
	5) tbb (thread building blocks)

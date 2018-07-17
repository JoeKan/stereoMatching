README  

Compilation:   
	1) cd build     
	2) cmake ..   
	3) make   

Run:    
	cd build
	
	For Primal Dual Algorithms run following:      
	./StereoMatching -m p   
	For Brute Force Algorithm run following:   
	./StereoMatching -m b    

Result:    
1) Primal Dual Algorithm: The depth image will be stored in build folder as "out_primaldual.png"      
2) Brute Force Algorithm: The depth image will be stored in build folder as "out_BruteForce.png"    


Libraries used:    
	1) openexr   
	2) ilmbase    
	3) FreeImage    
	4) Eigen    
	5) tbb (thread building blocks)   

Result:  
Sample result are also uploaded in the folder "Result"    

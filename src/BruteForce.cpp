#include "parameters.h"
#include "utility.h"
#include "Eigen.h"
#include "LoaderEXR.h"
#include "FreeImageHelper.h" 
#include "parameters.h"
#include "utility.h"


void brute_force_depth_calc(BYTE* colorFrame_r) {
    
    /*_d is the cost volume matrix */
    Eigen::MatrixXf _d(pixels, d_range);
    Eigen::VectorXf _d_Min(pixels);
    for(uint i = 0; i < pixels; ++i){
        for(uint j=0; j< d_range; j++)
            _d(i,j) = 0;
    }

   
    //std::cout << "inverse" << std::endl ;
    //std::cout << KInv << std::endl ;

    //BYTE* colorFrame_r = new BYTE[4* 640*480];
    //processNextFrame(8, colorFrame_r);
    BYTE* colorFrame_m = new BYTE[4* 640*480];

    for(uint frameNum = 1; frameNum<=m; ++frameNum){
        std::cout << "START: frameNum " << frameNum <<std::endl;
        processNextFrame(frameNum, colorFrame_m);
		// parallel version
		tbb::parallel_for(0, pixels, [&](int pixel){
			float step_depth = init_depth;
			for (uint d = 0; d<d_range; d++) {
				Eigen::Vector3f I_r = {
					(float)colorFrame_r[(pixel * 4)],
					(float)colorFrame_r[(pixel * 4) + 1],
					(float)colorFrame_r[(pixel * 4) + 2] };
				_d(pixel, d) += rho_r(pixel, step_depth, colorFrame_r, colorFrame_m, I_r, frameNum, 8);
				step_depth += inc_depth;
			}
		});
        std::cout << "END: frameNum " << frameNum << std::endl;
    }

    std::cout << "DONE" << std::endl;
	tbb::parallel_for(0, pixels, [&](int idx) {
		int index = 0;
		_d.row(idx).minCoeff(&index);
		_d_Min(idx) = init_depth + index * inc_depth;
	});

   std::cout << "_d_Min size = " << _d_Min.size() << std::endl;
   // now we need to visualize this _d_Min and comapre it with 0008.exr file
    std::cout << "Going to read exr file" << std::endl;
    std::vector <float> image_vir(640*480);
    std::string vir_filename = DATA_SYNTHETIC_DIR + "0008.exr";
    read_openexr(vir_filename,image_vir.data(), 640, 480, 1);
    /*for (std::vector<float>::const_iterator i = image_vir.begin(); i != image_vir.end(); ++i) {
        std::cout << *i << ' ';
    }*/
    std::cout << std::endl;

    //Eigen::VectorXf image_vir_Eigen(image_vir.data(), image_vir.size());
    Eigen::Map<Eigen::VectorXf> image_vir_Eigen(image_vir.data(),image_vir.size());
    std::cout << "image_vir_Eigen Size= " << image_vir_Eigen.size() << std::endl;
    Eigen::VectorXf Residue = _d_Min - image_vir_Eigen;
    std::cout << "Final L2 norm b/w exr file and calculated depth = " << Residue.lpNorm<2>() << std::endl;
	FreeImageB outImage(640, 480, 3);
	BYTE* outData = new BYTE[640 * 480 * 3];
	tbb::parallel_for(0, pixels, [&](int idx) {
        // 255 = white, 0 = black
        // close is white and far is black
        outData[idx * 3] = 255 - (_d_Min[idx] / max_depth) * 255;
        outData[idx * 3 + 1] = 255 - (_d_Min[idx] / max_depth) * 255;
        outData[idx * 3 + 2] = 255 - (_d_Min[idx] / max_depth) * 255;
	});
	outImage.data = outData;
	outImage.SaveImageToFile("out.png");

    delete colorFrame_r;
    
    return;
}
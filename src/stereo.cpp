#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <tbb/parallel_for.h>

#include "Eigen.h"
#include "LoaderEXR.h"
#include "FreeImageHelper.h" 
#include "parameters.h"
#include "utility.h"

Eigen::Vector4f pi_inverse(int pixelNum, float depth);
Eigen::Matrix4f T_mr(int m_frame, int r_frame);
Eigen::Vector2f pi(Eigen::Vector3f vec);

bool processNextFrame(int filenum, BYTE* colorFrame);
void brute_force_depth_calc(BYTE* colorFrame_r);
std::vector<Eigen::Matrix4f> poses;

Eigen::Matrix3f K;
Eigen::Matrix3f KInv;

#define uint unsigned int

int main(){
    BYTE* colorFrame_r = new BYTE[4* 640*480];
    
    K << 577.97f, 0.0f, 320.0f,
		 0.0f, 575.44f, 240.0f,
		 0.0f, 0.0f, 1.0f;

    KInv = K.inverse();
    load_all_matrices_from_n_files(poses);
    
    processNextFrame(8, colorFrame_r);
    brute_force_depth_calc(colorFrame_r);
 
    return 0;
}
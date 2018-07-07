#include "Eigen.h"
#include "parameters.h"
#include "utility.h"

bool processNextFrame(int filenum, BYTE* colorFrame);
void brute_force_depth_calc(BYTE* colorFrame_r);
std::vector<Eigen::Matrix4f> poses;

Eigen::Matrix3f K;
Eigen::Matrix3f KInv;

#define uint unsigned int

int main(){
    int current_ref_img = 8;
    BYTE* colorFrame_r = new BYTE[4* 640*480];
    read_images(current_ref_img, colorFrame_r);
    
    BYTE** colorFrames_b = new BYTE*[m];
    for(uint i=0;i<m;++i){
        colorFrames_b[i] = new BYTE[4* 640*480];
    }

    for(uint i=0;i<m;++i){
        read_images(current_ref_img - i - 1, colorFrames_b[i]);
    }
    
    K << 577.97f, 0.0f, 320.0f,
		 0.0f, 575.44f, 240.0f,
		 0.0f, 0.0f, 1.0f;

    KInv = K.inverse();
    load_all_matrices_from_n_files(poses);
    
    processNextFrame(8, colorFrame_r);
    //brute_force_depth_calc(colorFrame_r);

    //------Primal_Dual-------------
    Eigen::VectorXf d_(pixels);
    PrimalDual(current_ref_img, colorFrame_r, colorFrames_b, d_);

 
    return 0;
}
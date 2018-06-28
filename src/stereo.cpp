#include "FreeImageHelper.h" 
#include <string>
#include <iostream>
#include <cstring>
#include "Eigen.h"
#include <fstream>
#include <algorithm>
#include <tbb/parallel_for.h>
#include "LoaderEXR.h"


bool processNextFrame(int filenum, BYTE* colorFrame);
//bool ReadPoseFile(int filenum, Eigen::MatrixXf &Pose);
Eigen::Vector4f pi_inverse(int pixelNum, float depth);
void load_all_matrices_from_n_files(std::vector<Eigen::Matrix4f> &P);
Eigen::Matrix4f T_mr(int m_frame, int r_frame);
Eigen::Vector2f pi(Eigen::Vector3f vec);
float rho_r(int pixel, float depth, BYTE* colorFrame_r, BYTE* colorFrame_m, Eigen::Vector3f &I_r, int m_frame, int r_frame);

Eigen::Matrix3f K;
Eigen::Matrix3f KInv;

std::string DATA_SYNTHETIC_DIR;
const int pixelsWidth = 640;
const int pixelsHeight = 480;
const int pixels = pixelsWidth * pixelsHeight;
const int d_range = 30;
const int m = 7;

std::vector<Eigen::Matrix4f> poses;

#define uint unsigned int

int main(){
    
    DATA_SYNTHETIC_DIR = "./../output640x480/";
    load_all_matrices_from_n_files(poses);

    Eigen::MatrixXf _d(pixels, d_range);
    Eigen::VectorXf _d_Min(pixels);
    for(uint i = 0; i < pixels; ++i){
        for(uint j=0; j< d_range; j++)
            _d(i,j) = 0;
    }

    /*K <<    525.0f, 0.0f, 319.5f,
		    0.0f, 525.0f, 239.5f,
		    0.0f, 0.0f, 1.0f;*/
    K <<    577.97f, 0.0f, 320.0f,
		    0.0f, 575.44f, 240.0f,
		    0.0f, 0.0f, 1.0f;
    KInv = K.inverse();
    //std::cout << "inverse" << std::endl ;
    //std::cout << KInv << std::endl ;

    BYTE* colorFrame_r = new BYTE[4* 640*480];
    processNextFrame(8, colorFrame_r);
    BYTE* colorFrame_m = new BYTE[4* 640*480];

    for(uint frameNum = 1; frameNum<=m; ++frameNum){
        std::cout << "START: frameNum " << frameNum <<std::endl;
        processNextFrame(frameNum, colorFrame_m);
		// parallel version
		tbb::parallel_for(0, pixels, [&](int pixel){
			float step_depth = 0.1;
			for (uint d = 0; d<d_range; d++) {
				Eigen::Vector3f I_r = {
					(float)colorFrame_r[(pixel * 4)],
					(float)colorFrame_r[(pixel * 4) + 1],
					(float)colorFrame_r[(pixel * 4) + 2] };
				_d(pixel, d) += rho_r(pixel, step_depth, colorFrame_r, colorFrame_m, I_r, frameNum, 8);
				step_depth += 0.1;
			}
		});
        /*for(uint pixel = 0; pixel< pixels; pixel++){
			float step_depth = 0.1;
			for (uint d = 0; d<d_range; d++) {
				Vector3f I_r = {
					(float)colorFrame_r[(pixel * 4)],
					(float)colorFrame_r[(pixel * 4) + 1],
					(float)colorFrame_r[(pixel * 4) + 2] };
				_d(pixel, d) += rho_r(pixel, step_depth, colorFrame_r, colorFrame_m, I_r, frameNum, 8);
				step_depth += 0.1;
			}
        }*/
        std::cout << "END: frameNum " << frameNum << std::endl;
    }
    //-------------TEST CODE-------
    /*Eigen::MatrixXf mat(2,4);
    mat << 7.01, 2.01, 6.01, 9.01,
           3.01, 11.01, 7.01, 8.01;
    std::cout << mat.rowwise().minCoeff() << std::endl;*/

    //-----------------------
    std::cout << "DONE" << std::endl;
    // _d why everything inside _d is int ??
    //std::cout << _d << std::endl;
    std::cout << "_d size = " << _d.rowwise().minCoeff().size() << std::endl;
   _d_Min = _d.rowwise().minCoeff();
   std::cout << "_d_Min size = " << _d_Min.size() << std::endl;
   // now we need to visualize this _d_Min and comapre it with 0008.exr file

    std::cout << "Going to read exr file" << std::endl;
    std::vector <float> image_vir(640*480);
    std::string vir_filename = "/media/virendra/data/study/4_sem/3D_Scan/Project/stereoMatching/output640x480/0008.exr";
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


	// bruteforce depth map
	float* inverseDepth = new float[640 * 480];
	for (uint i = 0; i < pixels; ++i) {
		float min = 255.0f;
		uint d = 0;
		for (uint d = 0; d < d_range; ++d) {
			if (_d(i, d) < min) {
				min = _d(i, d);
			}
		}
		inverseDepth[i] = d;
	}

    ////manually check color values
    //for(unsigned int j = 0; j < 3; j++){
    //    std::cout<<I_r(j)<<std::endl;
    //}

    //std::cout<<pi_inverse(0, 0.1)<<std::endl;
    //std::cout<<pi_inverse(640,0.2)<<std::endl;

    //std::cout<<colorFrame;

    delete colorFrame_r;

return 0;
}

float rho_r(int pixel, float depth, BYTE* colorFrame_r, BYTE* colorFrame_m, Eigen::Vector3f &I_r, int m_frame, int r_frame){
    Eigen::Vector4f pi_inv = pi_inverse(pixel, depth);
    Eigen::Matrix4f temp = T_mr(m_frame, r_frame);
    Eigen::MatrixXf t_mr(3,4);
    t_mr = temp.block<3,4>(0,0);
    Eigen::Vector2f m_coordinate_f = pi(K * t_mr * pi_inv );
	// Clamp coordinates to actual image space, TODO: maybe skip if coordinates are out of range?
	const int coordX = std::max(0, std::min(pixelsWidth, (int)m_coordinate_f.x()));
	const int coordY = std::max(0, std::min(pixelsHeight, (int)m_coordinate_f.y()));
    //std::cout << "hello" << std::endl;
    Eigen::Vector3f I_m(0.0,0.0,0.0);
    int index_in_img =  coordY*pixelsWidth + coordX;
    //std::cout << m_coordinate_f <<std::endl;
    for(unsigned int j = 0; j < 3; j++){ //copying RGB values
        I_m(j) = colorFrame_m[(index_in_img* 4)+j];
    }

    Eigen::Vector3f Diff = I_r - I_m;
    //std::cout << "herere = " << Diff.lpNorm<1>() << std::endl;
    // Question: Won't L1 norm will make values very large?
    return Diff.lpNorm<1>();
}

Eigen::Vector4f pi_inverse(int pixelNum, float depth){
    Eigen::Vector3f u_dot = Eigen::Vector3f(pixelNum % 640,(int)(pixelNum/ 640), 1.0);
    //std::cout<<u_dot<<std::endl;
    Eigen::Vector3f tmp= (KInv * u_dot) * 1/depth;

    return Eigen::Vector4f(tmp.x(),tmp.y(),tmp.z(),1.0);
}

Eigen::Vector2f pi(Eigen::Vector3f vec){
    return Eigen::Vector2f(vec.x()/vec.z(), vec.y()/vec.z());
}

Eigen::Matrix4f T_mr(int m_frame, int r_frame){
    Eigen::Matrix4f temp_mw;
    temp_mw.setIdentity();
    temp_mw.block<3,3>(0,0) = poses[m_frame].block<3,3>(0,0).transpose();
    temp_mw.block<3,1>(0,3) = -1 * poses[m_frame].block<3,1>(0,3);
    return temp_mw*poses[r_frame];
}

bool processNextFrame(int filenum, BYTE* colorFrame){
    FreeImageB rgbImg;
    std::string filename;

    if((int)(filenum / 10) == 0){
        filename = "000" + std::to_string(filenum);
    }
    else if((int)(filenum / 100) == 0){
        filename = "00" + std::to_string(filenum);
    }
    else if(filenum == 0){
        filename = "0000";
    }
    else {
        filename = "0" + std::to_string(filenum);
    }

    rgbImg.LoadImageFromFile(DATA_SYNTHETIC_DIR + filename +".jpg");

    memcpy(colorFrame, rgbImg.data, 4 * 640 * 480);

    std::cout<<"./output640x480/" + filename<<std::endl;

    return true;
}

void load_all_matrices_from_n_files(std::vector<Eigen::Matrix4f> &P) {
    const int n_images = 120; /*N_FRAMES*100;*/
    P.resize(n_images);
    for (int i = 0; i < n_images; i++) {
        P[i].setIdentity();
        char buf[20];
        snprintf(buf, sizeof buf, "/%04d", i+ 1);
        std::string filename0 = std::string(DATA_SYNTHETIC_DIR) + std::string(buf) + ".dat";
        std::ifstream file(filename0);
        Eigen::Matrix4d Ptmp = Eigen::Matrix4d::Identity();
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 3; k++) {
                file >> Ptmp(k, j);
            }
        }
        //P[i] = Ptmp.inverse().eval().cast<float>();
        Ptmp.col(0).normalize();
        Ptmp.col(1).normalize();
        Ptmp.col(2).normalize();
        P[i] = Ptmp.inverse().eval().cast<float>();
    }
}

/*
bool ReadPoseFile(int filenum, Eigen::MatrixXf &Pose)
	{
        std::string filename;
        if((int)(filenum / 10) == 0){
            filename = "000" + std::to_string(filenum);
        }
        else if((int)(filenum / 100) == 0){
            filename = "00" + std::to_string(filenum);
        }
        else if(filenum == 0){
            filename = "0000";
        }
        else {
            filename = "0" + std::to_string(filenum);
        }

        filename = "./output640x480/" + filename + ".dat";

		std::ifstream file(filename, std::ios::in);
		if (!file.is_open()) return false;

		while (file.good())
		{
            file>>  Pose(0,0) >> Pose(0,1) >> Pose(0,2) >> Pose(0,3) >>
                    Pose(1,0) >> Pose(1,1) >> Pose(1,2) >> Pose(1,3) >>
                    Pose(2,0) >> Pose(2,1) >> Pose(2,2) >> Pose(2,3);
		}
		file.close();

		return true;
	}
*/

#include "FreeImageHelper.h" 
#include <string>
#include <iostream>
#include <cstring>
#include "Eigen.h"
#include <fstream>


bool processNextFrame(int filenum, BYTE* colorFrame);
//bool ReadPoseFile(int filenum, Eigen::MatrixXf &Pose);
Eigen::Vector4f pi_inverse(int pixelNum, float depth);
void load_all_matrices_from_n_files(std::vector<Eigen::Matrix4f> &P);
Eigen::Matrix4f T_mr(int m_frame, int r_frame);
Eigen::Vector2f pi(Vector3f vec);
float rho_r(int pixel, float depth, BYTE* colorFrame_r, BYTE* colorFrame_m, Vector3f &I_r, int m_frame, int r_frame);

Eigen::Matrix3f K;
Eigen::Matrix3f KInv;

std::string DATA_SYNTHETIC_DIR;
const int pixels = 640*480;
const int d_range = 30;
const int m = 7;

std::vector<Eigen::Matrix4f> poses;

int main(){

    DATA_SYNTHETIC_DIR = "./../output640x480/";
    load_all_matrices_from_n_files(poses);

    Eigen::MatrixXf _d(pixels, d_range);
    for(uint i = 0; i < pixels; ++i){
        for(uint j=0; j< d_range; j++)
            _d(i,j) = 0;
    }

    K <<    525.0f, 0.0f, 319.5f,
		    0.0f, 525.0f, 239.5f,
		    0.0f, 0.0f, 1.0f;
    KInv = K.inverse();
    std::cout << "inverse" << std::endl ;
    std::cout << KInv << std::endl ;

    BYTE* colorFrame_r = new BYTE[4* 640*480];
    processNextFrame(8, colorFrame_r);
    BYTE* colorFrame_m = new BYTE[4* 640*480];

    Vector3f I_r(0.0,0.0,0.0);
    for(uint j = 1; j<=m; ++j){
        processNextFrame(j, colorFrame_m);
        for(uint pixel = 0; pixel< pixels; pixel++){
            float step_depth = 0.1;
            for(uint d = 0; d<d_range; d++){
                for(unsigned int j = 0; j < 3; j++){
                    I_r(j) = colorFrame_r[(pixel * 4)+j];
                }
                _d(pixel, d) += rho_r(pixel, step_depth, colorFrame_r, colorFrame_m, I_r, j, 8);
                step_depth += 0.1;
            }
        }
    }

    //manually check color values
    for(unsigned int j = 0; j < 3; j++){
        std::cout<<I_r(j)<<std::endl;
    }

    std::cout<<pi_inverse(0, 0.1)<<std::endl;
    std::cout<<pi_inverse(640,0.2)<<std::endl;

    //std::cout<<colorFrame;

    delete colorFrame_r;

return 0;
}

float rho_r(int pixel, float depth, BYTE* colorFrame_r, BYTE* colorFrame_m, Vector3f &I_r, int m_frame, int r_frame){
    Eigen::Vector4f pi_inv = pi_inverse(pixel, depth);
    Eigen::Matrix4f temp = T_mr(m_frame, r_frame);
    Eigen::MatrixXf t_mr(3,4);
    t_mr = temp.block<3,4>(0,0);
    Eigen::Vector2f m_coordinate_f = pi(K * t_mr * pi_inv ) ;
    std::cout << "hello" << std::endl;
    Eigen::Vector3f I_m(0.0,0.0,0.0);
    int index_in_img =  (int)(m_coordinate_f.y())*640 + (int)(m_coordinate_f.x());
    std::cout << m_coordinate_f <<std::endl;
    for(unsigned int j = 0; j < 3; j++){ //copying RGB values
        I_m(j) = colorFrame_m[(index_in_img* 4)+j];
    }

    Eigen::Vector3f Diff = I_r - I_m;
    std::cout << "herere" << std::endl;
    return Diff.lpNorm<1>();
}

Eigen::Vector4f pi_inverse(int pixelNum, float depth){
    Eigen::Vector3f u_dot = Vector3f(pixelNum % 640,(int)(pixelNum/ 640), 1.0);
    std::cout<<u_dot<<std::endl;
    Eigen::Vector3f tmp= (KInv * u_dot) * 1/depth;

    return Eigen::Vector4f(tmp.x(),tmp.y(),tmp.z(),1.0);
}

Eigen::Vector2f pi(Vector3f vec){
    return Vector2f(vec.x()/vec.z(), vec.y()/vec.z());
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

    rgbImg.LoadImageFromFile("./output640x480/" + filename +".jpg");

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

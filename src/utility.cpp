/* Please refer the comment written in utility.h file once.
*/

#include "utility.h"
#include "parameters.h"
#include "LoaderEXR.h"
#include "FreeImageHelper.h" 

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

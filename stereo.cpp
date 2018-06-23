#include "FreeImageHelper.h" 
#include <string>
#include <iostream>
#include <cstring>
#include "Eigen.h"
#include <fstream>


bool processNextFrame(int filenum, BYTE* colorFrame);
bool ReadPoseFile(int filenum, Eigen::MatrixXf &Pose);
Vector3f pi_inverse(int pixelNum, float depth);

Eigen::Matrix3f K;
Eigen::Matrix3f KInv;

int main(){

    float _d[30];
    for(uint i = 0;i<30;++i){
        _d[i] = MINF;
    }

    K <<    525.0f, 0.0f, 319.5f,
		0.0f, 525.0f, 239.5f,
		0.0f, 0.0f, 1.0f;
    KInv = K.inverse();
    

    BYTE* colorFrame = new BYTE[4* 640*480];
    Eigen::MatrixXf Pose(3,4);
    processNextFrame(7, colorFrame);

    ReadPoseFile(7, Pose);
    //processNextFrame(29);
    //processNextFrame(120);

    std::cout<<Pose<<std::endl;


    Vector3f I_r(0.0,0.0,0.0);

    for(uint i = 0; i< 640*480;i++){
        for(unsigned int j = 0; j < 3; j++){
            I_r(j) = colorFrame[(i * 4)+j];
        }
    }

    //manually check color values
    for(unsigned int j = 0; j < 3; j++){
        std::cout<<I_r(j)<<std::endl;
    }

    std::cout<<pi_inverse(0, 0.1)<<std::endl;
    std::cout<<pi_inverse(640,0.2)<<std::endl;

    //std::cout<<colorFrame;

    delete colorFrame;

return 0;
}

Vector3f pi_inverse(int pixelNum, float depth){
    Vector3f u_dot = Vector3f((int)(pixelNum/ 640), pixelNum % 640, 1.0);
    //std::cout<<u_dot<<std::endl;
    return (KInv * u_dot) * 1/depth;
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
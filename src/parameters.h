/* All the constant terms should be declared here.
*/

#if !defined(STEREO_PARAMETERS_H)
#define STEREO_PARAMETERS_H 1

#include <iostream>
#include "Eigen.h"

const int pixelsWidth = 640;
const int pixelsHeight = 480;
const int pixels = pixelsWidth * pixelsHeight;
const int d_range = 30;
const int m = 7;

const std::string DATA_SYNTHETIC_DIR = "./../output640x480/";

extern Eigen::Matrix3f K;
extern Eigen::Matrix3f KInv;
extern std::vector<Eigen::Matrix4f> poses;

#endif
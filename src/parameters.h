/* All the constant terms should be declared here.
*/

#if !defined(STEREO_PARAMETERS_H)
#define STEREO_PARAMETERS_H 1

#include <iostream>
#include "Eigen.h"

const int pixelsWidth = 640;
const int pixelsHeight = 480;
const int pixels = pixelsWidth * pixelsHeight;
const int d_range = 100;
const int a_range = 30;
// const float init_depth = 0.5f;
const float inc_a = 0.1f;
const float min_depth = 0.5;
const float max_depth = 8/*init_depth + (d_range - 1) * inc_depth*/;
const float inc_depth = (max_depth - min_depth)/d_range;
const int m = 7;

const std::string DATA_SYNTHETIC_DIR = "./../output640x480/";

extern Eigen::Matrix4f K;
extern Eigen::Matrix4f KInv;
extern std::vector<Eigen::Matrix4f> poses;

#endif
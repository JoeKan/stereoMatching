/* This file contains code which will be used by many different files.
   For example: calculating error term: rho_r or calculating pi() etc.
   Please include utility.h in your cpp file for using these functions.
*/

#if !defined(STEREO_UTILITY_H)
#define STEREO_UTILITY_H 1

#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <tbb/parallel_for.h>
#include <FreeImage.h>

#include "Eigen.h"
#include "utility.h"

#define uint unsigned int

float rho_r(int pixel, float depth, BYTE* colorFrame_r, BYTE* colorFrame_m, Eigen::Vector3f &I_r, int m_frame, int r_frame);
Eigen::Vector4f pi_inverse(int pixelNum, float depth);
Eigen::Vector2f pi(Eigen::Vector3f vec);
Eigen::Matrix4f T_mr(int m_frame, int r_frame);
bool processNextFrame(int filenum, BYTE* colorFrame);
void load_all_matrices_from_n_files(std::vector<Eigen::Matrix4f> &P);
void read_images(int num, BYTE* colorFrame);
void cost_calc(BYTE* colorFrame_r, BYTE** colorFrames_b, int current_ref_img, Eigen::MatrixXf& _a);
void PrimalDual(int current_ref_img, BYTE* colorFrame_r, BYTE** colorFrames_b, Eigen::VectorXf& d_);

#endif
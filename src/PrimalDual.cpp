/* In this file we will implement the primal dual algorithms
   Please refer section 2.2.3 from dtam paper.
*/

#include <FreeImage.h>
#include "Eigen.h"
#include <cmath>
#include "parameters.h"
#include "utility.h"

// 6.2.3 of A. Chambolle and T. Pock. A first-order primal-dual
// algorithm for convex problems with applications to imaging.
float sigma_d(float epsilon, float theta) {

	float L = 2.0;

	float mu = 2.0*std::sqrt(epsilon / theta) / L;

	return mu / (2.0 / theta);
}
float sigma_q(float epsilon, float theta) {

	float L = 2.0;

	float mu = 2.0*std::sqrt(epsilon / theta) / L;

	return mu / (2.0*epsilon);
}
void computeG(float* g, BYTE* img, int width, int height, float alphaG, float betaG) {
	// equation 5 of the paper
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int idx = y * width + x;
			// forward difference gradient
			float gx = (x == width - 1) ? 0.0f : img[idx + 1] / 255.0f - img[idx] / 255.0f;
			float gy = (y == height - 1) ? 0.0f : img[idx + width] / 255.0f - img[idx] / 255.0f;
			g[idx] = expf(-alphaG * powf(sqrtf(gx*gx + gy * gy), betaG));
		}
	}
}
void updateQ(float* g, Eigen::VectorXf& a, Eigen::VectorXf& q, Eigen::VectorXf& d, int width, int height, float sigma_q, float sigma_d, float epsilon, float theta) {
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int idx = y * width + x;
			// gradient with forward difference
			// TODO use closed form
			float dd_x = (x == width - 1) ? 0.0f : d[idx + 1] - d[idx];
			float dd_y = (y == height - 1) ? 0.0f : d[idx + width] - d[idx];

			float qx = (q[idx] + sigma_q * g[idx] * dd_x) / (1.0f + sigma_q * epsilon);
			float qy = (q[idx + width * height] + sigma_q * g[idx] * dd_y) / (1.0f + sigma_q * epsilon);

			float maxq = std::fmaxf(1.0f, sqrtf(qx*qx + qy * qy));
			q[idx] = qx / maxq;
			q[idx + width * height] = qy / maxq;
		}
	}
}
void updateD(float* g, Eigen::VectorXf& a, Eigen::VectorXf& q, Eigen::VectorXf& d, int width, int height, float sigma_q, float sigma_d, float epsilon, float theta) {
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int idx = y * width + x;
			// gradient with backward difference
			// TODO use closed form
			float dqx_x = (x == 0) ? q[idx] - q[idx + 1] : q[idx] - q[idx - 1];
			float dqy_y = (y == 0) ? q[idx + width * height] - q[idx + width * height + width] : q[idx + width * height] - q[idx + width * height - width];
			float div_q = dqx_x + dqy_y;

			d[idx] = (d[idx] + sigma_d * (g[idx] * div_q + a[idx] / theta)) / (1.0f + sigma_d / theta);
		}
	}
}
void PrimalDual(int current_ref_img, BYTE* colorFrame_r, BYTE** colorFrames_b, Eigen::VectorXf& d_) {
    /* TODO */

	//_Eaux_Min = a_u for all pixels as written in the paper
	// d_ is d_u for all pixels per frame

	using namespace Eigen;
	float alphaG = 100.0f;
	float betaG = 1.6f;
	float theta_start = 0.2;
	float theta_min = 1.0e-4;
	float theta_step = 0.97;
	float epsilon = 0.00147;
	float lambda = 0.80;
	float theta = theta_start;

	float* g = new float[640 * 480];
	computeG(g, colorFrame_r, 640, 480, alphaG, betaG);
	Eigen::VectorXf q(640,480);

	Eigen::MatrixXf C_u_a(pixels, a_range);
	for(uint i = 0; i < pixels; ++i){
        for(uint j=0; j< a_range; j++)
            C_u_a(i,j) = 0;
    }

	cost_calc(colorFrame_r, colorFrames_b, current_ref_img, C_u_a);

	//initializing "_Eaux_Min" vector by doing argmin(cost ie C_a)
	Eigen::VectorXf _Eaux_Min(pixels);
	tbb::parallel_for(0, pixels, [&](int idx) {
		int index = 0;
		C_u_a.row(idx).minCoeff(&index);
		_Eaux_Min(idx) = init_depth + index * inc_a;
	});

	//first time initialization of d_u = a_u. check if this assignment makes a deep copy.
	d_ = _Eaux_Min;

	updateQ(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);
	updateD(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);

	
	int n = 1;
	Eigen::MatrixXf _Eaux(pixels, a_range);
	while (theta > theta_min) {
		create_Eaux(d_, _Eaux_Min, theta, lambda, C_u_a, _Eaux);

		tbb::parallel_for(0, pixels, [&](int idx) {
			int index = 0;
			_Eaux.row(idx).minCoeff(&index);
			_Eaux_Min(idx) = init_depth + index * inc_a;
		});

		updateQ(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);
		updateD(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);

		// TODO: update a - see equations 13 and 14 of the paper.

		float beta = (theta > 1e-3) ? 1e-3 : 1e-4;
		theta *= (1 - beta * n);
		n++;
	}
}

void create_Eaux(Eigen::VectorXf& d_, Eigen::VectorXf& a_, float theta, float lambda, Eigen::MatrixXf& C_u_a, Eigen::MatrixXf& _Eaux){
	tbb::parallel_for(0, pixels, [&](int pixel){
		for (uint i = 0; i < a_range; i++) {
			float E = (1/theta * (d_[pixel] - a_[pixel]) * (d_[pixel] - a_[pixel])) + lambda * C_u_a(pixel, i);
			_Eaux(pixel, i) = E;
		}
	});
}
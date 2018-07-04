/* In this file we will implement the primal dual algorithms
   Please refer section 2.2.3 from dtam paper.
*/

#include <FreeImage.h>
#include "Eigen.h"
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
			g[idx] = std::expf(-alphaG * std::powf(std::sqrtf(gx*gx + gy * gy), betaG));
		}
	}
}
void updateQ(float* g, float* a, float* q, float* d, int width, int height, float sigma_q, float sigma_d, float epsilon, float theta) {
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int idx = y * width + x;
			// gradient with forward difference
			// TODO use closed form
			float dd_x = (x == width - 1) ? 0.0f : d[idx + 1] - d[idx];
			float dd_y = (y == height - 1) ? 0.0f : d[idx + width] - d[idx];

			float qx = (q[idx] + sigma_q * g[idx] * dd_x) / (1.0f + sigma_q * epsilon);
			float qy = (q[idx + width * height] + sigma_q * g[idx] * dd_y) / (1.0f + sigma_q * epsilon);

			float maxq = std::fmaxf(1.0f, std::sqrtf(qx*qx + qy * qy));
			q[idx] = qx / maxq;
			q[idx + width * height] = qy / maxq;
		}
	}
}
void updateD(float* g, float* a, float* q, float* d, int width, int height, float sigma_q, float sigma_d, float epsilon, float theta) {
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
void PrimalDual(BYTE* colorFrame_r, Eigen::MatrixXf& d, Eigen::VectorXf& d_min) {
    /* TODO */
	using namespace Eigen;
	float alphaG = 100.0f;
	float betaG = 1.6f;
	float theta_start = 0.2;
	float theta_min = 1.0e-4;
	float theta_step = 0.97;
	float epsilon = 0.00147;
	float lambda = 0.80;

	float* g = new float[640 * 480];
	computeG(g, colorFrame_r, 640, 480, alphaG, betaG);
	float* q = new float[2 * 640 * 480];
	float* a = d_min.data(); // TODO: use minimized Eaux from paper
	float* d_ = new float[640 * 480];
	std::copy(a, a + 640 * 480, d_);
	float theta = theta_start;
	int n = 1;
	while (theta > theta_min) {
		updateQ(g, a, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);
		updateD(g, a, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);

		// TODO: update a - see equations 13 and 14 of the paper.

		float beta = (theta > 1e-3) ? 1e-3 : 1e-4;
		theta *= (1 - beta * n);
		n++;
	}
}
/* In this file we will implement the primal dual algorithms
   Please refer section 2.2.3 from dtam paper.
*/

#include <FreeImage.h>
#include "Eigen.h"
#include <cmath>
#include "parameters.h"
#include "utility.h"
#include "FreeImageHelper.h" 

void create_Eaux(Eigen::VectorXf& d_, float theta, float lambda, Eigen::MatrixXf& C_u_a, Eigen::MatrixXf& _Eaux);

// 6.2.3 of A. Chambolle and T. Pock. A first-order primal-dual
// algorithm for convex problems with applications to imaging.
float sigma_d(float epsilon, float theta) {

	float L = 4.0;

	float mu = 2.0*std::sqrt(epsilon / theta) / L;

	return mu / (2.0 / theta);
}
float sigma_q(float epsilon, float theta) {

	float L = 4.0;

	float mu = 2.0*std::sqrt(epsilon / theta) / L;

	return mu / (2.0*epsilon);
}
float toGray(float r, float g, float b) {
	return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}
void computeG(float* g, BYTE* img, int width, int height, float alphaG, float betaG) {
	// equation 5 of the paper
	tbb::parallel_for(0, pixels, [=](int i) {
		const int y = i / width;
		const int x = i % width;
		int idx = i * 4;
		const float gray = toGray(img[idx] / 255.0f, img[idx + 1] / 255.0f, img[idx + 2] / 255.0f);
		const float gray_right = (x == width - 1) ? 0.0f :
			toGray(img[idx + 4] / 255.0f,
				img[idx + 5] / 255.0f,
				img[idx + 6] / 255.0f);
		const float gray_left = (x == 0) ? 0.0f :
			toGray(img[idx - 4] / 255.0f,
				img[idx - 3] / 255.0f,
				img[idx - 2] / 255.0f);
		const float gray_down = (y == height - 1) ? 0.0f :
			toGray(img[idx + width * 4] / 255.0f,
				img[idx + width * 4 + 1] / 255.0f,
				img[idx + width * 4 + 2] / 255.0f);
		const float gray_up = (y == 0) ? 0.0f :
			toGray(img[idx - width * 4] / 255.0f,
				img[idx - width * 4 + 1] / 255.0f,
				img[idx - width * 4 + 2] / 255.0f);
		// forward difference gradient
		//float gx = gray_right - gray;
		//float gy = gray_down - gray;
		// central difference gradient
		float gx = gray_right - gray_left;
		float gy = gray_down - gray_up;
		g[i] = expf(-alphaG * powf(sqrtf(gx*gx + gy * gy), betaG));
	});
	//for (int y = 0; y < height; y++) {
	//	for (int x = 0; x < width; x++) {
	//		int idx = (y * width + x) * 4;
	//		// forward difference gradient
	//		const float gray = toGray(img[idx] / 255.0f, img[idx + 1] / 255.0f, img[idx + 2] / 255.0f);
	//		const float gray_right = (x == width - 1) ? 0.0f :
	//			toGray(img[idx + 4] / 255.0f,
	//				img[idx + 5] / 255.0f,
	//				img[idx + 6] / 255.0f);
	//		const float gray_down = (y == height - 1) ? 0.0f :
	//			toGray(img[idx + width * 4] / 255.0f,
	//				img[idx + width * 4 + 1] / 255.0f,
	//				img[idx + width * 4 + 2] / 255.0f);
	//		float gx = gray_right - gray;
	//		float gy = gray_down - gray;
	//		g[idx/4] = expf(-alphaG * powf(sqrtf(gx*gx + gy * gy), betaG));
	//	}
	//}
}
void updateQ(float* g, Eigen::VectorXf& a_, Eigen::VectorXf& q_, Eigen::VectorXf& d_, int width, int height, float sigma_q, float sigma_d, float epsilon, float theta) {
	float *a  = a_.data(); 
	float *q = q_.data();//new float[pixels * 2];
	float *d = d_.data();//new float[pixels];
	//Eigen::Map<Eigen::VectorXf>( a, a_.rows(), a_.cols() ) =  a_;
	//Eigen::Map<Eigen::VectorXf>( q, q_.rows(), q_.cols() ) =  q_;
	//Eigen::Map<Eigen::VectorXf>( d, d_.rows(), d_.cols() ) =  d_;

	tbb::parallel_for(0, pixels, [=](int i) {
		const int y = i / width;
		const int x = i % width;
		int idx = y * width + x;
		// gradient with forward difference
		// TODO use closed form
		float d_left = (x == 0) ? 0.0f : d[idx - 1];
		float d_right = (x == width - 1) ? 0.0f : d[idx + 1];
		//float dd_x = (x == width - 1) ? 0.0f : d[idx + 1] - d[idx];
		float dd_x = d_right - d_left;
		float d_up = (y == 0) ? 0.0f : d[idx - width];
		float d_down = (y == height - 1) ? 0.0f : d[idx + width];
		//float dd_y = (y == height - 1) ? 0.0f : d[idx + width] - d[idx];
		float dd_y = d_down - d_up;

		float qx = (q[idx] + sigma_q * g[idx] * dd_x) / (1.0f + sigma_q * epsilon);
		float qy = (q[idx + width * height] + sigma_q * g[idx] * dd_y) / (1.0f + sigma_q * epsilon);

		float maxq = std::fmaxf(1.0f, sqrtf(qx*qx + qy * qy));
		q[idx] = qx / maxq;
		q[idx + width * height] = qy / maxq;
	});
	//for (int y = 0; y < height; y++) {
	//	for (int x = 0; x < width; x++) {
	//		int idx = y * width + x;
	//		// gradient with forward difference
	//		// TODO use closed form
	//		float dd_x = (x == width - 1) ? 0.0f : d[idx + 1] - d[idx];
	//		float dd_y = (y == height - 1) ? 0.0f : d[idx + width] - d[idx];

	//		float qx = (q[idx] + sigma_q * g[idx] * dd_x) / (1.0f + sigma_q * epsilon);
	//		float qy = (q[idx + width * height] + sigma_q * g[idx] * dd_y) / (1.0f + sigma_q * epsilon);

	//		float maxq = std::fmaxf(1.0f, sqrtf(qx*qx + qy * qy));
	//		q[idx] = qx / maxq;
	//		q[idx + width * height] = qy / maxq;
	//	}
	//}
}
void updateD(float* g, Eigen::VectorXf& a_, Eigen::VectorXf& q_, Eigen::VectorXf& d_, int width, int height, float sigma_q, float sigma_d, float epsilon, float theta) {
	float *a = a_.data();
	float *q = q_.data();//new float[pixels * 2];
	float *d = d_.data();//new float[pixels];
	//Eigen::Map<Eigen::VectorXf>( a, a_.rows(), a_.cols() ) =  a_;
	//Eigen::Map<Eigen::VectorXf>( q, q_.rows(), q_.cols() ) =  q_;
	//Eigen::Map<Eigen::VectorXf>( d, d_.rows(), d_.cols() ) =  d_;

	tbb::parallel_for(0, pixels, [=](int i) {
		const int y = i / width;
		const int x = i % width;
		int idx = y * width + x;
		float qx_left = (x == 0) ? 0.0f : q[idx - 1];
		float qx_right = (x == width - 1) ? 0.0f : q[idx + 1];
		float qx_up = (y == 0) ? 0.0f : q[idx - width];
		float qx_down = (y == width - 1) ? 0.0f : q[idx + width];
		float qy_up = (y == 0) ? 0.0f : q[idx + width * height - width];
		float qy_down = (y == width - 1) ? 0.0f : q[idx + width * height + width];
		float qy_left = (x == 0) ? 0.0f : q[idx + width * height - 1];
		float qy_right = (x == width - 1) ? 0.0f : q[idx + width * height + 1];
		float dqx_x = (x == 0) ? q[idx] - q[idx + 1] : q[idx] - q[idx - 1];
		float dqy_y = (y == 0) ? q[idx + width * height] - q[idx + width * height + width] : q[idx + width * height] - q[idx + width * height - width];
		float div_q = dqx_x + dqy_y;
		//float dqx_x = qx_left - qx_right;
		//float dqx_y = qx_up - qx_down;
		//float dqy_y = qy_up - qy_down;
		//float dqy_x = qy_left - qy_right;
		//float div_q = sqrtf(dqx_x*dqx_x + dqx_y * dqx_y) + sqrtf(dqy_x*dqy_x + dqy_y * dqy_y);
		//div_q *= -1;
		d[idx] = (d[idx] + sigma_d * (g[idx] * div_q + a[idx] / theta)) / (1.0f + sigma_d / theta);
	});
	//for (int y = 0; y < height; y++) {
	//	for (int x = 0; x < width; x++) {
	//		int idx = y * width + x;
	//		// gradient with backward difference
	//		// TODO use closed form
	//		float dqx_x = (x == 0) ? q[idx] - q[idx + 1] : q[idx] - q[idx - 1];
	//		float dqy_y = (y == 0) ? q[idx + width * height] - q[idx + width * height + width] : q[idx + width * height] - q[idx + width * height - width];
	//		float div_q = dqx_x + dqy_y;

	//		d[idx] = (d[idx] + sigma_d * (g[idx] * div_q + a[idx] / theta)) / (1.0f + sigma_d / theta);
	//	}
	//}
}
void subsampleNewton(float* a, Eigen::MatrixXf& eaux, int width, int height) {
	// 2.2.5 from the paper
	tbb::parallel_for(0, pixels, [&](int idx) {
		int x = idx % width;
		int y = idx / width;
		int depth_index = (a[idx] - max_depth) / -inc_depth;
		if (depth_index == 0 || depth_index == d_range - 1) {
			return;
		}

		float a_energy = eaux(idx, depth_index + 1);
		float b_energy = a[idx];
		float c_energy = eaux(idx, depth_index - 1);
		float delta = ((a_energy + c_energy) == 2 * b_energy) ? 0.0f : ((a_energy - c_energy)*inc_depth) / (2 * (a_energy - 2 * b_energy + c_energy));
		delta = (fabsf(delta) > inc_depth) ? 0.0f : delta;
		a[idx] += delta;
	});
}

void saveDepthImage(const char* filename, Eigen::VectorXf& d_) {
	FreeImageB outImage(640, 480, 3);
	auto outData = new BYTE[640 * 480 * 3];
	tbb::parallel_for(0, pixels, [&](int idx) {
		// 255 = white, 0 = black
		// close is white and far is black
		outData[idx * 3] = 255 - (d_[idx] / max_depth) * 255;
		outData[idx * 3 + 1] = 255 - (d_[idx] / max_depth) * 255;
		outData[idx * 3 + 2] = 255 - (d_[idx] / max_depth) * 255;
	});
	outImage.data = outData;
	outImage.SaveImageToFile(filename);
}
void PrimalDual(int current_ref_img, BYTE* colorFrame_r, BYTE** colorFrames_b, Eigen::VectorXf& d_) {
	//_Eaux_Min = a_u for all pixels as written in the paper
	// d_ is d_u for all pixels per frame

	using namespace Eigen;
	const float off = d_range / 32.0f;
	float alphaG = 1.0f;//3.5f;//1.0f; // <- bad, i used 1
	float betaG = 0.1f;//1.0f;//0.1f; // <- bad, i used 0.1
	float theta_start = 0.2; // <- good
	float theta_min = 1.0e-4; // <- good
	float theta_step = 0.97; // <- good
	float epsilon = 0.01f * off;//0.00147; // <- good
	float lambda = 0.01f / off;//0.8f;//0.80; // <- good maybe 1.0
	float theta = theta_start;

	float* g = new float[640 * 480];
	computeG(g, colorFrame_r, 640, 480, alphaG, betaG);
	Eigen::VectorXf q(640*480*2);

	std::cout << __func__ << ":    " << "using Primal Dual Algo" << std::endl;

	Eigen::MatrixXf C_u_a(pixels, d_range);
	for(uint i = 0; i < pixels; ++i){
        for(uint j=0; j< d_range; j++)
            C_u_a(i,j) = 0;
    }

	cost_calc(colorFrame_r, colorFrames_b, current_ref_img, C_u_a);

	//initializing "_Eaux_Min" vector by doing argmin(cost ie C_a)
	Eigen::VectorXf _Eaux_Min(pixels);
	tbb::parallel_for(0, pixels, [&](int idx) {
		int index = 0;
		C_u_a.row(idx).minCoeff(&index);
		_Eaux_Min(idx) = max_depth - index * inc_depth;
	});
	tbb::parallel_for(0, pixels, [&](int idx) {
		int index = 0;
		C_u_a.row(idx).minCoeff(&index);
		d_(idx) = max_depth - index * inc_depth;
	});

	//first time initialization of d_u = a_u. check if this assignment makes a deep copy.
	//d_ = _Eaux_Min;
/*
	updateQ(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);
	updateD(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);*/

	saveDepthImage("start.png", d_);
	int n = 1;
	Eigen::MatrixXf _Eaux(pixels, d_range);
	while (theta > theta_min) {

		updateQ(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);
		updateD(g, _Eaux_Min, q, d_, 640, 480, sigma_q(epsilon, theta), sigma_d(epsilon, theta), epsilon, theta);

		if (n % 10 == 9) {
			std::stringstream ss;
			ss << "out_" << n << ".png";
			saveDepthImage(ss.str().c_str(), d_);
		}
		create_Eaux(d_, theta, lambda, C_u_a, _Eaux);

		tbb::parallel_for(0, pixels, [&](int idx) {
			int index = 0;
			_Eaux.row(idx).minCoeff(&index);
			_Eaux_Min(idx) = max_depth - index * inc_depth;
		});

		subsampleNewton(_Eaux_Min.data(), _Eaux, 640, 480);

		float beta = (theta > 1e-3) ? 1e-3 : 1e-4;
		theta *= (1 - beta * n);
		n++;
	}

}

void create_Eaux(Eigen::VectorXf& d_, float theta, float lambda, Eigen::MatrixXf& C_u_a, Eigen::MatrixXf& _Eaux){
	tbb::parallel_for(0, pixels, [&](int pixel){
		for (uint i = 0; i < d_range; i++) {
			float a_i = max_depth - i * inc_depth;
			float E = 1.0f/(2.0f*theta) * (d_[pixel] - a_i) * (d_[pixel] - a_i) + lambda * C_u_a(pixel, i);
			_Eaux(pixel, i) = E;
		}
	});
}
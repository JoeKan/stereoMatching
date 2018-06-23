void load_all_matrices_from_n_files(std::vector<Eigen::Matrix4f> &P) {
const int n_images = N_FRAMES*100;
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

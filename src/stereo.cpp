#include "Eigen.h"
#include "parameters.h"
#include "utility.h"
#include "FreeImageHelper.h" 

bool processNextFrame(int filenum, BYTE* colorFrame);
void brute_force_depth_calc(BYTE* colorFrame_r);
std::vector<Eigen::Matrix4f> poses;

Eigen::Matrix4f K;
Eigen::Matrix4f KInv;

//need an error heat map (red is very off, blue is good, yellow is middle ground)

#define uint unsigned int

int main(int argc, char* argv[]){

    char mode;

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "\t--mode,-m To use Primal dual(p) optimization or brute force(b), eg -m p or -m b" << std::endl;
        return 1;
    }
    std::vector <std::string> sources;
    std::string destination;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--mode" || std::string(argv[i]) == "-m") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                destination = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--mode option requires one argument." << std::endl;
                return 1;
            }  
        } else {
            std::cerr << "Usage: " << argv[0] << "\t--mode,-m To use Primal dual(p) optimization or brute force(b), eg -m p or -m b" << std::endl;
            return 1;
        }
    }

    mode = destination[0];

    BYTE* colorFrame_r = new BYTE[4* 640*480];
    read_images(current_ref_img, colorFrame_r);
    
    BYTE** colorFrames_b = new BYTE*[m];
    for(uint i=0;i<m;++i){
        colorFrames_b[i] = new BYTE[4* 640*480];
    }

    for(uint i=0;i<m;++i){
        read_images(current_ref_img - m + i, colorFrames_b[i]);
    }
    
    K << 577.97f, 0.0f, 320.0f, 0,
		 0.0f, 575.44f, 240.0f, 0,
		 0.0f, 0.0f, 1.0f, 0,
         0, 0, 0, 1;

    KInv = K.inverse();
    load_all_matrices_from_n_files(poses);
    
    //------Brute Force-------------
    if( mode == 'b'){
        processNextFrame(current_ref_img, colorFrame_r);
        brute_force_depth_calc(colorFrame_r);
    }

    //------Primal_Dual-------------
    if( mode == 'p'){
        Eigen::VectorXf d_(pixels);
        PrimalDual(current_ref_img, colorFrame_r, colorFrames_b, d_);
        FreeImageB outImage(640, 480, 3);
        BYTE* outData = new BYTE[640 * 480 * 3];
        tbb::parallel_for(0, pixels, [&](int idx) {
            // 255 = white, 0 = black
            // close is white and far is black
            outData[idx * 3] = 255 - (d_[idx] / max_depth) * 255;
            outData[idx * 3 + 1] = 255 - (d_[idx] / max_depth) * 255;
            outData[idx * 3 + 2] = 255 - (d_[idx] / max_depth) * 255;
        });
        outImage.data = outData;
        outImage.SaveImageToFile("out_primaldual.png");
    }

    if(mode != 'p' && mode != 'b'){
        std::cerr << "Usage: " << argv[0] << "\t--mode,-m To use Primal dual(p) optimization or brute force(b), eg -m p or -m b" << std::endl;
        std::cerr << "\t Use correct option b or p" << std::endl;
        return 1;
    }
    return 0;
}
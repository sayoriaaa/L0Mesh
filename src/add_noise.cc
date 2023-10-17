#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/edge_lengths.h>
#include <iostream>

Eigen::MatrixXd uniform(const Eigen::MatrixXd& V, double sigma) {
    // X ~ U(-sigma, sigma)
    Eigen::MatrixXd noise = sigma * Eigen::MatrixXd::Random(V.rows(), V.cols());
    return noise;
}

Eigen::MatrixXd gaussian(const Eigen::MatrixXd& V, double sigma) {
    // X ~ U(-sigma, sigma)
    Eigen::MatrixXd noise = Eigen::MatrixXd::Zero(V.rows(), V.cols());
    for(int i=0; i<4; i++){
        noise += uniform(V, sigma);
    }
    noise *= 0.25;
    return noise;
}

int main(int argc, char *argv[])
{   
    double factor = 0.3;
    bool showHelp = false;
    int type = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-f" || arg == "--factor") {
            if (i + 1 < argc) {
                factor = std::stod(argv[i + 1]);
                i++; // Skip the next argument (number)
            } else {
                std::cerr << "-f, --factor requires an argument." << std::endl;
                return 1;
            }
        } else if (arg == "-u" || arg == "--uniform") {
            type = 1;
        } else if (arg == "-g" || arg == "--gaussian") {
            type = 0;
        } else if (arg == "-h" || arg == "--help") {
            showHelp = true;
        }
    }

    if (showHelp) {
        std::cout << "Usage: " << argv[0] << " [options]\n"
                  << "Options:\n"
                  << "  -f, --factor <number> set sigma as f*average edge length\n"
                  << "  -g, --gaussian use gaussian noise\n"
                  << "  -u, --uniform use uniform noise\n"
                  ;
        return 0;
    }

    // end parsing

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd L;

    // read mesh 
    igl::readOBJ( "../examples/cube.obj", V, F);
    igl::edge_lengths(V, F, L);

    // calc average edge length
    double average_length = L.mean();
    std::cout << "Average edge length: " << average_length << std::endl;

    // build noise matrix (avoid the usage of for)
    Eigen::MatrixXd noise;
    double sigma = average_length * factor;
    if(type==0) noise = uniform(V, sigma);
    else if(type==1) noise = gaussian(V, sigma);
    // add noise
    V += noise;
    igl::writeOBJ("../examples/cube_.obj", V, F);

}

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/cotmatrix.h>
#include <Eigen/Sparse>
#include <igl/edges.h>
#include <igl/edge_lengths.h>
#include <Eigen/Cholesky>
#include <iostream>
#include <string>

double average_dihedral(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    Eigen::MatrixXi E; 
    igl::edges(F, E);

    double sum = 0;

    for (int i = 0; i < E.rows(); i++)
    {
        int vertex1 = E(i, 0); 
        int vertex2 = E(i, 1); 

        int vertex3 = -1, vertex4 = -1;

        // find vertex3 vertex4
        for (int j = 0; j < F.rows(); j++)
        {
            int v0 = F(j, 0);
            int v1 = F(j, 1);
            int v2 = F(j, 2);

            if ((v0 == vertex1 && v1 == vertex2) || (v0 == vertex2 && v1 == vertex1) ||
                (v1 == vertex1 && v2 == vertex2) || (v1 == vertex2 && v2 == vertex1) ||
                (v2 == vertex1 && v0 == vertex2) || (v2 == vertex2 && v0 == vertex1))
            {
                int find_v;
                if ((v0 == vertex1 && v1 == vertex2) || (v0 == vertex2 && v1 == vertex1)) find_v = v2;
                if ((v1 == vertex1 && v2 == vertex2) || (v1 == vertex2 && v2 == vertex1)) find_v = v0;
                if ((v2 == vertex1 && v0 == vertex2) || (v2 == vertex2 && v0 == vertex1)) find_v = v1;
                
                if (vertex3 == -1) vertex3 = find_v;
                else {
                    vertex4 = find_v;
                    break;
                }
            }
        }

        // calc weights
        //       v1
        //     /  |  \
        //    /   |   \
        //   v3   |    v4
        //   \    |    /
        //    \   |   /
        //        v2
        Eigen::Vector3d p1 = V.row(vertex1);
        Eigen::Vector3d p2 = V.row(vertex2);
        Eigen::Vector3d p3 = V.row(vertex3);
        Eigen::Vector3d p4 = V.row(vertex4);

        Eigen::Vector3d norm1 = ((p1-p2).cross(p3-p2)).normalized();
        Eigen::Vector3d norm2 = ((p4-p2).cross(p1-p2)).normalized();

        double cos_theta = norm1.dot(norm2);
        double angle_rad = std::acos(cos_theta);
        sum += angle_rad;
    }
    return sum/E.rows();
}

void cotEdge(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& L) {
    Eigen::MatrixXi E; 
    igl::edges(F, E);

    // init Laplacian operator
    L.resize(E.rows(),V.rows());
    L.reserve(10*E.rows()); // 3-simplex, following setting of libigl

    std::vector<Eigen::Triplet<double> > tripletList;

    for (int i = 0; i < E.rows(); i++)
    {
        int vertex1 = E(i, 0); 
        int vertex2 = E(i, 1); 

        int vertex3 = -1, vertex4 = -1;

        // find vertex3 vertex4
        for (int j = 0; j < F.rows(); j++)
        {
            int v0 = F(j, 0);
            int v1 = F(j, 1);
            int v2 = F(j, 2);

            if ((v0 == vertex1 && v1 == vertex2) || (v0 == vertex2 && v1 == vertex1) ||
                (v1 == vertex1 && v2 == vertex2) || (v1 == vertex2 && v2 == vertex1) ||
                (v2 == vertex1 && v0 == vertex2) || (v2 == vertex2 && v0 == vertex1))
            {
                int find_v;
                if ((v0 == vertex1 && v1 == vertex2) || (v0 == vertex2 && v1 == vertex1)) find_v = v2;
                if ((v1 == vertex1 && v2 == vertex2) || (v1 == vertex2 && v2 == vertex1)) find_v = v0;
                if ((v2 == vertex1 && v0 == vertex2) || (v2 == vertex2 && v0 == vertex1)) find_v = v1;
                
                if (vertex3 == -1) vertex3 = find_v;
                else {
                    vertex4 = find_v;
                    break;
                }
            }
        }

        // calc weights
        //       v1
        //     /  |  \
        //    /   |   \
        //   v3   |    v4
        //   \    |    /
        //    \   |   /
        //        v2
        Eigen::Vector3d p1 = V.row(vertex1);
        Eigen::Vector3d p2 = V.row(vertex2);
        Eigen::Vector3d p3 = V.row(vertex3);
        Eigen::Vector3d p4 = V.row(vertex4);

        double cot312 = ((p3-p1).dot(p2-p1)) / ((p3-p1).cross(p2-p1)).norm();
        double cot412 = ((p4-p1).dot(p2-p1)) / ((p4-p1).cross(p2-p1)).norm();
        double cot321 = ((p3-p2).dot(p1-p2)) / ((p3-p2).cross(p1-p2)).norm();
        double cot421 = ((p4-p2).dot(p1-p2)) / ((p4-p2).cross(p1-p2)).norm();

        double coef1 = -1 * (cot321 + cot421);
        double coef2 = -1 * (cot312 + cot412);
        double coef3 = cot321 + cot312;
        double coef4 = cot421 + cot412;

        tripletList.push_back(Eigen::Triplet<double> (i, vertex1, coef1));
        tripletList.push_back(Eigen::Triplet<double> (i, vertex2, coef2));
        tripletList.push_back(Eigen::Triplet<double> (i, vertex3, coef3));
        tripletList.push_back(Eigen::Triplet<double> (i, vertex4, coef4));

    }

    L.setFromTriplets(tripletList.begin(), tripletList.end());
}

int main(int argc, char *argv[])
{   
    bool showHelp = false;
    std::string inputFile = "../examples/cube_.obj";
    std::string outputFile = "../examples/denoised.obj";

    double lambda = -1;
    double kappa = 1.414;
    double beta_max = 1000;
    int type = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) {
                inputFile = argv[i + 1];
                i++; // Skip the next argument (file name)
            } else {
                std::cerr << "-i, --input requires an argument." << std::endl;
                return 1;
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                outputFile = argv[i + 1];
                i++; // Skip the next argument (file name)
            } else {
                std::cerr << "-o, --output requires an argument." << std::endl;
                return 1;
            }
        } else if (arg == "-l" || arg == "--lambda") {
            if (i + 1 < argc) {
                lambda = std::stod(argv[i + 1]);
                i++; // Skip the next argument (number)
            } else {
                std::cerr << "-l, --lambda requires an argument." << std::endl;
                return 1;
            }
        } else if (arg == "-k" || arg == "--kappa") {
            if (i + 1 < argc) {
                kappa = std::stod(argv[i + 1]);
                i++; // Skip the next argument (number)
            } else {
                std::cerr << "-k, --kappa requires an argument." << std::endl;
                return 1;
            }
        } else if (arg == "-bm" || arg == "--beta_max") {
            if (i + 1 < argc) {
                beta_max = std::stod(argv[i + 1]);
                i++; // Skip the next argument (number)
            } else {
                std::cerr << "-bm, --beta_max requires an argument." << std::endl;
                return 1;
            }
        } else if (arg == "-e" || arg == "--edge") {
            type = 1;
        } else if (arg == "-v" || arg == "--vertex") {
            type = 0;
        } else if (arg == "-h" || arg == "--help") {
            showHelp = true;
        }
    }

    if (showHelp) {
        std::cout << "Usage: " << argv[0] << " [options]\n"
                  << "Options:\n"
                  << "  -i, --input <file>  Input file\n"
                  << "  -o, --output <file>  Output file\n"
                  << "  -l, --lambda <number> control balance between L0 and similarity\n"
                  << "  -k, --kappa <number> control convergence speed\n"
                  << "  -bm, --beta_max <number> control convergence up thres\n"
                  << "  -v, --vertex use vertex based Laplacian\n"
                  << "  -e, --edge use edge based Laplacian\n"
                  ;
        return 0;
    }
    // end parsing

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> D; // Laplacian operator

    // read mesh 
    igl::readOBJ(inputFile, V, F);

    Eigen::SparseMatrix<double> I(V.rows(), V.rows());
    I.setIdentity();

    if(lambda==-1){
        std::cout << "set lambda as default: 0.02*l^2_e*gamma\n";
        double gamma = average_dihedral(V, F);
        Eigen::MatrixXd L;
        igl::edge_lengths(V, F, L);
        double average_length = L.mean();
        lambda = 0.02 * average_length * average_length * gamma;
        std::cout << "Average dihedral angle: " << (gamma * 180.0 / 3.1415) << "\n"
                  << "Average edge length: " << average_length << "\n"
                  << "auto lambda: " << lambda << "\n"
                  ;
    }

    // start opt
    double beta = 1.0e-3;
    Eigen::MatrixXd p = V;
    while(beta < beta_max){
        // build Laplacian operator
        std::cout << "beta: " << beta << std::endl;
        if(type==0) igl::cotmatrix(p, F, D);
        else if(type==1) cotEdge(p, F, D);
        // local optimization
        Eigen::MatrixXd delta = D * p;
        for (int i = 0; i < delta.rows(); ++i) {
            if (delta.row(i).squaredNorm() < lambda / beta) {
                delta.row(i).setZero();
            }
        }
        // global optimization
        Eigen::SparseMatrix<double> A = I + beta * D.transpose() * D;
        Eigen::MatrixXd b = V + beta * (D.transpose() * delta);

        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);

        if(solver.info() != Eigen::Success) {
            std::cerr << "Cholesky decomposition failed" << std::endl;
            return -1;
        }


        for (int i = 0; i < b.cols(); ++i) {
            Eigen::VectorXd col_b = b.col(i);
            Eigen::VectorXd col_x = solver.solve(col_b);
            p.col(i) = col_x;
        }

        if(solver.info() != Eigen::Success) {
            std::cerr << "Cholesky solve failed" << std::endl;
            return -1;
        }

        // update parameter
        beta *= kappa;

    }

    // write mesh
    igl::writeOBJ(outputFile, p, F);

}

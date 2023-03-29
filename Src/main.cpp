#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "IO.hpp"
#include "micromechanics.hpp"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char *argv[]) {
    if (argc != 8) {
        std::cerr << "command line arguments: <input_material> <initial_strain> <output_file> <e/s> <xx/yy/zz/xy/yx/xz/zx/yz/zy> <tolerance> <maximum number of iterations>" << std::endl;
        return 1;
    }

    std::string in_mat_filename = argv[1];
    std::string in_strain_filename = argv[2];
    std::string out_filename = argv[3];
    std::string e_or_s = argv[4];
    std::string index = argv[5];
    double tol = std::stod(argv[6]);
    int maxit = std::stoi(argv[7]);

    std::ifstream in_file(in_mat_filename); 
    if (!in_file) {
        std::cerr << "could not find input_material file" << std::endl;
        return 1;
    }
    in_file.close();

    in_file.open(in_strain_filename);
    if (!in_file) {
        std::cerr << "could not find initial_strain file" << std::endl;
        return 1;
    }
    in_file.close();

    if(e_or_s != "e" && e_or_s != "s") {
        std::cerr << "<e/s> should be e or s" << std::endl;
        return 1;
    } 

    if(index != "xx" && index != "yy" && index != "zz" && index != "xy" && index != "yx" && index != "xz" && index != "zx" && index != "yz" && index != "zy") {
        std::cerr << "the last argument should be xx/yy/zz/xy/yx/xz/zx/yz/zy" << std::endl;
        return 1;
    } 


    auto my_args = read_material(in_mat_filename);
    Eigen::Array<double, 6, 1> E = read_strain(in_strain_filename);


    Eigen::ArrayXXd c = my_args.c;
    Eigen::Tensor<int, 3> mat = my_args.mat;
    Eigen::Array<double, 3, 1> prds = my_args.prds;

    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);
    Eigen::Tensor<double, 4> t;

    m.iteration();

    if (e_or_s == "e") {
        t = m.getStrain();
    } else if (e_or_s == "s") {
        t = m.getStress();
    } 

    int idx;

    if (index == "xx") {
        idx = 0;
    } else if (index == "yy") {
        idx = 1;
    } else if (index == "zz") {
        idx = 2;
    } else if (index == "yz" || index == "zy") {
        idx = 3;
    } else if (index == "xz" || index == "zx") {
        idx = 4;
    } else if (index == "xy" || index == "yx") {
        idx = 5;
    }
  
    write_material(prds, t, idx, out_filename);

    return 0;
}
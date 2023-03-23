#ifndef read_material_hpp
#define read_material_hpp

#include <string>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"

// read in the Lamé coefficients of every material present in the cuboid and read in the material for every voxel in the cuboid 
auto read_material(std::string const & filename) {
    std::string line, lambda, mu, t;
    int num_mat, xdim, ydim, zdim;
    
    std::ifstream file(filename);

    // read in number of materials present in the cuboid and the Lamé coefficients for every material
    std::getline(file, line);
    std::istringstream split_line(line);
    split_line >> num_mat;

    Eigen::ArrayXXd c(2, num_mat);

    for (int i = 0; i < num_mat; i++) {
        std::getline(file, line);
        split_line = std::istringstream(line);
        split_line >> lambda;
        split_line >> mu;
        c(0, i) = std::stod(lambda);
        c(1, i) = std::stod(mu);
    }

    // read in the dimensions of the cuboid and the material index for every voxel
    std::getline(file, line);
    split_line = std::istringstream(line);
    split_line >> xdim;
    split_line >> ydim;
    split_line >> zdim;

    Eigen::Tensor<int, 3> mat(xdim, ydim, zdim);

    for (int k = 0; k < zdim; k++) {
        for (int j = 0; j < ydim; j++) {
            for (int i = 0; i < xdim; i++) {
                std::getline(file, line);
                mat(i, j, k) = std::stoi(line);
            }
        }
    }

    // read in the period for every dimension of the cuboid
    Eigen::ArrayXd prds(3);

    std::getline(file, line);
    split_line = std::istringstream(line);
    for (int i = 0; i < 3; i++) {
        split_line >> t;
        prds(i) = std::stod(t); 
    }

    file.close();

    // return all arguments in a custom struct
    struct args {
        Eigen::ArrayXXd c;
        Eigen::Tensor<int, 3> mat;
        Eigen::ArrayXd prds;
    };

    return args{c, mat, prds};
}

#endif

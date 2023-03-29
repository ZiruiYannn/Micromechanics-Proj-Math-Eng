#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

typedef double precision_type;
typedef int size_type;
typedef Eigen::Array<size_type, 3, 1> Array3;
typedef Eigen::Array<precision_type, 3 ,1> Vec3;
typedef Eigen::Array<precision_type, 6 ,1> Vec6;



#endif
#ifndef MICROMECHANICS_HPP_
#define MICROMECHANICS_HPP_

#include "types.hpp"

namespace mme {

    template <typename type_mat, typename type_c, typename maxit_type>
    class micromechanics {
        private:
            Array3 dims_;
            
            Vec6 strain_0;
            type_mat mat_;
            type_c c_;
            precision_type lamda_ref;
            precision_type mu_ref;
            Vec3 prds_;
            
            precision_type tol_;
            maxit_type maxit_;
            

        public:            
            micromechanics(Vec6 E, type_mat mat, type_c c, Vec3 prds, precision_type tol, maxit_type maxit): 
            strain_0(E), mat_(mat), c_(c), prds_(prds), tol_(tol), maxit_(maxit) {
                // mat(matrial number, x, y, z)
                auto dims_mat = mat.dimensions();
                for (size_type i=0; i<3; i++) {
                    dims_[i] = dims_mat[i+1];
                }
                // c is a matrix, the first row is lamda, the second row is mu
                auto dims_c = c.dimensions();
                lamda_ref = 0;
                mu_ref = 0;
                for (size_type i=0; i<dims_c[0];i++) {
                    lamda_ref += c(0,i);
                }
                lamda_ref /= dims_c[0];
                for (size_type i=0; i<dims_c[1];i++) {
                    mu_ref += c(1,i);
                }
                mu_ref /= dims_c[1];

            }

            ~micromechanics() {}
        
        public:
            typedef typename Eigen::Tensor<precision_type, 6, dims_[0], dims_[1], dims_[2]> TVec6;
            typedef typename Eigen::Tensor<precision_type, 3, dims_[0], dims_[1], dims_[2]> TVec3;
            typedef typename Eigen::Tensor<precision_type, 3, 3, 3, 3> TGreen;

        private:
            Vec6 stressCompute(Vec6 eps, precision_type lamda, precision_type mu) {
                precision_type sum =  eps[0] + eps[1] + eps[2];
                Vec6 sig;
                sig = eps;

                return sig;
            }

            precision_type error(TVec6 sig) {
                return 
            }

            Vec3 waveVec(Array3 inds) {

            }

            TGreen greenOp(Vec3 xi) {

            }

            Vec6 polarization(Vec6 sig, Vec6 eps) {

            }
        
        private:
            TVec6 stress_;
            TVec6 strain_;
            
        public:
            //do the iteration and compute stress_ and strain_
            void iteration() {

            }

            TVec6 getStress() const {
                return stress_;
            }

            TVec6 getStress() const {
                return strain_;
            }

    };

};

#endif
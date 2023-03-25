#ifndef MICROMECHANICS_HPP_
#define MICROMECHANICS_HPP_


#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"


namespace mme {

    template <typename Precision>
    class micromechanics {
        public:
            Eigen::Array<int,3> dims_;
            
            Eigen::Array<Precision, 6> strain_0;
            Eigen::Tensor<Precision, 4> mat_;
            Eigen::Tensor<Precision, 2> c_;
            Precision lamda_ref;
            Precision mu_ref;
            Eigen::Array<Precision, 3> prds_;
            
            precision tol_;
            int maxit_;

            Eigen::Tensor<Precision, 4> strain_;
            Eigen::Tensor<Precision, 4> stress_;
            

        public:            
            micromechanics(Eigen::Array<Precision, 6> E, Eigen::Tensor<Precison,4> mat, Eigen::Tensor<Precision, 2> c, \
            Eigen::Array<Precision, 3> prds, Precision tol, int maxit): strain_0(E), mat_(mat), c_(c), prds_(prds), tol_(tol), \
            maxit_(maxit) {

                //initial dims_
                // mat(matrial number, x, y, z)
                auto dims_mat = mat.dimensions();
                for (int i=0; i<3; i++) {
                    dims_(i) = dims_mat(i+1);
                }

                //initial lamda_ref and mu_ref
                // c is a matrix, the first row is lamda, the second row is mu
                Eigen::array<int, 1> reduce_dim={0}; 
                Eigen::Tensor<int, 1> max_c = c.maximum(reduce_dim); //store the max value of each row
                Eigen::Tensor<int, 1> min_c = c.minimum(reduce_dim); //store the min value of each row

                lamda_ref = (max_c(0) + min_c(0))/2;
                mu_ref = (max_c(1) + min_c(1))/2;

                //initial strain_ and stress_
                Eigen::Tensor<Precision,4> eps(6,dims_mat(0),dims_mat(1),dims_mat(2));
                Eigen::Tensor<Precision,4> sig(6,dims_mat(0),dims_mat(1),dims_mat(2));

                for (int k = 0; k < dims_(2); k++) {
                    for (int j = 0; j < dims_(1); j++) {
                        for (int i = 0; i < dims_(0); i++) {
                            for (int l=0; l < 6; l++) {
                                eps(l,i,j,k) = E(l);
                            }
                        }
                    }
                }

                strain_ = eps;

                for (int k = 0; k < dims_(2); k++) {
                    for (int j = 0; j < dims_(1); j++) {
                        for (int i = 0; i < dims_(0); i++) {
                            for (int l=0; l < 6; l++) {
                                sig(l,i,j,k) = 0;
                            }
                        }
                    }
                }

                stress_ = sig;

            }

            ~micromechanics() {}

        private:
            Eigen::Array<Precision, 6> stressCompute(Eigen::Array<Precision, 6> epsVec, Precision lamda, Precision mu) {
                Precision sum =  epsVec[0] + epsVec[1] + epsVec[2];
                Vec6 sig;
                sig = eps;

                return sig;
            }

            Precision error(Eigen::Tensor<Precision,4> sig) {
                return 
            }

            Eigen::Array<Precision,3> waveVec(Eigen::Array<int,3> inds) {

            }

            Eigen::Tensor<Precision, 4> greenOp(Eigen::Array<Precision,3> xi) {

            }

            Eigen::Array<Precision, 6> polarization(Eigen::Array<Precision, 6> sigVec, Eigen::Array<Precision, 6> epsVec) {

            }
            
        public:
            //do the iteration and compute stress_ and strain_
            void iteration() {

            }

            Eigen::Tensor<Precision, 4> getStress() const {
                return stress_;
            }

            Eigen::Tensor<Precision, 4> getStress() const {
                return strain_;
            }

    };

};

#endif
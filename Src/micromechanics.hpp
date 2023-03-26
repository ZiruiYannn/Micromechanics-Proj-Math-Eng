#ifndef MICROMECHANICS_HPP_
#define MICROMECHANICS_HPP_


#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include <complex>
#include <cmath>
#include <fftw3.h>


namespace mme {

    template <typename Precision>
    class micromechanics {
        public:
            Eigen::Array<int,3, 1> dims_;
            
            Eigen::Array<Precision, 6, 1> strain_0;
            Eigen::Tensor<Precision, 3> mat_;
            Eigen::Tensor<Precision, 2> c_;
            Precision lamda_ref;
            Precision mu_ref;
            Eigen::Array<Precision, 3, 1> prds_;
            
            precision tol_;
            int maxit_;

            Eigen::Tensor<Precision, 4> strain_;
            Eigen::Tensor<Precision, 4> stress_;
            
            bool reach_maxit = false;
        public:            
            micromechanics(Eigen::Array<Precision, 6, 1> E, Eigen::Tensor<Precison,3> mat, Eigen::Tensor<Precision, 2> c, \
            Eigen::Array<Precision, 3, 1> prds, Precision tol, int maxit): strain_0(E), mat_(mat), c_(c), prds_(prds), tol_(tol), \
            maxit_(maxit) {

                //initial dims_
                auto dims_mat = mat.dimensions();
                for (int i=0; i<3; i++) {
                    dims_(i) = dims_mat(i);
                }

                //initial lamda_ref and mu_ref
                // c is a matrix, the first row is lamda, the second row is mu
                Eigen::array<int, 1, 1> reduce_dim={0}; 
                Eigen::Tensor<int, 1> max_c = c.maximum(reduce_dim); //store the max value of each row
                Eigen::Tensor<int, 1> min_c = c.minimum(reduce_dim); //store the min value of each row

                lamda_ref = (max_c(0) + min_c(0))/2;
                mu_ref = (max_c(1) + min_c(1))/2;

                //initial strain_ and stress_
                Eigen::Tensor<Precision,4> eps(6,dims_mat(0),dims_mat(1),dims_mat(2));
                Eigen::Tensor<Precision,4> sig(6,dims_mat(0),dims_mat(1),dims_mat(2));
                
                Eigen::Array<Precision, 6, 1> sigVec;

                for (int k = 0; k < dims_(2); k++) {
                    for (int j = 0; j < dims_(1); j++) {
                        for (int i = 0; i < dims_(0); i++) {
                            for (int l=0; l < 6; l++) {
                                eps(l,i,j,k) = E(l,0);
                            }
                        }
                    }
                }

                strain_ = eps;

                for (int k = 0; k < dims_(2); k++) {
                    for (int j = 0; j < dims_(1); j++) {
                        for (int i = 0; i < dims_(0); i++) {
                            for (int l=0; l < 6; l++) {
                                sig(l,i,j,k) = stressCompute(E, c(0, mat(i , j, k)), c(1, mat(i , j, k)));
                            }
                        }
                    }
                }

                stress_ = sig;

            }

            ~micromechanics() {}

        public:
            Eigen::Array<Precision, 6, 1> stressCompute(Eigen::Array<Precision, 6, 1> epsVec, Precision lamda, Precision mu) {
                Precision sum =  epsVec[0] + epsVec[1] + epsVec[2];
                Eigen::Array<Precision, 6, 1> sigVec;
                sigVec = 2*mu*epsVec;
                for (int i = 0; i<3; i++) {
                    sig(i,0) += lamda*sum;
                }

                return sig;
            }

            Precision error(Eigen::Tensor<Precision,4> sig) {
                Precision e;

                return e;
            }

            Eigen::Array<Precision,3, 1> waveVec(Eigen::Array<int, 3, 1> inds) {
                
            }

            Eigen::Tensor<Precision, 4> greenOp(Eigen::Array<Precision, 3, 1> xi) {
                Eigen::Tensor<Precision, 4> gam(3,3,3,3);
                gam.setZero();
                Precision coef_a, coef_b, xi_square;
                Eigen::Tensor<Precision, 4> temp;
                xi_square = xi(0,0) * xi(0,0) + xi(1,0) * xi(1,0) + xi(2,0) * xi(2,0);
                coef_a = 4 * mu_ref * xi_square;
                coef_b = lamda_ref + mu_ref;
                coef_b /= (mu_ref * (lamda_ref + 2*mu_ref)) * xi_square * xi_square;

                for (j=0 ; j<3; j++) {
                    for (i=0; i<3; i++) {
                        for (h=0; h<3; h++) {
                            for (k=0; k<3; k++) {
                                if (k == i) {
                                    gam(k, h ,i, j) += xi(h, 0) * xi(j, 0);                            
                                }
                                if (h == i) {
                                    gam(k, h ,i, j) += xi(k, 0) * xi(j, 0);                            
                                }
                                if (k == j) {
                                    gam(k, h ,i, j) += xi(h, 0) * xi(i, 0);                            
                                }
                                if (h == j) {
                                    gam(k, h ,i, j) += xi(k, 0) * xi(i, 0);                            
                                }
                                gam(k ,h, i ,j) /= a; 
                                gam(k, h ,i, j) -= b * xi(h, 0) * xi(j, 0) * xi(k, 0) * xi(i, 0);
                            }
                        }
                    }
                }
                return gam;
            }

            Eigen::Array<Precision, 6, 1> polarization(Eigen::Array<Precision, 6, 1> sigVec, Eigen::Array<Precision, 6, 1> epsVec) {
                Eigen::Array<Precision, 6, 1> tauVec, c_eps;
                c_eps = stressCompute(epsVec, lamda_ref, mu_ref);
                tauVec = sigVec - c_eps;
                return tauVec;
            }

            Eigen::Array<Precision, 6, 1> ten2vec(Eigen::Tensor<Precision, 2> T) {
                Eigen::Array<Precision, 6, 1> V;
                V(0,0) = T(0,0);
                V(1,0) = T(1,1);
                V(2,2) = T(2,2);
                V(3,0) = T(0,1); //xy
                V(4,0) = T(0,2); //xz
                V(5,0) = T(1,2); //yz
                
                return V; 
            }

            Eigen::Tensor<Precision, 2> vec2ten(Eigen::Array<Precision, 6, 1> V) {
                Eigen::Tensor<Precision, 2> T(3,3);
                for (int i=0; i<3; i++) {
                    T(i,i) = V(i,0);
                }

                T(0,1) = V(3,0);
                T(1,0) = V(3,0);
                T(0,2) = V(4,0);
                T(2,0) = V(4,0);
                T(1,2) = V(5,0);
                T(2,1) = V(5,0);

                return T;
            }


            //real2freq
            Eigen::Tensor<std::complex<Precision>, 4> r2f(Eigen::Tensor<Precision, 4> data_r) {
                typedef std::complex<Precision> fft_complex_p;

                int N_in = dims_(0,0) * dims_(1,0) * dims_(2,0);
                int N_out = dims_(0,0) * dims_(1,0) * (dims_(2,0)/2+1);

                Precision* in = (Precision*)fftw_malloc(sizeof(Precision) * N_in);
                fftw_complex_p* out = (fftw_complex_p*)fftw_malloc(sizeof(fft_complex_p) * N_out);               
                
                Eigen::Tensor<std::complex<std::complex<Precision>, 4> data_f(6,dims_(0,0),dims_(1,0),dims_(2,0)/2+1);

                fftw_plan plan = fftw_plan_dft_r2c_3d(dims_(0,0), dims_(1,0), dims_(2,0),\
                 in, out, FFTW_ESTIMATE);

                //eigen col-major, fftw row-major => j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l
                for (int i=0; i<6; i++) {
                    for (int l=0; l<dims_(2,0); l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; i<dims_(0,0); j++) {
                                in[j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l] = data_r(i,j,k,l);
                            }
                        }
                    }

                    fftw_execute(plan);

                    for (int l=0; l<(dims_(2,0)/2 + 1); l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; i<dims_(0,0); j++) {
                                data_f(i,j,k,l) = out[j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l];
                            }
                        }
                    }
                }

                fftw_destroy_plan(plan);
                fftw_free(in);
                fftw_free(out);

                return data_f;
            } 

            //freq2real
            Eigen::Tensor<Precision, 4> f2r(Eigen::Tensor<std::complex<Precision>, 4> data_f) {
                typedef std::complex<Precision> fft_complex_p;

                int N_out = dims_(0,0) * dims_(1,0) * dims_(2,0);
                int N_in = dims_(0,0) * dims_(1,0) * (dims_(2,0)/2+1);

                Precision* out = (Precision*)fftw_malloc(sizeof(Precision) * N_out);
                fftw_complex_p* in = (fftw_complex_p*)fftw_malloc(sizeof(fft_complex_p) * N_in);               
                
                Eigen::Tensor<std::complex<std::complex<Precision>, 4> data_r(6,dims_(0,0),dims_(1,0),dims_(2,0));

                fftw_plan plan = fftw_plan_dft_c2r_3d(dims_(0,0), dims_(1,0), dims_(2,0),\
                 in, out, FFTW_ESTIMATE);

                for (int i=0; i<6; i++) {
                    for (int l=0; l<(dims_(2,0)/2 + 1); l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; i<dims_(0,0); j++) {
                                in[j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l] = data_f(i,j,k,l);
                            }
                        }
                    }

                    fftw_execute(plan);

                    for (int l=0; l<dims_(2,0); l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; i<dims_(0,0); j++) {
                                data_r(i,j,k,l) = out[j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l];
                            }
                        }
                    }
                }

                fftw_destroy_plan(plan);
                fftw_free(in);
                fftw_free(out);

                return data_r;
            }
            
        public:
            //do the iteration and compute stress_ and strain_
            void iteration() {
                Eigen::Array<Precision, 6, 1> tau()

                int count = 0;
                Precision err;
                err = error(stress_);

                while (err > tol && count < maxit_) {
                    
                }

                if (count == maxit_)


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
#ifndef MICROMECHANICS_HPP_
#define MICROMECHANICS_HPP_

#include <complex>
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "fftw3.h"


namespace mme {

    template <typename Precision>
    class micromechanics {
        public:
            Eigen::Array<int,3, 1> dims_;
            
            Eigen::Array<Precision, 6, 1> strain_0;
            Eigen::Tensor<int, 3> mat_;
            Eigen::Array<Precision, 2, Eigen::Dynamic> c_;
            Precision lamda_ref;
            Precision mu_ref;
            Eigen::Array<Precision, 3, 1> prds_;
            
            Precision tol_;
            int maxit_;

            Eigen::Tensor<Precision, 4> strain_;
            Eigen::Tensor<Precision, 4> stress_;
            
            bool reach_maxit = false;
        public:
            micromechanics() = default;

            micromechanics(Eigen::Array<Precision, 6, 1> E, Eigen::Tensor<int, 3> mat, Eigen::Array<Precision, 2, Eigen::Dynamic> c, \
            Eigen::Array<Precision, 3, 1> prds, Precision tol, int maxit)
            : strain_0(E), mat_(mat), c_(c), prds_(prds), tol_(tol), maxit_(maxit) 
            {
                //initial dims_
                for (int i=0; i<3; i++) {
                    dims_(i) = mat.dimension(i);
                }

                //initial lamda_ref and mu_ref
                // c is a matrix, the first row is lamda, the second row is mu
                auto max_c = c_.rowwise().maxCoeff();
                auto min_c = c_.rowwise().minCoeff();
                lamda_ref = (max_c(0) + min_c(0))/2;
                mu_ref = (max_c(1) + min_c(1))/2;

                //initial strain_ and stress_
                Eigen::Tensor<Precision,4> eps(6,dims_(0),dims_(1),dims_(2));
                Eigen::Tensor<Precision,4> sig(6,dims_(0),dims_(1),dims_(2));

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
                                array2tensor4d(sig, stressCompute(E, c(0, mat(i , j, k)), c(1, mat(i , j, k))), i, j, k);
                        }
                    }
                }

                stress_ = sig;
            }

            ~micromechanics() {}

        template<typename P>
            Eigen::Array<P, 6, 1> tensor4d2array(const Eigen::Tensor<P, 4, Eigen::ColMajor>& t, int ind1, int ind2, int ind3) const{
                Eigen::Array<P, 6, 1> a;
                for (int i=0; i<6; i++) {
                    a(i) = t(i, ind1, ind2, ind3);
                }

                return a;
            }

            template<typename P>
            void array2tensor4d(Eigen::Tensor<P, 4, Eigen::ColMajor>& t, Eigen::Array<P, 6, 1> a, int ind1, int ind2, int ind3) {
                for (int i=0; i<6; i++) {
                    t(i, ind1, ind2, ind3) = a(i);
                }
            }

            Eigen::Array<Precision, 6, 1> stressCompute(const Eigen::Array<Precision, 6, 1>& epsVec, Precision lamda, Precision mu) const{
                Precision sum =  epsVec[0] + epsVec[1] + epsVec[2];
                Eigen::Array<Precision, 6, 1> sigVec;
                sigVec = 2*mu*epsVec;
                for (int i = 0; i<3; i++) {
                    sigVec(i,0) += lamda*sum;
                }

                return sigVec;
            }

            Precision error(const Eigen::Tensor<Precision,4, Eigen::ColMajor>& sig) const{
                Precision err;
                Eigen::Tensor<std::complex<Precision>, 1> wave_ten(3);
                Eigen::Tensor<std::complex<Precision>, 4> sig_freq(sig.dimension(0), sig.dimension(1), sig.dimension(2), sig.dimension(3));
                Eigen::Tensor<std::complex<Precision>, 2> sig_ten(3,3);
                Eigen::Array<Precision,3, 1> wave_vec;
                Eigen::Array<int,3, 1> inds;
                Eigen::Tensor<std::complex<Precision>, 1> stress_dot_wave(3);

                sig_freq = r2f(sig);
                err = 0;
                for (int k=0; k<dims_(2)/2+1; k++) {
                    for (int j=0; j<dims_(1); j++) {
                        for (int i=0; i<dims_(0); i++) {
                            inds(0) = i;
                            inds(1) = j;
                            inds(2) = k;
                            wave_vec = waveVec(inds);
                            wave_ten(0) = wave_vec(0);
                            wave_ten(1) = wave_vec(1);
                            wave_ten(2) = wave_vec(2);
                            sig_ten = vec2ten(tensor4d2array(sig_freq, i, j, k));
                            stress_dot_wave = sig_ten.contract(wave_ten, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(1, 0)});
                            err += std::norm(stress_dot_wave(0)) + std::norm(stress_dot_wave(1)) + std::norm(stress_dot_wave(2));
                        }
                    }
                }

                err /= dims_(0) * dims_(1) *(dims_(2)/2+1);
                err = std::sqrt(err);
                sig_ten = vec2ten(tensor4d2array(sig_freq, 0, 0, 0)); 
                Precision sig0 = 0;
                for (int j=0; j<3; j++) {
                    for (int i=0; i<3; i++) {
                        sig0 += std::norm(sig_ten(i, j));
                    }
                }
                err /= sig0;
                
                return err;
            }

            Eigen::Array<Precision,3, 1> waveVec(const Eigen::Array<int, 3, 1>& inds) const{
                Eigen::Array<Precision,3, 1> xi;
                for(int i=0; i<3; i++) {
                    int bound = dims_[i] % 2 == 0 ? dims_[i] / 2 - 1 : dims_[i] / 2;
                    if (inds(i) <= bound) {
                        xi(i) = inds(i) /  prds_(i);
                    } 
                    else {
                        xi(i) = (-dims_(i) + inds(i)) / prds_(i);
                    }
                }

                return xi;
            }

            Eigen::Tensor<Precision, 4> greenOp(const Eigen::Array<Precision, 3, 1>& xi) const{
                Eigen::Tensor<Precision, 4> gam(3,3,3,3);
                gam.setZero();
                Precision coef_a, coef_b, xi_square;
                Eigen::Tensor<Precision, 4> temp;
                xi_square = xi(0,0) * xi(0,0) + xi(1,0) * xi(1,0) + xi(2,0) * xi(2,0);
                coef_a = 4 * mu_ref * xi_square;
                coef_b = lamda_ref + mu_ref;
                coef_b /= (mu_ref * (lamda_ref + 2*mu_ref)) * xi_square * xi_square;

                for (int j=0 ; j<3; j++) {
                    for (int i=0; i<3; i++) {
                        for (int h=0; h<3; h++) {
                            for (int k=0; k<3; k++) {
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
                                gam(k ,h, i ,j) /= coef_a; 
                                gam(k, h ,i, j) -= coef_b * xi(h, 0) * xi(j, 0) * xi(k, 0) * xi(i, 0);
                            }
                        }
                    }
                }
                return gam;
            }

            Eigen::Array<Precision, 6, 1> polarization(const Eigen::Array<Precision, 6, 1>& sigVec, const Eigen::Array<Precision, 6, 1>& epsVec) const{
                Eigen::Array<Precision, 6, 1> tauVec, c_eps;
                c_eps = stressCompute(epsVec, lamda_ref, mu_ref);
                tauVec = sigVec - c_eps;
                return tauVec;
            }

            template<typename P>
            Eigen::Array<P, 6, 1> ten2vec(const Eigen::Tensor<P, 2, Eigen::ColMajor>& T) const{
                Eigen::Array<P, 6, 1> V;
                V(0,0) = T(0,0);
                V(1,0) = T(1,1);
                V(2,0) = T(2,2);
                V(3,0) = T(1,2); //yz
                V(4,0) = T(0,2); //xz
                V(5,0) = T(0,1); //xy
                
                return V; 
            }

            template<typename P>
            Eigen::Tensor<P, 2> vec2ten(const Eigen::Array<P, 6, 1>& V) const{
                Eigen::Tensor<P, 2> T(3,3);
                for (int i=0; i<3; i++) {
                    T(i,i) = V(i,0);
                }

                T(0,1) = V(5,0);
                T(1,0) = V(5,0);
                T(0,2) = V(4,0);
                T(2,0) = V(4,0);
                T(1,2) = V(3,0);
                T(2,1) = V(3,0);

                return T;
            }


            //real2freq
            Eigen::Tensor<std::complex<Precision>, 4> r2f(const Eigen::Tensor<Precision, 4, Eigen::ColMajor>& data_r) const {
                int N_in = dims_(0,0) * dims_(1,0) * dims_(2,0);
                int N_out = dims_(0,0) * dims_(1,0) * (dims_(2,0)/2+1);

                Precision* in = (Precision*)fftw_malloc(sizeof(Precision) * N_in);
                std::complex<Precision>* out = (std::complex<Precision>*)fftw_malloc(sizeof(std::complex<Precision>) * N_out);
                
                Eigen::Tensor<std::complex<Precision>, 4> data_f(6,dims_(0,0),dims_(1,0),dims_(2,0)/2+1);

                fftw_plan plan = fftw_plan_dft_r2c_3d(dims_(0,0), dims_(1,0), dims_(2,0),\
                 in, reinterpret_cast<fftw_complex*>(out), FFTW_ESTIMATE);

                //eigen col-major, fftw row-major => j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l
                for (int i=0; i<6; i++) {
                    for (int l=0; l<dims_(2,0); l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; j<dims_(0,0); j++) {
                                in[j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l] = data_r(i,j,k,l);
                            }
                        }
                    }

                    fftw_execute(plan);

                    for (int l=0; l<dims_(2,0)/2 + 1; l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; j<dims_(0,0); j++) {
                                data_f(i,j,k,l) = out[j * dims_(1,0) *(dims_(2,0)/2+1) + k*(dims_(2,0)/2+1) + l];
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
            Eigen::Tensor<Precision, 4> f2r(const Eigen::Tensor<std::complex<Precision>, 4, Eigen::ColMajor>& data_f) const {
                int N_out = dims_(0,0) * dims_(1,0) * dims_(2,0);
                int N_in = dims_(0,0) * dims_(1,0) * (dims_(2,0)/2+1);

                Precision* out = (Precision*)fftw_malloc(sizeof(Precision) * N_out);
                std::complex<Precision>* in = (std::complex<Precision>*)fftw_malloc(sizeof(std::complex<Precision>) * N_in);               
                
                Eigen::Tensor<Precision, 4> data_r(6,dims_(0,0),dims_(1,0),dims_(2,0));

                fftw_plan plan = fftw_plan_dft_c2r_3d(dims_(0,0), dims_(1,0), dims_(2,0),\
                 reinterpret_cast<fftw_complex*>(in), out, FFTW_ESTIMATE);

                for (int i=0; i<6; i++) {
                    for (int l=0; l<dims_(2,0)/2 + 1; l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; j<dims_(0,0); j++) {
                                in[j * dims_(1,0) *(dims_(2,0)/2+1) + k*(dims_(2,0)/2+1) + l] = data_f(i,j,k,l);
                            }
                        }
                    }

                    fftw_execute(plan);

                    for (int l=0; l<dims_(2,0); l++) {
                        for (int k=0; k<dims_(1,0); k++) {
                            for (int j=0; j<dims_(0,0); j++) {
                                data_r(i,j,k,l) = out[j * dims_(1,0) *dims_(2,0) + k*dims_(2,0) + l] / N_out;
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
                Eigen::Tensor<Precision, 4> tau(6,dims_(0,0),dims_(1,0),dims_(2,0));
                Eigen::Tensor<std::complex<Precision>, 4> tau_f(6,dims_(0,0),dims_(1,0),dims_(2,0)/2 + 1);
                Eigen::Tensor<std::complex<Precision>, 4> eps_f(6,dims_(0,0),dims_(1,0),dims_(2,0)/2 + 1);
                Eigen::Tensor<Precision, 4> gam(2,2,2,2);
                Eigen::Tensor<std::complex<Precision>, 2> temp(3,3);
                Eigen::Array<Precision, 3, 1> xi;
                Eigen::Array<int, 3, 1> inds;
                Eigen::Array<std::complex<Precision>, 6, 1> scaled_strain_0;

                scaled_strain_0.setZero();
                for (int i = 0; i < 6; i++) {
                    scaled_strain_0(i).real(strain_0(i)*dims_(0)*dims_(1)*dims_(2)); 
                }

                int count = 0;
                Precision err;
                err = error(stress_);

                while (err > tol_ && count < maxit_) {
                    for (int k=0; k < dims_(2,0); k++) {
                        for (int j=0; j < dims_(1,0); j++) {
                            for (int i=0; i < dims_(0,0); i++) {
                                array2tensor4d(tau, polarization(tensor4d2array(stress_, i, j, k), tensor4d2array(strain_, i, j, k)), i, j, k);
                            }
                        }
                    }
                    // tau = stress_ -  c0_eps;
                    tau_f = r2f(tau);
                    for (int k=0; k < dims_(2,0)/2 + 1; k++) {
                        for (int j=0; j < dims_(1,0); j++) {
                            for (int i=0; i < dims_(0,0); i++) {
                                inds(0) = i;
                                inds(1) = j;
                                inds(2) = k;
                                xi = waveVec(inds);
                                gam = greenOp(xi);
                                temp = -gam.contract(vec2ten(tensor4d2array(tau_f, i, j, k)), Eigen::array<Eigen::IndexPair<int>, 2>{Eigen::IndexPair<int>(2, 0), Eigen::IndexPair<int>(3, 1)});
                                array2tensor4d(eps_f, ten2vec(temp), i, j, k); 
                            }
                        }
                    }
                    array2tensor4d(eps_f, scaled_strain_0, 0, 0, 0);

                    strain_ = f2r(eps_f);

                    for (int k = 0; k < dims_(2); k++) {
                        for (int j = 0; j < dims_(1); j++) {
                            for (int i = 0; i < dims_(0); i++) {
                                array2tensor4d(stress_, stressCompute(tensor4d2array(strain_, i, j, k), c_(0, mat_(i , j, k)), c_(1, mat_(i , j, k))), i, j, k);
                            }
                        }
                    }

                    count++;
                    err = error(stress_);
                }

                if (count >= maxit_) {
                    reach_maxit = true;
                }
            }


            Eigen::Tensor<Precision, 4> getStress() const {
                return stress_;
            }

            Eigen::Tensor<Precision, 4> getStrain() const {
                return strain_;
            }

    };

};

#endif
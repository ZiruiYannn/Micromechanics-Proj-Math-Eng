#include <complex>
#include "gtest/gtest.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "unsupported/Eigen/CXX11/Tensor"
#include "fftw3.h"
#include "micromechanics.hpp"

class FunctionTest: public ::testing::Test {
    protected:
    void SetUp() override {
        E << 1.0, -0.5, -0.5, 0.0, 0.0, 0.0;
        mat = Eigen::Tensor<int, 3>(2, 3, 4);
        mat.setConstant(0);
        c << 1.0, 3.0, 
             2.0, 4.0;
        prds << 2.0, 3.0, 4.0;
        tol = 1e-4;
        maxit = 30;
    }

    Eigen::Array<double, 6, 1> E;
    Eigen::Tensor<int, 3> mat;
    Eigen::Array<double, 2, 2> c;
    Eigen::Array<double, 3, 1> prds;
    double tol;
    int maxit;
};

TEST_F(FunctionTest, Initialization) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    ASSERT_EQ(m.dims_.size(), 3);
    ASSERT_EQ(m.dims_(0), 2);
    ASSERT_EQ(m.dims_(1), 3);
    ASSERT_EQ(m.dims_(2), 4);

    ASSERT_EQ(m.strain_.dimension(0), 6);
    ASSERT_EQ(m.strain_.dimension(1), 2);
    ASSERT_EQ(m.strain_.dimension(2), 3);
    ASSERT_EQ(m.strain_.dimension(3), 4);
    Eigen::ArrayXd epsVec(6);
    for (int i = 0; i < 6; i++) {
        epsVec(i) = m.strain_(i, 1, 2, 3);
    }
    EXPECT_TRUE(m.strain_(epsVec.isApprox(E)));

    ASSERT_EQ(m.stress_.dimension(0), 6);
    ASSERT_EQ(m.stress_.dimension(1), 2);
    ASSERT_EQ(m.stress_.dimension(2), 3);
    ASSERT_EQ(m.stress_.dimension(3), 4);

    EXPECT_DOUBLE_EQ(m.lamda_ref, 2.0);
    EXPECT_DOUBLE_EQ(m.mu_ref, 3.0);
}

TEST_F(FunctionTest, tensor4d2array) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Tensor<double, 4> t(6, 2, 2, 2);
    t(0, 0, 0, 0) = 1.0; t(0, 1, 0, 0) = 6.0;
    t(1, 0, 0, 0) = 2.0; t(1, 1, 0, 0) = 5.0;
    t(2, 0, 0, 0) = 3.0; t(2, 1, 0, 0) = 4.0;
    t(3, 0, 0, 0) = 4.0; t(3, 1, 0, 0) = 3.0;
    t(4, 0, 0, 0) = 5.0; t(4, 1, 0, 0) = 2.0;
    t(5, 0, 0, 0) = 6.0; t(5, 1, 0, 0) = 1.0;

    Eigen::Array<double, 6, 1> v = m.tensor4d2array(t, 0, 0, 0);

    EXPECT_DOUBLE_EQ(v(0), 1.0);
    EXPECT_DOUBLE_EQ(v(1), 2.0);
    EXPECT_DOUBLE_EQ(v(2), 3.0);
    EXPECT_DOUBLE_EQ(v(3), 4.0);
    EXPECT_DOUBLE_EQ(v(4), 5.0);
    EXPECT_DOUBLE_EQ(v(5), 6.0);

    v = m.tensor4d2array(t, 1, 0, 0);

    EXPECT_DOUBLE_EQ(v(0), 6.0);
    EXPECT_DOUBLE_EQ(v(1), 5.0);
    EXPECT_DOUBLE_EQ(v(2), 4.0);
    EXPECT_DOUBLE_EQ(v(3), 3.0);
    EXPECT_DOUBLE_EQ(v(4), 2.0);
    EXPECT_DOUBLE_EQ(v(5), 1.0);
}

TEST_F(FunctionTest, array2tensor4d) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Tensor<double, 4> t(6, 2, 2, 2);
    t.setConstant(0.0);

    Eigen::Array<double, 6, 1> v;
    v << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    m.array2tensor4d(t, v, 0, 0, 0);

    EXPECT_DOUBLE_EQ(t(0, 0, 0, 0), 1.0);
    EXPECT_DOUBLE_EQ(t(1, 0, 0, 0), 2.0);
    EXPECT_DOUBLE_EQ(t(2, 0, 0, 0), 3.0);
    EXPECT_DOUBLE_EQ(t(3, 0, 0, 0), 4.0);
    EXPECT_DOUBLE_EQ(t(4, 0, 0, 0), 5.0);
    EXPECT_DOUBLE_EQ(t(5, 0, 0, 0), 6.0);

    v << 6.0, 5.0, 4.0, 3.0, 2.0, 1.0;

    m.array2tensor4d(t, v, 0, 1, 1);

    EXPECT_DOUBLE_EQ(t(0, 0, 1, 1), 6.0);
    EXPECT_DOUBLE_EQ(t(1, 0, 1, 1), 5.0);
    EXPECT_DOUBLE_EQ(t(2, 0, 1, 1), 4.0);
    EXPECT_DOUBLE_EQ(t(3, 0, 1, 1), 3.0);
    EXPECT_DOUBLE_EQ(t(4, 0, 1, 1), 2.0);
    EXPECT_DOUBLE_EQ(t(5, 0, 1, 1), 1.0);
}

TEST_F(FunctionTest, stressCompute) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Array<double, 6, 1> strain;
    strain << 1.0, 2.0, 3.0, 0.75, 0.5, 0.25;
    double lamda = 3.0;
    double mu = 4.0;
    Eigen::Array<double, 6, 1> stress = m.stressCompute(strain, lamda, mu);
    
    EXPECT_DOUBLE_EQ(stress(0), 26.0);
    EXPECT_DOUBLE_EQ(stress(1), 34.0);
    EXPECT_DOUBLE_EQ(stress(2), 42.0);
    EXPECT_DOUBLE_EQ(stress(3), 6.0);
    EXPECT_DOUBLE_EQ(stress(4), 4.0);
    EXPECT_DOUBLE_EQ(stress(5), 2.0);
}

TEST_F(FunctionTest, error) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    m.stress_.setConstant(1.0);
    for (int i = 0; i < 6; i++) {
        m.stress_(i, 0, 0, 0) = 2.0;
    }

    double e = m.error(m.stress_);

    EXPECT_NEAR(e, 0.00015108205966200843, 1e-5);
}

TEST_F(FunctionTest, waveVec) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Array<int, 3, 1> inds;
    inds << 0, 0, 0;
    
    Eigen::Array<double, 3, 1> xi = m.waveVec(inds);

    EXPECT_DOUBLE_EQ(xi(0), 0.0);
    EXPECT_DOUBLE_EQ(xi(1), 0.0);
    EXPECT_DOUBLE_EQ(xi(2), 0.0);

    inds << 1, 1, 1;

    xi = m.waveVec(inds);

    EXPECT_DOUBLE_EQ(xi(0), -0.5);
    EXPECT_DOUBLE_EQ(xi(1), 1.0/3.0);
    EXPECT_DOUBLE_EQ(xi(2), 0.25);

    inds << 1, 2, 2;
    
    xi = m.waveVec(inds);

    EXPECT_DOUBLE_EQ(xi(1), -1.0/3.0);
    EXPECT_DOUBLE_EQ(xi(2), -0.5);
}

TEST_F(FunctionTest, greenOp) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Array<double, 3, 1> xi;
    xi << -0.5, 1.0/3.0, 0.25;

    Eigen::Tensor<double, 4> gam = m.greenOp(xi);

    ASSERT_EQ(gam.dimension(0), 3);
    ASSERT_EQ(gam.dimension(1), 3);
    ASSERT_EQ(gam.dimension(2), 3);
    ASSERT_EQ(gam.dimension(3), 3);

    EXPECT_NEAR(gam(0, 0, 0, 0), 0.12416017199677509, 1e-5);
    EXPECT_NEAR(gam(1, 1, 1, 1), 0.07309862940069875, 1e-5);
    EXPECT_NEAR(gam(2, 2, 2, 2), 0.0446452566514378, 1e-5);
    EXPECT_NEAR(gam(0, 1, 2, 0), -0.007793603869927439, 1e-5);

    Eigen::Tensor<double, 4> temp = -gam;

    EXPECT_NEAR(temp(0, 0, 0, 0), -0.12416017199677509, 1e-5);
}

TEST_F(FunctionTest, polarization) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Array<double, 6, 1> stress = m.stressCompute(m.strain_0, m.c_(0, 0), m.c_(1, 0));
    Eigen::Array<double, 6, 1> tau = m.polarization(stress, m.strain_0);

    EXPECT_DOUBLE_EQ(tau(0), -2.0);
    EXPECT_DOUBLE_EQ(tau(1), 1.0);
    EXPECT_DOUBLE_EQ(tau(2), 1.0);
    EXPECT_DOUBLE_EQ(tau(3), 0.0);
    EXPECT_DOUBLE_EQ(tau(4), 0.0);
    EXPECT_DOUBLE_EQ(tau(5), 0.0);
}

TEST_F(FunctionTest, ten2vec) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Tensor<double, 2> t(3, 3);
    t(0, 0) = 1.0;
    t(1, 1) = 2.0;
    t(2, 2) = 3.0;
    t(1, 2) = 4.0; t(2, 1) = 4.0;
    t(0, 2) = 5.0; t(2, 0) = 5.0;
    t(0, 1) = 6.0; t(1, 0) = 6.0;

    Eigen::Array<double, 6, 1> v = m.ten2vec(t);

    EXPECT_DOUBLE_EQ(v(0), 1.0);
    EXPECT_DOUBLE_EQ(v(1), 2.0);
    EXPECT_DOUBLE_EQ(v(2), 3.0);
    EXPECT_DOUBLE_EQ(v(3), 4.0);
    EXPECT_DOUBLE_EQ(v(4), 5.0);
    EXPECT_DOUBLE_EQ(v(5), 6.0);
}

TEST_F(FunctionTest, vec2ten) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Array<double, 6, 1> v;
    v << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    Eigen::Tensor<double, 2> t = m.vec2ten(v);

    EXPECT_DOUBLE_EQ(t(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(t(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(t(2, 2), 3.0);
    EXPECT_DOUBLE_EQ(t(1, 2), 4.0); EXPECT_DOUBLE_EQ(t(2, 1), 4.0);
    EXPECT_DOUBLE_EQ(t(0, 2), 5.0); EXPECT_DOUBLE_EQ(t(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(t(0, 1), 6.0); EXPECT_DOUBLE_EQ(t(1, 0), 6.0);
}

TEST_F(FunctionTest, r2f) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    for (int k = 2; k < 4; k++) {
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 2; i++) {
                m.strain_(0, i, j, k) = 0.0;
                m.strain_(1, i, j, k) = 0.0;
                m.strain_(2, i, j, k) = 0.0;
                m.strain_(3, i, j, k) = 0.2;
                m.strain_(4, i, j, k) = 0.3;
                m.strain_(5, i, j, k) = 0.4;
            }
        }
    }

    Eigen::Tensor<std::complex<double>, 4> strain_freq(6, 2, 3, 3);
    strain_freq = m.r2f(m.strain_);

    EXPECT_DOUBLE_EQ(std::real(strain_freq(0, 0, 0, 0)), 12.0); EXPECT_DOUBLE_EQ(std::imag(strain_freq(0, 0, 0, 0)), 0.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(1, 0, 0, 0)), -6.0); EXPECT_DOUBLE_EQ(std::imag(strain_freq(1, 0, 0, 0)), 0.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(2, 0, 0, 0)), -6.0); EXPECT_DOUBLE_EQ(std::imag(strain_freq(2, 0, 0, 0)), 0.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(3, 0, 0, 0)), 2.4); EXPECT_DOUBLE_EQ(std::imag(strain_freq(3, 0, 0, 0)), 0.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(4, 0, 0, 0)), 3.6); EXPECT_DOUBLE_EQ(std::imag(strain_freq(4, 0, 0, 0)), 0.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(5, 0, 0, 0)), 4.8); EXPECT_DOUBLE_EQ(std::imag(strain_freq(5, 0, 0, 0)), 0.0);

    EXPECT_DOUBLE_EQ(std::real(strain_freq(0, 0, 0, 1)), 6.0); EXPECT_DOUBLE_EQ(std::imag(strain_freq(0, 0, 0, 1)), -6.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(1, 0, 0, 1)), -3.0); EXPECT_DOUBLE_EQ(std::imag(strain_freq(1, 0, 0, 1)), 3.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(2, 0, 0, 1)), -3.0); EXPECT_DOUBLE_EQ(std::imag(strain_freq(2, 0, 0, 1)), 3.0);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(3, 0, 0, 1)), -1.2); EXPECT_DOUBLE_EQ(std::imag(strain_freq(3, 0, 0, 1)), 1.2);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(4, 0, 0, 1)), -1.8); EXPECT_DOUBLE_EQ(std::imag(strain_freq(4, 0, 0, 1)), 1.8);
    EXPECT_DOUBLE_EQ(std::real(strain_freq(5, 0, 0, 1)), -2.4); EXPECT_DOUBLE_EQ(std::imag(strain_freq(5, 0, 0, 1)), 2.4);
}

TEST_F(FunctionTest, f2r) {
    mme::micromechanics<double> m(E, mat, c, prds, tol, maxit);

    Eigen::Tensor<std::complex<double>, 4> strain_freq(6, 2, 3, 3);
    strain_freq.setConstant(0.0);
    strain_freq(0, 0, 0, 0).real(-10.8);
    strain_freq(1, 0, 0, 0).real(-21.6);
    strain_freq(2, 0, 0, 0).real(-32.4);
    strain_freq(3, 0, 0, 0).real(-43.2);
    strain_freq(4, 0, 0, 0).real(-54.0);
    strain_freq(5, 0, 0, 0).real(-64.8);
    strain_freq(0, 0, 0, 1).real(6.6); strain_freq(0, 0, 0, 1).imag(-6.6);
    strain_freq(1, 0, 0, 1).real(13.2); strain_freq(1, 0, 0, 1).imag(-13.2);
    strain_freq(2, 0, 0, 1).real(19.8); strain_freq(2, 0, 0, 1).imag(-19.8);
    strain_freq(3, 0, 0, 1).real(26.4); strain_freq(3, 0, 0, 1).imag(-26.4);
    strain_freq(4, 0, 0, 1).real(33.0); strain_freq(4, 0, 0, 1).imag(-33.0);
    strain_freq(5, 0, 0, 1).real(39.6); strain_freq(5, 0, 0, 1).imag(-39.6);

    Eigen::Tensor<double, 4> strain(6, 2, 3, 4);
    strain = m.f2r(strain_freq);

    EXPECT_DOUBLE_EQ(strain(0, 0, 0, 0), 0.1);
    EXPECT_DOUBLE_EQ(strain(1, 0, 0, 0), 0.2);
    EXPECT_DOUBLE_EQ(strain(2, 0, 0, 0), 0.3);
    EXPECT_DOUBLE_EQ(strain(3, 0, 0, 0), 0.4);
    EXPECT_DOUBLE_EQ(strain(4, 0, 0, 0), 0.5);
    EXPECT_DOUBLE_EQ(strain(5, 0, 0, 0), 0.6);

    EXPECT_DOUBLE_EQ(strain(0, 0, 0, 2), -1.0);
    EXPECT_DOUBLE_EQ(strain(1, 0, 0, 2), -2.0);
    EXPECT_DOUBLE_EQ(strain(2, 0, 0, 2), -3.0);
    EXPECT_DOUBLE_EQ(strain(3, 0, 0, 2), -4.0);
    EXPECT_DOUBLE_EQ(strain(4, 0, 0, 2), -5.0);
    EXPECT_DOUBLE_EQ(strain(5, 0, 0, 2), -6.0);
}



class AlgorithmTest: public ::testing::Test {
    protected:
    void SetUp() override {
        lambda1 = mod1*nu / ((1+nu)*(1-2*nu));
        mu1 = mod1 / (2*(1+nu));

        lambda2 = mod2*nu / ((1+nu)*(1-2*nu));
        mu2 = mod2 / (2*(1+nu));

        lambda3 = mod3*nu / ((1+nu)*(1-2*nu));
        mu3 = mod3 / (2*(1+nu));

        Ex << E0, -nu*E0, -nu*E0, 0.0, 0.0, 0.0;
        Ez << -nu*E0, -nu*E0, E0, 0.0, 0.0, 0.0;
    }

    double avg(Eigen::Tensor<double, 4> const & field, int const & idx, Eigen::Array<int, 2, 3> const & range) {
        double sum = 0.0;

        for (int k = range(0, 0); k < range(1, 0); k++) {
            for (int j = range(0, 1); j < range(1, 1); j++) {
                for (int i = range(0, 2); k < range(1, 2); i++) {
                    sum += field(idx, i, j, k);
                }
            }
        }

        return sum / (field.dimension(0)*field.dimension(1)*field.dimension(2));
    }

    double nu = 0.28;
    double mod1 = 210e9;
    double mod2 = 0.4*mod1;
    double mod3 = 0.6*mod1;
    double lambda1;
    double mu1;
    double lambda2;
    double mu2;
    double lambda3;
    double mu3;
    double E0 = 0.01;
    Eigen::Array<double, 6, 1> Ex;
    Eigen::Array<double, 6, 1> Ez;
    double tol = 1e-4;
    int maxit = 30;
};

TEST_F(AlgorithmTest, avg) {
    Eigen::Tensor<double, 4> t(3, 3, 3, 3);
    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
                t(0, i, j, k) = 1;
                t(1, i, j, k) = k+1;
                t(2, i, j, k) = 0;
            }
        }
    }
    Eigen::Array<int, 2, 3> range;
    range << 0, 0, 0,
             3, 3, 3;
    
    EXPECT_DOUBLE_EQ(avg(t, 0, range), 1.0);
    EXPECT_DOUBLE_EQ(avg(t, 1, range), 2.0);
    EXPECT_DOUBLE_EQ(avg(t, 2, range), 0.0);

    range << 0, 0, 1,
             3, 3, 3;

    EXPECT_DOUBLE_EQ(avg(t, 1, range), 2.5);
}

TEST_F(AlgorithmTest, Homogenous) {
    Eigen::Tensor<int, 3> mat(10, 10, 10);
    mat.setConstant(0);
    Eigen::Array<double, 2, 1> c;
    c << lambda1,
         mu1;
    Eigen::Array<double, 3, 1> prds;
    prds << 10.0, 10.0, 10.0;

    mme::micromechanics<double> m(Ez, mat, c, prds, tol, maxit);
    m.iteration();

    Eigen::Array<int, 2, 3> range;
    range << 0, 0, 0,
             10, 10, 10;
    double avge = avg(m.getStrain(), 2, range);
    double avgs = avg(m.getStress(), 2, range);

    EXPECT_DOUBLE_EQ(avgs/avge, mod1);
}

TEST_F(AlgorithmTest, Series2) {
    Eigen::Tensor<int, 3> mat(10, 10, 10);
    mat.setConstant(0);
    for (int k = 5; k < 10; k++) {
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < 10; i++) {
                mat(i, j, k) = 1;
            }
        }
    }
    Eigen::Array<double, 2, 2> c;
    c << lambda1, lambda2,
         mu1, mu2;
    Eigen::Array<double, 3, 1> prds;
    prds << 10.0, 10.0, 10.0;

    mme::micromechanics<double> m(Ez, mat, c, prds, tol, maxit);
    m.iteration();

    Eigen::Array<int, 2, 3> range1;
    range1 << 0, 0, 0,
             10, 10, 5;
    Eigen::Array<int, 2, 3> range2;
    range2 << 0, 0, 5,
             10, 10, 10;
    double avg1 = avg(m.getStrain(), 2, range1);
    double avg2 = avg(m.getStrain(), 2, range2);

    EXPECT_NEAR(avg1/avg2, 0.5, 0.05);

    avg1 = avg(m.getStress(), 2, range1);
    avg2 = avg(m.getStress(), 2, range2);

    EXPECT_NEAR(avg2/avg1, 1.0, 1e-5);
}

TEST_F(AlgorithmTest, Parallel2) {
    Eigen::Tensor<int, 3> mat(10, 10, 10);
    mat.setConstant(0);
    for (int k = 5; k < 10; k++) {
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < 10; i++) {
                mat(i, j, k) = 1;
            }
        }
    }
    Eigen::Array<double, 2, 2> c;
    c << lambda1, lambda2,
         mu1, mu2;
    Eigen::Array<double, 3, 1> prds;
    prds << 10.0, 10.0, 10.0;

    mme::micromechanics<double> m(Ex, mat, c, prds, tol, maxit);
    m.iteration();

    Eigen::Array<int, 2, 3> range1;
    range1 << 0, 0, 0,
             10, 10, 5;
    Eigen::Array<int, 2, 3> range2;
    range2 << 0, 0, 5,
             10, 10, 10;
    double avg1 = avg(m.getStress(), 0, range1);
    double avg2 = avg(m.getStress(), 0, range2);

    EXPECT_NEAR(avg1/avg2, 0.4, 1e-5);

    avg1 = avg(m.getStrain(), 0, range1);
    avg2 = avg(m.getStrain(), 0, range2);

    EXPECT_NEAR(avg2/avg1, 1.0, 1e-5);
}

TEST_F(AlgorithmTest, Series3) {
    Eigen::Tensor<int, 3> mat(10, 10, 15);
    mat.setConstant(0);
    for (int k = 5; k < 10; k++) {
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < 10; i++) {
                mat(i, j, k) = 1;
            }
        }
    }
    for (int k = 10; k < 15; k++) {
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < 10; i++) {
                mat(i, j, k) = 2;
            }
        }
    }
    Eigen::Array<double, 2, 3> c;
    c << lambda1, lambda2, lambda3,
         mu1, mu2, mu3;
    Eigen::Array<double, 3, 1> prds;
    prds << 10.0, 10.0, 15.0;

    mme::micromechanics<double> m(Ez, mat, c, prds, tol, maxit);
    m.iteration();

    Eigen::Array<int, 2, 3> range1;
    range1 << 0, 0, 0,
             10, 10, 5;
    Eigen::Array<int, 2, 3> range2;
    range2 << 0, 0, 5,
             10, 10, 10;
    Eigen::Array<int, 2, 3> range3;
    range3 << 0, 0, 10,
             10, 10, 15;  
    double avg1 = avg(m.getStrain(), 2, range1);
    double avg2 = avg(m.getStrain(), 2, range2);
    double avg3 = avg(m.getStrain(), 2, range3);

    EXPECT_NEAR(avg1/avg2, 0.5, 0.05);
    EXPECT_NEAR(avg1/avg3, 0.6, 0.05);

    avg1 = avg(m.getStress(), 2, range1);
    avg2 = avg(m.getStress(), 2, range2);
    avg3 = avg(m.getStress(), 2, range3);

    EXPECT_NEAR(avg2/avg1, 1.0, 0.5);
    EXPECT_NEAR(avg3/avg1, 1.0, 0.5);
}

TEST_F(AlgorithmTest, Parallel3) {
    Eigen::Tensor<int, 3> mat(10, 10, 15);
    mat.setConstant(0);
    for (int k = 5; k < 10; k++) {
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < 10; i++) {
                mat(i, j, k) = 1;
            }
        }
    }
    for (int k = 10; k < 15; k++) {
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < 10; i++) {
                mat(i, j, k) = 2;
            }
        }
    }
    Eigen::Array<double, 2, 3> c;
    c << lambda1, lambda2, lambda3,
         mu1, mu2, mu3;
    Eigen::Array<double, 3, 1> prds;
    prds << 10.0, 10.0, 15.0;

    mme::micromechanics<double> m(Ex, mat, c, prds, tol, maxit);
    m.iteration();

    Eigen::Array<int, 2, 3> range1;
    range1 << 0, 0, 0,
             10, 10, 5;
    Eigen::Array<int, 2, 3> range2;
    range2 << 0, 0, 5,
             10, 10, 10;
    Eigen::Array<int, 2, 3> range3;
    range3 << 0, 0, 10,
             10, 10, 15;  
    double avg1 = avg(m.getStress(), 2, range1);
    double avg2 = avg(m.getStress(), 2, range2);
    double avg3 = avg(m.getStress(), 2, range3);

    EXPECT_NEAR(avg2/avg1, 0.4, 1e-5);
    EXPECT_NEAR(avg3/avg1, 0.6, 1e-5);

    avg1 = avg(m.getStrain(), 2, range1);
    avg2 = avg(m.getStrain(), 2, range2);
    avg3 = avg(m.getStrain(), 2, range3);

    EXPECT_NEAR(avg1/avg2, 1.0, 1e-5);
    EXPECT_NEAR(avg1/avg3, 1.0, 1e-5);
}




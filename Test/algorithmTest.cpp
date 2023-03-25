#include "gtest/gtest.h"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "micromechanics.hpp"

class FunctionTest: public ::testing::Test {
    protected:
    void SetUp() override {
        mat = Eigen::Tensor<int, 3>(2, 3, 4);
        E << 1.0, -0.5, -0.5, 0.0, 0.0, 0.0;
        mat.setConstant(0);
        c << 1.0, 3.0, 
             2.0, 4.0;
        prds << 2.0, 3.0, 4.0;
        tol = 1e-4;
        maxit = 50;
        m = new mme::micromechanics<double>(E, mat, c, prds, tol, maxit);
    }

    Eigen::Array<double, 6, 1> E;
    Eigen::Tensor<int, 3> mat;
    Eigen::Array<double, 2, 2> c;
    Eigen::Array<double, 3, 1> prds;
    double tol;
    int maxit;
    mme::micromechanics<double> m;
};

TEST_F(FunctionTest, Initialization) {
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

    // polarization

    // frequencies 

    EXPECT_DOUBLE_EQ(m.lamda_ref, 2.0);
    EXPECT_DOUBLE_EQ(m.mu_ref, 3.0);
}

TEST_F(FunctionTest, stressCompute) {
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
    m.stress_.setConstant(1.0);
    for (int i = 0; i < 6; i++) {
        m.stress_(i, 0, 0, 0) = 2.0;
    }

    double e = m.error(m.stress_);

    EXPECT_NEAR(e, 0.00015108205966200843, 1e-5);
}

TEST_F(FunctionTest, waveVec) {
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
}

TEST_F(FunctionTest, polarization) {
    Eigen::Array<double, 6, 1> stress = m.stressCompute(m.strain_0, m.c_(0, 0), m.c_(1, 0));
    Eigen::Array<double, 6, 1> tau = m.polarization(stress, m.strain_0);

    EXPECT_DOUBLE_EQ(tau(0), -2.0);
    EXPECT_DOUBLE_EQ(tau(1), 1.0);
    EXPECT_DOUBLE_EQ(tau(2), 1.0);
    EXPECT_DOUBLE_EQ(tau(3), 0.0);
    EXPECT_DOUBLE_EQ(tau(4), 0.0);
    EXPECT_DOUBLE_EQ(tau(6), 0.0);
}



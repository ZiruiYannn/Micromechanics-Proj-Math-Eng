#include "gtest/gtest.h"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "fftw3.h"
#include "micromechanics.hpp"

class AlgorithmTest: public ::testing::Test {
    protected:
    void SetUp() override {
        Eigen::ArrayXd E(6);
        E << 1.0, -0.5, -0.5, 0.0, 0.0, 0.0;
        Eigen::Tensor<int, 3> mat(2, 3, 4);
        mat.setConstant(0);
        Eigen::ArrayXXd c(2, 2);
        c << 1.0, 3.0, 
             2.0, 4.0;
        Eigen::ArrayXd prds(3);
        prds << 2.0, 3.0, 4.0;
        double tol = 1e-4;
        int maxit = 50;

        m = new mme::micromechanics<>(E, mat, c, prds, tol, maxit);
    }

    mme::micromechanics<double> m;
};

TEST_F(AlgorithmTest, Initialization) {
    ASSERT_EQ(m.dims_.size(), 3);
    ASSERT_EQ(m.dims_(0), 2);
    ASSERT_EQ(m.dims_(1), 3);
    ASSERT_EQ(m.dims_(2), 4);

    ASSERT_EQ(m.strain_.dimension(0), 6);
    ASSERT_EQ(m.strain_.dimension(1), 2);
    ASSERT_EQ(m.strain_.dimension(2), 3);
    ASSERT_EQ(m.strain_.dimension(3), 4);
    ASSERT_TRUE(m.strain_(all, last, last, last).isApprox(E));

    ASSERT_EQ(m.stress_.dimension(0), 6);
    ASSERT_EQ(m.stress_.dimension(1), 2);
    ASSERT_EQ(m.stress_.dimension(2), 3);
    ASSERT_EQ(m.stress_.dimension(3), 4);

    // polarization

    // frequencies 

    ASSERT_DOUBLE_EQ(m.lamda_ref, 2.0);
    ASSERT_DOUBLE_EQ(m.mu_ref, 3.0);
}

TEST_F(AlgorithmTest, stressCompute) {
    Eigen::ArrayXd eps(6);
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

TEST_F(AlgorithmTest, error) {
    m.stress_.setConstant(1.0);
    for (i = 0; i < 6; i++) {
        m.stress_(i, 0, 0, 0) = 2.0;
    }

    double e = m.error();

    EXPECT_NEAR(e, 0.00015108205966200843, 1e-5);
}

TEST_F(AlgorithmTest, waveVec) {
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

TEST_F(AlgorithmTest, greenOp) {
    // check if function greenOp works
}

TEST_F(AlgorithmTest, polarization) {
    // check if function polarization works
}



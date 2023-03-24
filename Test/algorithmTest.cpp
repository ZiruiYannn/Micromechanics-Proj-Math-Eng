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
        c << 1.0, 3.0, 2.0, 4.0;
        Eigen::ArrayXd prds(3);
        prds << 2.0, 3.0, 4.0;
        double tol = 1e-4;
        int maxit = 50;

        m = new mme::micromechanics<>(E, mat, c, prds, tol, maxit);
    }

    mme::micromechanics<double> m;
};

TEST_F(AlgorithmTest, Initialization) {
    // check if micromechanics instance is correctly initialized
}

TEST_F(AlgorithmTest, stressCompute) {
    // check if function stressCompute works
}

TEST_F(AlgorithmTest, error) {
    // check if function error works
}

TEST_F(AlgorithmTest, waveVec) {
    // check if function waveVec works
}

TEST_F(AlgorithmTest, greenOp) {
    // check if function greenOp works
}

TEST_F(AlgorithmTest, polarization) {
    // check if function polarization works
}



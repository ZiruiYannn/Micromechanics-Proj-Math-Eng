#include "gtest/gtest.h"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "IO.hpp"

class InTest: public ::testing::Test {
    protected:
    void SetUp() override {
        auto [c_in, mat_in, prds_in] = read_material("../Input/test_mat.in"); 
        c = c_in;
        mat = mat_in;
        prds = prds_in;
    }

    Eigen::ArrayXXd c;
    Eigen::Tensor<int, 3> mat;
    Eigen::ArrayXd prds;
};

TEST_F(InTest, Dimensions) {
    ASSERT_EQ(c.rows(), 2);
    ASSERT_EQ(c.cols(), 3);

    ASSERT_EQ(mat.rank(), 3);
    ASSERT_EQ(mat.dimension(0), 2);
    ASSERT_EQ(mat.dimension(1), 2);
    ASSERT_EQ(mat.dimension(2), 3);

    ASSERT_EQ(prds.size(), 3);
}

TEST_F(InTest, Contents) {
    double nu = 0.28;
    double mod = 210e9;
    double lambda = mod*nu / ((1+nu)*(1-2*nu));
    double mu = mod / (2*(1+nu)); 
    
    EXPECT_DOUBLE_EQ(c(0, 0), lambda);
    EXPECT_DOUBLE_EQ(c(1, 0), mu);

    EXPECT_EQ(mat(0, 0, 0), 0);
    EXPECT_EQ(mat(0, 0, 1), 1);
    EXPECT_EQ(mat(0, 0, 2), 2);

    EXPECT_DOUBLE_EQ(prds(0), 2.0);
    EXPECT_DOUBLE_EQ(prds(1), 2.0);
    EXPECT_DOUBLE_EQ(prds(2), 3.0);
}

TEST(InTest, strain) {
    Eigen::Array<double, 6, 1> eps = read_strain("../Input/test_e.in");

    EXPECT_DOUBLE_EQ(eps(0), 1.0);
    EXPECT_DOUBLE_EQ(eps(1), 0.0);
    EXPECT_DOUBLE_EQ(eps(2), 0.0);
    EXPECT_DOUBLE_EQ(eps(3), 0.0);
    EXPECT_DOUBLE_EQ(eps(4), 0.0);
    EXPECT_DOUBLE_EQ(eps(5), 0.0);
}

TEST(OutTest, Working) {
    Eigen::Array<double, 3, 1> prds;
    prds << 2.0, 2.0, 3.0;
    Eigen::Tensor<double, 4> t(2, 2, 2, 3);
    t.setConstant(0);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            t(0, i, j, 1) = 1.0;
        }
    }
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            t(0, i, j, 2) = 2.0;
        }
    }

    write_material(prds, t, 0, "../Output/test.out");
}
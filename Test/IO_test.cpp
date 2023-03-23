#include "gtest/gtest.h"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "IO.hpp"

class InTest: public ::testing::Test {
    protected:
    void SetUp() override {
        auto [c_in, mat_in, prds_in] = read_material("../input/test.in"); 
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
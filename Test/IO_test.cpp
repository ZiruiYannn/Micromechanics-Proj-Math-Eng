#include "gtest/gtest.h"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "read_material.hpp"

TEST(InTest, BasicTest) {
    double nu = 0.28;
    double mod = 210e9;
    double lambda = mod*nu / ((1+nu)*(1-2*nu));
    double mu = mod / (2*(1+nu));

    auto [c, mat, prds] = read_material("input/test.in"); 
    
    // check if dimensions are as expected
    ASSERT_EQ(c.rows(), 1);
    ASSERT_EQ(c.cols(), 2);

    ASSERT_EQ(mat.rank(), 3);
    for (int i = 0; i < 3; i++){
        ASSERT_EQ(mat.dimension(i), 2);
    }

    ASSERT_EQ(prds.size(), 3);


    // check if contents are as expected
    EXPECT_DOUBLE_EQ(c(0, 0), lambda);
    EXPECT_DOUBLE_EQ(c(0, 1), mu);

    EXPECT_EQ(mat(1, 1, 1), 0);

    EXPECT_DOUBLE_EQ(prds(2), 2.0);
}
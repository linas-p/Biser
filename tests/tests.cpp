/* Copyright 2016
 * Linas Petkevicius
 * Vilnius University
 * GNU General Public License
 * */

#include <gtest/gtest.h>

#include <iostream>
#include <ctime>
#include <vector>
#include <string>

#include "BiserLikeModel/utils.h"
#define LARGE    1000000
#define EPSILON  1e-10
#define EPSILION 1e-2


using namespace BiserLikeModel;

TEST(MM2, usage) {
    /*Critical cases*/
    ASSERT_DOUBLE_EQ(MM2(1, 0), 0);
    ASSERT_DOUBLE_EQ(MM2(0, 1), 0);
    ASSERT_DOUBLE_EQ(MM2(0, 0), 0);

    /*Limit cases*/
    EXPECT_NEAR(MM2(LARGE, 1), 1, EPSILION);
    EXPECT_NEAR(MM2(-LARGE, 1), 1, EPSILION);
    EXPECT_NEAR(MM2(1, LARGE), 0.5, EPSILION);
    EXPECT_NEAR(MM2(1, -LARGE), 0.5, EPSILION);

    /*Usage*/
    EXPECT_NEAR(MM2(1, 1), 1./3., EPSILION);
}

TEST(MM, usage) {

    /*Critical cases*/
    ASSERT_DOUBLE_EQ(MM(0, 1, 1), 0);
    ASSERT_DOUBLE_EQ(MM(1, 0, 1), 0);

    /*Limit cases*/
    EXPECT_NEAR(MM(LARGE, 2, 1), 2., EPSILION);
    EXPECT_NEAR(MM(1, 1, LARGE), 0., EPSILION);

    /*Usage*/
    EXPECT_NEAR(MM(1, 1, 1), 1./2., EPSILION);
}


TEST(LaplacePolar0, usage) {

    /*Usage*/
    double array1[2] = {0, 1};
    EXPECT_NEAR(LaplacePolar0(array1, 1.), 2, EPSILION);

    double array2[2] = {1, 0};
    EXPECT_NEAR(LaplacePolar0(array2, 1.), -2, EPSILION);


}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}


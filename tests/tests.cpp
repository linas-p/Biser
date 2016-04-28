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
    ASSERT_DOUBLE_EQ(MM2(1, 0), 0);
    ASSERT_DOUBLE_EQ(MM2(0, 1), 0);
    ASSERT_DOUBLE_EQ(MM2(0, 0), 0);

    EXPECT_NEAR(MM2(LARGE, 1), 1, EPSILION);
    EXPECT_NEAR(MM2(-LARGE, 1), 1, EPSILION);
    EXPECT_NEAR(MM2(1, LARGE), 0.5, EPSILION);
    EXPECT_NEAR(MM2(1, -LARGE), 0.5, EPSILION);

    EXPECT_NEAR(MM2(1, 1), 1./3., EPSILION);
}



int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}


#include "test_fmm.h"
#include "gtest/gtest.h"

TEST(FMMTest, Accuracy) {
  EXPECT_GT(1e-3, test_fmm(10));
  EXPECT_GT(1e-6, test_fmm(20));
}

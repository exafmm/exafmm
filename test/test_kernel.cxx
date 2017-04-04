#include "test_kernel.h"
#include "gtest/gtest.h"

TEST(KernelTest, Accuracy) {
  EXPECT_GT(1e-3, test_kernel(10));
  EXPECT_GT(1e-6, test_kernel(20));
  EXPECT_GT(1e-9, test_kernel(30));
}

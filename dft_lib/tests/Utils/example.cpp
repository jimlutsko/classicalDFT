#include "dft_lib/utils/example.h"

#include <gtest/gtest.h>

TEST(utilsTest, testMultiplyByTwo) {
    //arrange
    //act
    //assert
    EXPECT_EQ(dft_core::MultiplyByTwo(0), 0);
    EXPECT_EQ(dft_core::MultiplyByTwo(10), 20);
    EXPECT_EQ(dft_core::MultiplyByTwo(50), 100);
}

TEST(utilsTest, testSalute)
{
    auto name = std::string("Miguel");
    ASSERT_EQ(dft_core::Salute(name), "Hello Miguel");
    ASSERT_NE(dft_core::Salute(name), "miguel");
}

TEST(utilsTest, testGetType)
{
    ASSERT_EQ(dft_core::GetType("asda"), "PKc");
    ASSERT_EQ(dft_core::GetType((double) 1.0), "d");
    ASSERT_EQ(dft_core::GetType(1), "i");
}
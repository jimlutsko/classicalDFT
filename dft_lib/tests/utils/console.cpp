#include <gtest/gtest.h>

#include "dft_lib/utils/console.h"

//region Methods
TEST(console, write_works_ok)
{
  testing::internal::CaptureStdout();

  console::write("test");
  std::string output = testing::internal::GetCapturedStdout();

  auto expected_str = std::string("test");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}

TEST(console, write_line_works_ok)
{
  testing::internal::CaptureStdout();

  console::write_line("test");
  std::string output = testing::internal::GetCapturedStdout();

  auto expected_str = std::string("test\n");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}

TEST(console, new_line_works_ok)
{
  testing::internal::CaptureStdout();

  console::new_line();
  std::string output = testing::internal::GetCapturedStdout();

  auto expected_str = std::string("\n");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}
//endregion


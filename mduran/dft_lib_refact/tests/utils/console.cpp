#include <gtest/gtest.h>

#include "dft_lib/utils/console.h"

//region Methods
TEST(console, write_works_ok)
{
  testing::internal::CaptureStdout();

  console::Write("test");
  std::string output = testing::internal::GetCapturedStdout();

  auto expected_str = std::string("test");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}

TEST(console, write_line_works_ok)
{
  testing::internal::CaptureStdout();

  console::WriteLine("test");
  std::string output = testing::internal::GetCapturedStdout();

  auto expected_str = std::string("test\n");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}

TEST(console, write_line_intializer_list_works_ok)
{
  testing::internal::CaptureStdout();

  console::WriteLine({"test", "one", "two"});
  std::string output = testing::internal::GetCapturedStdout();

  auto expected_str = std::string("test\none\ntwo\n");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}

TEST(console, new_line_works_ok)
{
  testing::internal::CaptureStdout();

  console::NewLine();
  std::string output = testing::internal::GetCapturedStdout();

  auto expected_str = std::string("\n");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}
//endregion


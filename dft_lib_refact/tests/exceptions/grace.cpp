#include <gtest/gtest.h>

#include "dft_lib/exceptions/grace_exception.h"

TEST(grace_exceptions, grace_exception_cttor_test)
{
  std::string msg = "new exception";
  auto exception = dft_core::exception::GraceException(msg);
  EXPECT_THROW(throw exception, dft_core::exception::GraceException);
  ASSERT_STREQ(exception.error_message().c_str(), msg.c_str());
  ASSERT_STREQ(exception.what(), msg.c_str());
}

TEST(grace_exceptions, grace_no_open_cttor_test)
{
  std::string msg = "No grace subprocess currently connected.";
  auto exception = dft_core::exception::GraceNotOpenedException();
  EXPECT_THROW(throw exception, dft_core::exception::GraceNotOpenedException);
  ASSERT_STREQ(exception.error_message().c_str(), msg.c_str());
  ASSERT_STREQ(exception.what(), msg.c_str());
}

TEST(grace_exceptions, grace_communication_failed_cttor_default_test)
{
  std::string msg = "There was a problem while communicating with Grace.";
  auto exception = dft_core::exception::GraceCommunicationFailedException();
  EXPECT_THROW(throw exception, dft_core::exception::GraceCommunicationFailedException);
  ASSERT_STREQ(exception.error_message().c_str(), msg.c_str());
  ASSERT_STREQ(exception.what(), msg.c_str());
}

TEST(grace_exceptions, grace_communication_failed_cttor_test)
{
  std::string msg = "example";
  auto exception = dft_core::exception::GraceCommunicationFailedException(msg);
  EXPECT_THROW(throw exception, dft_core::exception::GraceCommunicationFailedException);
  ASSERT_STREQ(exception.error_message().c_str(), msg.c_str());
  ASSERT_STREQ(exception.what(), msg.c_str());
}
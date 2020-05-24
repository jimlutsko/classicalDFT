#include "dft_lib/utils/config_parser.h"

#include <gtest/gtest.h>

#include <boost/property_tree/detail/info_parser_error.hpp>
#include <boost/property_tree/detail/xml_parser_error.hpp>
#include <boost/range/combine.hpp>

#include "dft_lib/utils/console.h"

using namespace dft_core;

//region Cttors:

TEST(config_parser, default_cttor_works_ok)
{
  std::string expected_file_path = "config.ini";
  auto expected_file_type = dft_core::config_parser::FileType::INI;

  auto config = dft_core::config_parser::ConfigParser();
  ASSERT_STREQ(config.config_file_path().c_str(), expected_file_path.c_str());
  ASSERT_EQ(config.config_file_type(), expected_file_type);
}

TEST(config_parser, specific_cttor_works_ok)
{
  std::vector<dft_core::config_parser::FileType> types {
      dft_core::config_parser::FileType::INI,
      dft_core::config_parser::FileType::JSON,
      dft_core::config_parser::FileType::XML,
      dft_core::config_parser::FileType::INFO
  };

  std::vector<std::string> files {
      "config.ini",
      "config.json",
      "config.xml",
      "config.info"
  };

  std::string f;
  dft_core::config_parser::FileType t;

  for (const auto& tuple : boost::combine(types, files)) {
    boost::tie(t, f) = tuple;
    auto config = dft_core::config_parser::ConfigParser(f, t);
    auto expected_file_type = t;
    auto expected_file_path = f;
    ASSERT_STREQ(config.config_file_path().c_str(), expected_file_path.c_str());
    ASSERT_EQ(config.config_file_type(), expected_file_type);
  }
}

TEST(config_parser, tree_works_ok)
{
  std::vector<dft_core::config_parser::FileType> types {
      dft_core::config_parser::FileType::INI,
      dft_core::config_parser::FileType::JSON,
      dft_core::config_parser::FileType::XML,
      dft_core::config_parser::FileType::INFO
  };

  std::vector<std::string> files {
      "config.ini",
      "config.json",
      "config.xml",
      "config.info"
  };

  std::string f;
  dft_core::config_parser::FileType t;

  for (const auto& tuple : boost::combine(types, files)) {
    boost::tie(t, f) = tuple;
    auto config = dft_core::config_parser::ConfigParser(f, t);
    auto expected_file_type = t;
    auto expected_file_path = f;

    double expected_double_value = 10;
    std::string expected_string_value = "a_text_string";

    std::string section = (config.config_file_path().find("info") == std::string::npos) ? "default." : "";
    auto actual_string_value = config.tree().get<std::string>(section + "StringValue");
    auto actual_double_value = config.tree().get<double>(section + "DoubleValue");

    ASSERT_DOUBLE_EQ(actual_double_value, expected_double_value);
    ASSERT_STREQ(actual_string_value.c_str(), expected_string_value.c_str());
  }
}

//endregion

//region Exceptions:

TEST(config_parser, cttor_reader_throws_ini_parser_error)
{
  auto file_type = config_parser::FileType::INI;
  std::string file_path = "config_not.ini";

  EXPECT_THROW(
      dft_core::config_parser::ConfigParser(file_path, file_type),
      boost::property_tree::ini_parser::ini_parser_error
  );
}

TEST(config_parser, cttor_reader_throws_json_parser_error)
{
  auto file_type = config_parser::FileType::JSON;
  std::string file_path = "config_not.json";

  EXPECT_THROW(
      dft_core::config_parser::ConfigParser(file_path, file_type),
      boost::property_tree::json_parser::json_parser_error
  );
}

TEST(config_parser, cttor_reader_throws_xml_parser_error)
{
  auto file_type = config_parser::FileType::XML;
  std::string file_path = "config_not.xml";

  EXPECT_THROW(
      dft_core::config_parser::ConfigParser(file_path, file_type),
      boost::property_tree::xml_parser::xml_parser_error
  );
}

TEST(config_parser, cttor_reader_throws_info_parser_error)
{
  auto file_type = config_parser::FileType::INFO;
  std::string file_path = "config_not.info";

  EXPECT_THROW(
      dft_core::config_parser::ConfigParser(file_path, file_type),
      boost::property_tree::info_parser::info_parser_error
  );
}

//endregion



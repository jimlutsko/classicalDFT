#include "dft_lib/utils/config_parser.h"

#include <utility>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "dft_lib/utils/console.h"

namespace dft_core
{
namespace config_parser
{
ConfigParser::ConfigParser()
{
  this->ReadConfigFile(this->config_file_path(), this->config_file_type());
}

ConfigParser::ConfigParser(std::string config_file, const FileType& file_type)
    : config_file_path_(std::move(config_file)), config_file_type_(file_type)
{
  this->ReadConfigFile(this->config_file_path(), this->config_file_type());
}

void ConfigParser::SetConfigFilePath(const std::string &config_file)
{
  config_file_path_ = config_file;
}

void ConfigParser::SetConfigFileType(const FileType &file_type)
{
  config_file_type_ = file_type;
}

void ConfigParser::ReadConfigFile(const std::string& config_file, const FileType& file_type)
{
  this->SetConfigFilePath(config_file);
  this->SetConfigFileType(file_type);

  switch(file_type) {
    case FileType::INI:
      boost::property_tree::ini_parser::read_ini(this->config_file_path(), this->config_tree_);
      break;
    case FileType::JSON:
      boost::property_tree::json_parser::read_json(this->config_file_path(), this->config_tree_);
      break;
    case FileType::INFO:
      boost::property_tree::info_parser::read_info(this->config_file_path(), this->config_tree_);
      break;
    case FileType::XML:
      boost::property_tree::xml_parser::read_xml(this->config_file_path(), this->config_tree_);
      break;
  }

}
}
}
#ifndef CLASSICALDFT_CONFIG_PARSER_H
#define CLASSICALDFT_CONFIG_PARSER_H

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser/error.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>

namespace dft_core
{
namespace config_parser
{

/**
 * The various configuration formats supported
 */
enum class FileType {
  INI = 0,
  JSON,
  INFO,
  XML
};

/// The default name of the configuration file:
const std::string DEFAULT_CONFIG_FILE_NAME = "config.ini";
/// The default type of the configuration file:
const FileType DEFAULT_CONFIG_FILE_TYPE = FileType::INI;

/**
 * @brief A configuration-file parser based on boost::property_tree::ptree
 *
 * @detailed The ConfigParser class is a wrapper of the boost::property_tree::ptree class
 * whilst utilising `boost::property_tree::read_ini` and `boost::property_tree::write_ini` methods.
 * This allows us to ingest/export configuration values from/to an external file.
 * The standard file templates compatible with this wrapper class are the industry standards, i.e.
 * `INI`, `JSON`, `XML` and `INFO`.
 */
class ConfigParser
{
 private:
  /// The configuration file full path (default = `config.ini`)
  std::string config_file_path_ = DEFAULT_CONFIG_FILE_NAME;
  /// The configuration file type (default = `config_parser::INI`)
  FileType config_file_type_ = DEFAULT_CONFIG_FILE_TYPE;
  /// The actual 'boost::property_tree::ptree' object which has all the parsing functionality
  boost::property_tree::ptree config_tree_;

 public:
  //region Cttors:

  /**
   * @brief Wrapper of the boost::property_tree functionality. The default constructor will assume
   * we are consuming a `config.ini` file
   *
   * @throw boost::property_tree::ini_parser_error in case of error deserialising the property tree
   */
  ConfigParser();

  /**
   * @brief Wrapper of the boost::property_tree functionality
   *
   * @param config_file the path to the configuration file to parse
   * @param file_type must be one of the supported files, specified in the enum `FileType`
   * @throw boost::property_tree::ini_parser_error In case of error deserialising the property tree (if FileType::INI)
   * @throw boost::property_tree::json_parser_error In case of error deserialising the property tree (if FileType::JSON)
   * @throw boost::property_tree::xml_parser_error In case of error deserialising the property tree (if FileType::XML)
   * @throw boost::property_tree::info_parser_error In case of error deserialising the property tree (if FileType::INFO)
   */
  explicit ConfigParser(std::string config_file, const config_parser::FileType& file_type = DEFAULT_CONFIG_FILE_TYPE);

  //endregion

  //region Mutators:

  /**
   * @brief Setter of the `config_file_path_` private attribute
   *
   * @param[in] config_file std::string with the full path to find the configuration file
   */
  void SetConfigFilePath(const std::string& config_file);

  /**
   * #brief Setter of the  `config_file_type_` private attribute
   *
   * @param[in] file_type is one of the possible `config_parser::FileType` values
   */
  void SetConfigFileType(const FileType& file_type);

  /**
   * @brief Wrapper which decides the correct parser to use from `boost::property_tree` based
   * on the `file_type`
   *
   * @param[in] config_file the full path to find the configuration file
   * @param[in] file_type one of the possible `config_parser::FileType` values
   */
  void ReadConfigFile(const std::string& config_file, const FileType& file_type);

  //endregion

  //region Inspectors:

  /// Gets the config_file_path_ private attribute
  const std::string& config_file_path() const { return config_file_path_; }
  /// Gets the config_file_type_ private attribute
  const FileType& config_file_type() const { return config_file_type_; }
  /// Gets the config_tree_ private attribute
  boost::property_tree::ptree tree() const { return config_tree_; };

  //endregion
};
}}
#endif  // CLASSICALDFT_CONFIG_PARSER_H

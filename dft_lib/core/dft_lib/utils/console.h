#ifndef CLASSICALDFT_CONSOLE_H
#define CLASSICALDFT_CONSOLE_H

#include <chrono>
#include <iostream>
#include <iomanip>

namespace console
{
  /**
   * @brief Colors to be used when writing in the terminal
   */
  namespace color
  {
    const std::string red {"\033[0;31m"};
    const std::string green {"\033[1;32m"};
    const std::string yellow {"\033[1;33m"};
    const std::string cyan {"\033[0;36m"};
    const std::string magenta {"\033[0;35m"};
    const std::string reset {"\033[0m"};
  }

  /**
   * @brief Methods to format the text when writing in the terminal
   */
  namespace format
  {
    /**
     * @brief Returns the input message string as a bold string
     * @param[in] msg the message string to be formatted
     * @return std::string
     */
    static std::string bold(const std::string& msg)
    {
      return "\x1b[1m" + msg + "\x1b[0m";
    }

    /**
     * @brief Returns the input message string as a blinking string
     * @param[in] msg the message string to be formatted
     * @return std::string
     */
    static std::string blink(const std::string& msg)
    {
      return "\033[33;5;7m" + msg + "\033[0m";
    }
  }

  /**
   * @brief Wrapper of the std::cout method
   * @param[in] string message to write in the terminal
   */
  static void write(const std::string& msg)
  {
    std::cout << msg;
  }

  /**
   * @brief Wrapper of the std::cout + std::endl method
   * @param[in] string message to write in the terminal
   */
  static void write_line(const std::string& msg)
  {
    std::cout << msg << std::endl;
  }

  /**
   * @brief Wrapper of std::endl
   */
  static void new_line()
  {
    std::cout << std::endl;
  }

  /**
   * @brief Pauses the terminal until any character is introduced
   */
  static void pause()
  {
    std::cin.ignore();
  }

  /**
   * @brief Wrapper of the std::in method which returns whatever written in terminal as a std::string
   * @returns std::string the written characters in terminal
   */
  static std::string read_line()
  {
    std::string out;
    std::cin >> out;
    return out;
  }

  /**
   * @brief Pauses the terminal until any character is introduced
   */
  static void wait()
  {
    console::write_line("Press enter to continue...");
    std::cin.ignore();
  }

  /**
   * @brief Returns the current time as std::string in format YYYY-MM-DD H:M:S
   */
  static std::string _now_str()
  {
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(19, '\0');
    std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return s;
  }

  /**
   * @brief Writes in terminal a given message in green color specifying the time and with
   * a decorator text "[i] Info" in front of the message
   * @param[in] msg the message to write as "information" message
   */
  static void info(const std::string& msg)
  {
    std::string message = color::green + _now_str() + " | " + "[i] Info: "  + msg + color::reset;
    console::write_line(message);
  }

  /**
   * @brief Writes in terminal a given message in yellow color specifying the time and with
   * a decorator text "[?] Warning" in front of the message
   * @param[in] msg the message to write as "warning" message
   */
  static void warning(const std::string& msg)
  {
    std::string message = color::yellow + _now_str() + " | " + "[?] Warning: " + msg + color::reset;
    console::write_line(message);
  }

  /**
   * @brief Writes in terminal a given message in red color specifying the time and with
   * a decorator text "[!] Error" in front of the message
   * @param[in] msg the message to write as "error" message
   */
  static void error(const std::string& msg)
  {
    std::string message = color::red + _now_str() + " | " +  "[!] Error: " + msg + color::reset;
    console::write_line(message);
  }

  /**
   * @brief Writes in terminal a given message in cyan color specifying the time and with
   * a decorator text "[+] Debug" in front of the message
   * @param[in] msg the message to write as "warning" message
   */
  static void debug(const std::string& msg)
  {
    std::string message = color::cyan + _now_str() + " | " + "[+] Debug: " + msg + color::reset;
    console::write_line(message);
  }
}

#endif //CLASSICALDFT_CONSOLE_H
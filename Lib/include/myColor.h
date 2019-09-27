/************************************************************************
 * Copyright 2019 James F. Lutsko
 * jim@lutsko.com
 * http://www.lutsko.com
 *
 *  This file is part of AMECRYS_kMC.
 *
 *  AMECRYS_kMC is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AMECRYS_kMC is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AMECRYS_kMC.  If not, see <https://www.gnu.org/licenses/>.
 **************************************************************************/




#ifndef __LUTSKO__MY_COLOR__
#define __LUTSKO__MY_COLOR__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <cstring>

extern "C" {
#include <unistd.h>
}

enum class style : unsigned char {
  Reset     = 0,
    bold      = 1,
    dim       = 2,
    italic    = 3,
    underline = 4,
    blink     = 5,
    reversed  = 6,
    conceal   = 7,
    crossed   = 8
    };
enum class fg : unsigned char {
  def     = 39,
    black   = 30,
    red     = 31,
    green   = 32,
    yellow  = 33,
    blue    = 34,
    magenta = 35,
    cyan    = 36,
    gray    = 37
    };
enum class bg : unsigned char {
  def     = 49,
    black   = 40,
    red     = 41,
    green   = 42,
    yellow  = 43,
    blue    = 44,
    magenta = 45,
    cyan    = 46,
    gray    = 47
    };

class myColor
{
 public:
 myColor(fg code) : code_(static_cast<int>(code)){isAllowed_ = isTerminal() && supportsColor() ? true : false;}
 myColor(style code) : code_(static_cast<int>(code)){isAllowed_ = isTerminal() && supportsColor() ? true : false;}
 myColor(bg code) : code_(static_cast<int>(code)){isAllowed_ = isTerminal() && supportsColor() ? true : false;}

  bool isAllowed() const { return isAllowed_;}
  int code() const { return code_;}

  friend std::ostream &operator<<(std::ostream& os, const myColor& c)
  {
    os.flush(); 
    if(c.isAllowed()) std::cout  << "\e[" << static_cast<int>(c.code()) << "m";
    return os;
  }

  const static myColor BLACK;
  const static myColor RED;
  const static myColor GREEN;
  const static myColor YELLOW;
  const static myColor BLUE;
  const static myColor MAGENTA;
  const static myColor CYAN;
  const static myColor GRAY;
  const static myColor DEF;    

  const static myColor RESET;
  const static myColor BOLD;    
  
 private:
  int code_ = 0;//fg::def;
  bool isAllowed_ = false;

  bool isTerminal()
  {
    return isatty(STDOUT_FILENO);
  }

  bool supportsColor()
  {
    if(const char *env_p = std::getenv("TERM")) {
      const char *const term[8] = {
	"xterm", "xterm-256", "xterm-256color", "vt100",
	"color", "ansi",      "cygwin",         "linux"};
      for(unsigned int i = 0; i < 8; ++i) {
	if(std::strcmp(env_p, term[i]) == 0) return true;
      }
    }
    return false;
  }  
};

#endif

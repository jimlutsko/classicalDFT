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




#ifndef __LUTSKO__LOG__
#define __LUTSKO__LOG__

#include <chrono>
#include <ctime>
#include <sstream>

#include "options.h"
#include "config.h"
#include "GitSHA1.h"


// Write a stream buffer that prefixes each line with Plop
class LogStreamBuf: public std::stringbuf
{
 public:
  LogStreamBuf(ofstream& log) : log_(log) {}
 ~LogStreamBuf(){ if(pbase() != pptr()) putOutput();}

  // When we sync the stream with the output. 
  // 1) Output Plop then the buffer
  // 2) Reset the buffer
  // 3) flush the actual output stream we are using.
  virtual int sync() { putOutput(); return 0;}

  void putOutput()
  {
    //    cout << std::boolalpha;
    // Called by destructor.
    // destructor can not call virtual methods.
    const string &s = str();
    cout << s;

    if(atNewLine_)
      {
	std::wstring::size_type pos = s.find_first_not_of(" \t");
	if(pos != string::npos)
	  {
	    char c = s[pos];
	    if(c != '-' && !isdigit(c))
	      {
		log_ << "#";
	      }
	  }
      }
    log_ << s;
    if(s.back() == '\n') atNewLine_ = true;
    else atNewLine_ = false;
    
    str("");
    cout.flush();
    log_.flush();
  }

 private:
  ofstream &log_;
  bool atNewLine_ = true;
};


/**
  *  @Brief A utility object that sends output to a log file
  */  
class Log: public std::ostream
{
 public:

  /**
   *   @brief  Constructor - registers the name of the logfile and opens the output stream. Also initializes the log file with the current time.
   *
   *   @param name is the name of the log file.
   *   @param Major is the major version number
   *   @param Minor is the minor version number
   *   @param prog is the name of the program (null means the name and version number are not output)
   *   @param numtasks is the number of MPI tasks (-1 means no mpi)
   *   @param isNew : will create new log file if true, will try to append to existing file if false
   */  
 Log(string &name, int Major = -1, int Minor = -1, char *prog = NULL, int numtasks = -1, bool isNew = true) : Log(name.c_str(), Major, Minor,prog, numtasks, isNew) {}

  /**
   *   @brief  Default Constructor: default name of log file is "log.dat".  Also initializes the log file with the current time.
   *
   */    
 Log() : Log("log.dat"){}

  /**
   *   @brief  Constructor - registers the name of the logfile and opens the output stream.  Also initializes the log file with the current time.
   *
   *   @param name is the name of the log file.
   *   @param Major is the major version number
   *   @param Minor is the minor version number
   *   @param prog is the name of the program (null means the name and version number are not output)
   *   @param numtasks is the number of MPI tasks (-1 means no mpi)
   *   @param isNew : will create new log file if true, will try to append to existing file if false
   */    
  Log(const char *name, int Major = -1, int Minor = -1, char *prog = NULL, int numtasks = -1, bool isNew = true, bool verbose = true)
   : log_(name,(isNew ? ios::trunc : ios::app)), buffer(log_), ostream(&buffer)    
    {
      auto now = std::chrono::system_clock::now();
      std::time_t now_time = std::chrono::system_clock::to_time_t(now);

      if(verbose) *this << "*****************************************************************" << endl;
      if(verbose) if(prog != NULL) *this << prog << " version " << Major << "." << Minor << endl;
      if(verbose) *this << std::ctime(&now_time) << " " << endl;
      if(verbose) *this << "#Using:\tLib " << PROJECT_NAME << endl
	    << "\tversion: " << PROJECT_VER << endl
	    << "\tgit revision: " << g_GIT_SHA1 << endl;
      //      *this << "Library built " << _TIMEZ_  << endl;
      if(verbose) *this << "\tLinear Algebra library: " << DFT_Vec::get_library_name() << endl;
      if(verbose) if(numtasks > 0) *this << " MPI: numtasks = " << numtasks << endl;
      if(verbose) *this << "*****************************************************************" << endl  << endl;      
    }

  /**
   *   @brief  Destructor - flushes and closes the log file. 
   */      
  ~Log() {log_.flush(); log_.close();}

  /**
   *   @brief  Write the options object to the log file (effectively records the input file).
   *  
   *   @param options is the Options object to copy
   *
   *   @returns none
   */
  void write(Options &options)
  {
    *this << "#Input parameters:" << endl <<  "#" << endl;
    options.write(*this);
    *this << "=================================" << endl;
    this->flush();
  }

 private:
  ofstream log_;           ///< The stream in which the log is recorded.
  LogStreamBuf buffer;      ///< The buffer that controls writing
};

#endif // __LUTSKO__LOG__

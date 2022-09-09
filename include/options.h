#ifndef __OPIONS_H__
#define __OPIONS_H__

#include <map>
#include <stdexcept>
#include <vector>

//typedef std::pair<char const *, int *> int_option_pair;
//typedef std::pair<char const *, long *> long_option_pair;
//typedef std::pair<char const *, double *> d_option_pair;


/**
  *  @brief  UTILITY: reads and stores values from an input file.
  *
  *  @detailed A simple class for reading an input file. The basic use is like this:  <br>
  *            int main(char argc, char** argv) <br>
  *            { <br>
  *                 Options options;  <br>
  *                 double x;  <br>  <br>
  *
  *                 options.addOption("Xvar", &x);  <br> <br>
  *              
  *                 options.read(argc, argv);  <br> <br>
  *
  *                 ofstream log("log.dat");  <br>
  *                 options.write(log); <br> <br>
  *
  *                 etc. <br>
  *            } <br> <br>
  *               
  *       This reads the value of x from an input file passed as the first command line argument and having the format
  *
  *       #Comments start with hash <br>
  *       Xvar = 3.1 <br>
  *      
  *       The values of all variables registered with option is then written to the output stream "log".
  */


class Options
{
 public:
  Options() {}
  ~Options(){}

  void addOption(char const * name, int * place)
    {
	string s(name);
	intOptions_[s] = place;
    }
  void addOption(char const * name, long * place)
    {
      longOptions_[name] = place;
    }
  void addOption(char const * name, double * place)
    {
     dOptions_[name] = place;
    }
  void addOption(char const * name, string * place)
    {
     cOptions_[name] = place;
    }
  void addOption(char const * name, bool * place)
    {
     bOptions_[name] = place;
    }
  void addOption(char const * name, vector<double> * place)
    {
     vdOptions_[name] = place;
    }
  void addOption(char const * name, vector<int> * place)
    {
     viOptions_[name] = place;
    }  

  void read(int argc, char ** argv, bool bPrint = true);
  void read(char const * file, bool bPrint  = true);
  void write(ostream &of) const;
  
  int getIntOption(string&  name) const 
    {
      std::map<string, int*>::const_iterator i = intOptions_.find(name);
      if (i == intOptions_.end ())
	throw std::runtime_error("Key not found in options");
      return *(i->second);
    }
  int getLongOption(string&  name) const 
    {
      std::map<string, long*>::const_iterator i = longOptions_.find(name);
      if (i == longOptions_.end ())
	throw std::runtime_error("Key not found in options");
      return *(i->second);
    }
  double getDoubleOption(string& name) const 
    {
      std::map<string, double*>::const_iterator i = dOptions_.find(name);
      if (i == dOptions_.end ())
	throw std::runtime_error("Key not found in options");
      return *(i->second);
    }
  string* getCharOption(string& name) const 
    {
      std::map<string, string*>::const_iterator i = cOptions_.find(name);
      if (i == cOptions_.end ())
	throw std::runtime_error("Key not found in options");
      return (i->second);
    }
  bool getBoolOption(string& name) const 
    {
      std::map<string, bool*>::const_iterator i = bOptions_.find(name);
      if (i == bOptions_.end ())
	throw std::runtime_error("Key not found in options");
      return *(i->second);
    }
  std::vector<double> getVectorDoubleOption(string& name) const 
    {
      std::map<string, vector<double>*>::const_iterator i = vdOptions_.find(name);
      if (i == vdOptions_.end ())
	throw std::runtime_error("Key not found in options");
      return *(i->second);
    }
  std::vector<int> getVectorIntOption(string& name) const 
    {
      std::map<string, vector<int>*>::const_iterator i = viOptions_.find(name);
      if (i == viOptions_.end())
	throw std::runtime_error("Key not found in options");
      return *(i->second);
    }  
    
 private:

  std::map<string, int*>    intOptions_;
  std::map<string, long*>   longOptions_;
  std::map<string, double*> dOptions_;
  std::map<string, string*> cOptions_;
  std::map<string, bool*>   bOptions_;
  std::map<string, vector<double>*>   vdOptions_;
  std::map<string, vector<int>*>      viOptions_;
};

#endif //__OPIONS_H__

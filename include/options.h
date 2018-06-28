#ifndef __OPIONS_H__
#define __OPIONS_H__

#include <map>
#include <stdexcept>

//typedef std::pair<char const *, int *> int_option_pair;
//typedef std::pair<char const *, long *> long_option_pair;
//typedef std::pair<char const *, double *> d_option_pair;

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

  void read(int argc, char ** argv, bool bPrint = true);
  void read(char const * file, bool bPrint  = true);
  void write(ofstream &of) const;

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
  double getBoolOption(string& name) const 
    {
      std::map<string, bool*>::const_iterator i = bOptions_.find(name);
      if (i == bOptions_.end ())
	throw std::runtime_error("Key not found in options");
      return *(i->second);
    }
 private:

  std::map<string, int*>    intOptions_;
  std::map<string, long*>   longOptions_;
  std::map<string, double*> dOptions_;
  std::map<string, string*> cOptions_;
  std::map<string, bool*>   bOptions_;
};

#endif //__OPIONS_H__

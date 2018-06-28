#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>

using namespace std;

#include "options.h"

void Options::read(int argc, char ** argv, bool bPrint)
{
  if(argc < 2)
  {
      stringstream message(stringstream::in | stringstream::out);
      message <<  "Usage: " << argv[0] << "  <input.dat>";
      throw std::runtime_error(message.str().c_str());
  }
  if(bPrint) cout << "input file is : " << argv[1] << endl;
  read(argv[1],bPrint);
}


void Options::read(char const * file, bool bPrint)
{
  if(bPrint) cout << endl << "++++++++++++++++++++ reading parameters from " << file << " +++++++++++++++++++" << endl << endl;

  ifstream f(file,ios::in);
  

  if(!f.good()) 
      throw std::runtime_error("Could not open input file ... ");

  const char * delim = "= ";

  while (f.good())
    {
      char buf[256];
      f.getline(buf,256);

      if(buf[0] == '#') continue;
      
      char * pch = strtok (buf,delim);
      if(!pch) continue;

      if(bPrint) cout << pch << " = ";

      if(intOptions_.find(pch) != intOptions_.end())
	{
	  int *place = intOptions_[pch];
	  pch = strtok(NULL,delim);
	  *place = atoi(pch);
	  if(bPrint) cout << *place << endl;
	} else if(longOptions_.find(pch) != longOptions_.end()) {
	  long *place = longOptions_[pch];
	  pch = strtok(NULL,delim);
	  *place = atol(pch);
	  if(bPrint) cout << *place << endl;
	} else if(dOptions_.find(pch) != dOptions_.end()) {
	  double *place = dOptions_[pch];
	  pch = strtok(NULL,delim);
	  *place = atof(pch);
	  if(bPrint) cout << *place << endl;
	} else if(cOptions_.find(pch) != cOptions_.end()) {
	  string *place = cOptions_[pch];
	  pch = strtok(NULL,delim);
	  if(pch != NULL) place->assign(pch);
	  if(bPrint) cout << *place << endl;
	} else if(bOptions_.find(pch) != bOptions_.end()) {
	  bool *place = bOptions_[pch];
	  pch = strtok(NULL,delim);
	  if(strcmp(pch,"true") == 0) *place = true;
	  else if(strcmp(pch,"false") == 0) *place = false;
	  else {
	      if(bPrint) cout << "Input was : " << pch << endl;
	      throw std::runtime_error("Unrecognized boolean input in util//options.cpp");
	  }
	  if(bPrint) cout << boolalpha << *place << endl;
	} else {
	  if(bPrint) cout <<  "<not a parameter>" << endl;
	}
    }
  f.close();

  if(bPrint) cout << endl << "++++++++++++++++++++ end parameters +++++++++++++++++++" << endl << endl;
}

void Options::write(ofstream &of) const
{
  map<string, int*>::const_iterator iter;
  for (iter=intOptions_.begin(); iter != intOptions_.end(); ++iter) 
    of << "#" <<  iter->first << " = " << *(iter->second) << endl;

  map<string, long*>::const_iterator iterl;
  for (iterl=longOptions_.begin(); iterl != longOptions_.end(); ++iterl) 
    of << "#" <<  iterl->first << " = " << *(iterl->second) << endl;

  map<string, double*>::const_iterator iter1;
  for (iter1=dOptions_.begin(); iter1 != dOptions_.end(); ++iter1) 
    of << "#" <<  iter1->first << " = " << *(iter1->second) << endl;

  map<string, string*>::const_iterator iter2;
  for (iter2=cOptions_.begin(); iter2 != cOptions_.end(); ++iter2) 
    of << "#" <<  iter2->first << " = " << (iter2->second == NULL ? "" : *(iter2->second)) << endl;
    
  map<string, bool*>::const_iterator iterb;
  for (iterb=bOptions_.begin(); iterb != bOptions_.end(); ++iterb) 
    of << "#" <<  iterb->first << " = " << *(iterb->second) << endl;

  of << "#" <<  endl;
}

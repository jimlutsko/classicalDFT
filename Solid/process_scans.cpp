#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <iterator>
#include <dirent.h>

using namespace std;

void removeSubstring(string &s, string &toRemove)
{
  size_t pos = std::string::npos;
  while ((pos  = s.find(toRemove) )!= std::string::npos)
    s.erase(pos, toRemove.length());
}

template <class Container>
void split4(const std::string& str, Container& cont,
              const std::string& delims = " ")
{
    std::size_t current, previous = 0;
    current = str.find_first_of(delims);
    while (current != std::string::npos) {
        cont.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find_first_of(delims, previous);
    }
    cont.push_back(str.substr(previous, current - previous));
}

void quadraticFit(double x1, double x2, double x3, double y1, double y2, double y3, double &a0, double &a1, double &a2)
{
  a2 = (x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))/((x1-x2)*(x1-x3)*(x2-x3));
  a1 = ((y2-y1)/(x2-x1))-a2*(x1+x2);
  a0 = y1-a2*x1*x1-a1*x1;
}

void parseHeader(string &line, vector<double> &header)
{
  string delim(" =");
  std::vector<std::string> words;
  split4(line, words, delim);
  for(auto &x : words)
    {
      try {
	double d = stod(x);
	header.push_back(d);
      } catch(...){}
    }
}

void processFile(string &s, vector<double> &header1, vector<double> &header2, double &Dsol, double &Osol, double &Cvac)
{
  cout << "File " << s << endl;
  string::size_type sz;
  ifstream in(s);
  string line;

  // parse the first line to get the thermodynamic variables  
  getline(in,line); parseHeader(line, header1);
  getline(in,line); parseHeader(line, header2);
  
  // skip the next line if there was a second header
  if(header2.size() != 0)
    getline(in,line);

  // Read values from each line
  vector< vector<double> > vals; 
  while(getline(in,line))
    {
      stringstream ss(line);
      vector<double> v(7);
      ss >> v[0]  >> v[1]  >> v[2]  >> v[3]  >> v[4]  >> v[5]  >> v[6];
      vals.push_back(v);
    }
  // Find minimum energy
  if(vals.size() < 3) throw std::runtime_error("Not enough values");
  int p = 0;
  if(vals[p][4] < vals[p+1][4]) {stringstream ss; ss << "No minimum in file " << s; throw std::runtime_error(ss.str().c_str());}
  for(int i=3;i<vals.size() && vals[p+2][4] < vals[p+1][4];i++)
    p++;
  if(vals[p+2][4] < vals[p+1][4]) throw std::runtime_error("No minimum");
  // do simple quadratic fit
  double a0,a1,a2;
  quadraticFit(vals[p][0],vals[p+1][0], vals[p+2][0],vals[p][4], vals[p+1][4], vals[p+2][4],a0,a1,a2);
  double x = -a1/(2*a2);
  Osol = a0+a1*x+a2*x*x;

  quadraticFit(vals[p][0],vals[p+1][0], vals[p+2][0],vals[p][3], vals[p+1][3], vals[p+2][3],a0,a1,a2);  
  Dsol = a0+a1*x+a2*x*x;

  quadraticFit(vals[p][0],vals[p+1][0], vals[p+2][0],vals[p][5], vals[p+1][5], vals[p+2][5],a0,a1,a2);  
  Cvac = a0+a1*x+a2*x*x;  
}

int main()
{
  string prefix("scan_mu_");
  string suffix(".dat");


  vector< vector<double> > thermo1;
  vector< vector<double> > thermo2;
  DIR *dirp = opendir(".");
  struct dirent *dp;  
  while ((dp = readdir(dirp)) != NULL)
    {
      string s(dp->d_name);

      double Mu, Dliq, Oliq, Dsol, Osol, Cvac;
      
      if(s.rfind(prefix, 0) == 0)
	{
	  vector<double> header1;
	  vector<double> header2;
	  
	  processFile(s, header1, header2, Dsol, Osol, Cvac);

	  vector<double> v;
	  v.push_back(header1[0]);
	  v.push_back(header1[1]);
	  v.push_back(header1[2]);
	  v.push_back(Dsol);
	  v.push_back(Osol);
	  v.push_back(Cvac);
	  thermo1.push_back(v);
	  
	  if(header2.size() > 0)
	    {
	      vector<double> v2(v);
	      v2.push_back(header2[0]);
	      v2.push_back(header2[1]);
	      v2.push_back(header2[2]);
	      thermo2.push_back(v2);	      
	    }
	      
	}
    }
  std::sort(thermo1.begin(), thermo1.end(), [](const std::vector<double>& a, const std::vector<double>& b) {return a[0] < b[0];});
  if(thermo2.size() > 0)
    std::sort(thermo2.begin(), thermo2.end(), [](const std::vector<double>& a, const std::vector<double>& b) {return a[0] < b[0];});

  for(int i = 0; i < 2; i++)
    {
      if(i == 1 && thermo2.size() < 1) continue;

      vector<vector<double>> &thermo = (i == 0 ? thermo1 : thermo2);
  
      int nw = 12;
      int p = -1;
      int sgn = 1;
      for(int i=0;i<thermo.size();i++)
	{
	  cout << setw(nw) << "Mu   = " << setw(nw) << thermo[i][0]
	       << setw(nw) << "Dliq = " << setw(nw) << thermo[i][2]
	       << setw(nw) << "Oliq = " << setw(nw) << thermo[i][1]
	       << setw(nw) << "Dsol = " << setw(nw) << thermo[i][3]
	       << setw(nw) << "Osol = " << setw(nw) << thermo[i][4]
	       << setw(nw) << "Cvac = " << setw(nw) << thermo[i][5]
	       << endl;
	  int sgn2 = (thermo[i][1] - thermo[i][4] > 0 ? 1 : -1);
	  if(i == 0) sgn = sgn2;
	  else if(sgn*sgn2 < 0 && p < 0) p = i-1;
	}
      if(p < 0) throw std::runtime_error("No coexistence found");
      double x0 = p-(thermo[p][1]-thermo[p][4])/(thermo[p+1][1]-thermo[p+1][4]-thermo[p][1]+thermo[p][4]);
      p = x0;

      cout << endl;
      cout << "Coexistence: " << endl;
      cout << " Mu   = " << thermo[p][0]+(thermo[p+1][0]-thermo[p][0])*(x0-p) << endl
	   << " Dliq = " << thermo[p][2]+(thermo[p+1][2]-thermo[p][2])*(x0-p) << endl
	   << " Oliq = " << thermo[p][1]+(thermo[p+1][1]-thermo[p][1])*(x0-p) << endl
	   << " Dsol = " << thermo[p][3]+(thermo[p+1][3]-thermo[p][3])*(x0-p) << endl
	   << " Osol = " << thermo[p][4]+(thermo[p+1][4]-thermo[p][4])*(x0-p) << endl
	   << " Cvac = " << thermo[p][5]+(thermo[p+1][5]-thermo[p][5])*(x0-p) << endl
	   << endl;
    }
  
  
}

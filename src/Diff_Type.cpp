#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#include "Minimizer.h"

DiffType get_diff_type_from_string(string diff_type_string)
{
  for(unsigned n = 0; DiffType_names(n); n++)
    if(diff_type_string == DiffType_names(n))
      return DiffType_values(n);
  throw std::runtime_error("get_diff_type_from_string could not find difference type");
}

string get_diff_type_name(DiffType diff)
{
  for (unsigned n = 0; DiffType_names(n); n++)
    {
      if(diff == DiffType_values(n))
	return DiffType_names(n);
    }
  throw std::runtime_error("get_diff_type_name could not find difference type");
}

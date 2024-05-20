#ifndef __LUTSKO_TABLE__
#define __LUTSKO_TABLE__

#include <fstream>
#include <vector>
#include <iomanip>

/**
  *  @brief  UTILITY: reads a table from an input file. The format is assumed to be one row per line with entries separated by white space.
  */
class Table
{
public:
  Table(std::ifstream &in);

  ~Table(){}

  bool isEmpty() const { return (data_.size() == 0 || data_[0].size() == 0);}

  double val(int i, int j) const { return data_[i][j];}
  double get(int i, int j) const { return val(i,j);}
  
  int nRows() const { return data_.size();}
  int nCols() const { return data_[0].size();}
  double colMax(int col) const;
  double colMin(int col) const;

  std::vector<double> getColumn(int j) const;
  
  void write(std::ostream &out) const;
  void write_gnuplot(std::ostream &out) const;

protected:
  std::vector< std::vector<double> > data_;
};

#endif // __LUTSKO_TABLE__ sentinal

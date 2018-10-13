#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

#include "Table.h"

using namespace std;



Table::Table(ifstream &in)
{
    if(!in.good())
	throw std::runtime_error("input file cannot be opened");


    for(string buf; getline(in,buf);)
    {
	vector<double> line;

	if(!isdigit(buf[0]) && buf[0] != '-')
	    continue;

	stringstream st(buf);
	double d;
	while(st >> d) {line.push_back(d);}
	data_.push_back(line);
    }
}


double Table::colMax(int col) const
{
    if(isEmpty())
	throw std::runtime_error("Table is empty");

    if(col < 0 || col >= nCols())
	throw std::runtime_error("col out of range in Table::colMax");

    double tmax = val(0,col);

    for(int i=1;i<nRows();i++)
	tmax = max(tmax, val(i,col));

    return tmax;
}

double Table::colMin(int col) const
{
    if(isEmpty())
	throw std::runtime_error("Table is empty");

    if(col < 0 || col >= nCols())
	throw std::runtime_error("col out of range in Table::colMax");

    double tmin = val(0,col);

    for(int i=1;i<nRows();i++)
	tmin = min(tmin, val(i,col));

    return tmin;
}


void Table::write(ostream &out) const
{
    if(isEmpty())
	throw std::runtime_error("Table is empty");

    for(const vector<double> &line: data_)
    {
	for (const double  &d: line)
	    out << d << "\t";
	out << endl;
    }
}



void Table::write_gnuplot(ostream &out) const
{
    if(isEmpty())
	throw std::runtime_error("Table is empty");

    if(nCols() < 3)
	throw std::runtime_error("Table::write_gnuplot requires at least  three columns");

    if(nRows() < 2)
	throw std::runtime_error("Table::write_gnuplot requires at least 2 rows");

    int comp = 0;

    if(fabs(data_[0][0] -data_[1][0]) > fabs(data_[0][1] - data_[1][1])) comp = 1;

    double xold = data_.front()[comp];

    for(const vector<double> &line: data_)
    {
	if(fabs(line[comp]-xold) > 1e-8)
	{
	    xold = line[comp];
	    out << endl;
	}

	out << line[comp] << "\t" 
	    << line[1-comp] << "\t" 
	    << line[2] << endl;
    }
}

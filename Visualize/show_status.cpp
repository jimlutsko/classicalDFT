
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <unistd.h>

using namespace std;

extern char   __BUILD_DATE;
extern char   __BUILD_NUMBER;

#include "Grace.h"
#include "options.h"
#include "TimeStamp.h"
#include "Table.h"


int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;

  Grace g(800,600,1);

  g.pause();

  bool first = true;
  
  for(int pos = 0; pos > -1; pos++)
    {
      stringstream ss;
      ss << "status_" << pos << ".dat";

      cout << "opening " << ss.str() << endl;

      ifstream in(ss.str());
      if(!in.good()) {pos = -2; continue;}
      Table t(in);

      g.deleteDataSet(0);

      for(int i=0;i<t.nRows();i++)
	{
	  g.addPoint(i,t.val(i,2));
	  if(first) g.addPoint(i,t.val(i,2),1);
	}
      g.redraw();

      if(t.nRows() < 1) pos = -2;

      usleep(100000);
      first = false;
    }
  g.pause();
  g.close();
  return 1;
}

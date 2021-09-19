#ifndef __GRACE__
#define __GRACE__

#include <unistd.h>
#include <grace_np.h>

#include <vector>

/**
  *  @brief  UTILITY: wrapper for xmgrace graphics program.
  *
  *  @detailed This class provides an interface for the xmgrace graphics program which can be called programatically to provide a cheap solution  for displaying dynamically updated line graphs etc.
  */

class Grace
{
 public:
  explicit Grace(int sizex = 800, int sizey = 600, int Ngraph = 1, bool show = true);
  ~Grace(){}

  void close() { GraceClose();}
  void store() const;
  void store(std::string &file) const;

  void kill() const { sendCommand("Kill G0");}

  void addPoint(double x, double y, int N=0, int G=0) const;

  int  addDataSet(std::vector<double> const &x, std::vector<double> const &y);
  void replaceDataSet(std::vector<double> const &x, std::vector<double> const &y, int N = 0) const;
  void deleteDataSet(int N = 0, int G=0) const;
  void setLegend(const char *s, int n);
  void setColor(int color, int nm, int G=0) const;

  void setXAxisLabel(const char *s, int Graph = -1);
  void setYAxisLabel(const char *s, int Graph = -1);

  void setCharSize(double s, int Graph = -1);

  void setTitle(const char *s);
  void setSubTitle(const char *s);

  void redraw(int autoScale = 1, int Graph = -1) const;

  void setLimits(double xmin, double xmax, double ymin, double ymax, int Graph = -1);

  void setTicks(double dx, double dy, int Graph = -1);
  void setNumberXMinorTicks(int n, int Graph = -1);
  void setNumberYMinorTicks(int n, int Graph = -1);
    
  void setSymbol(int set, int symbol, int Graph = -1);
  void setSymbolFill(int set, int color, int Graph =-1);

  void pause() const {char cc; std::cout << "Enter any character to continue: "; std::cin >> cc;}

  void noLine(int dataSet, int Graph = -1);  
  void Line(int dataSet, int Graph = -1);  
  void circle(int dataSet, int Graph = -1);  
  void square(int dataSet, int Graph = -1);  

  void size(int dataSet, double size, int Graph = -1);

  static int CIRCLE;
  static int SQUARE;


  static int JPG;
  static int PNG;
  static int EPS1;
  static int PS;
  static int PDF;

  void printToFile(std::string &file, int format);

  void autotic(int Graph = -1) const;

 private:

  void sendCommand(const std::string &s) const {if(GraceIsOpen()) GracePrintf(s.c_str());}

  double xmin_, xmax_;
  double ymin_, ymax_;
  int nMaxDataSet_;
  int Ngraph_;
};


#endif // __GRACE__ sentinal

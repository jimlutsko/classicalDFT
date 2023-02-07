#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include "Grace.h"


int Grace::CIRCLE = 1;
int Grace::SQUARE = 2;



int Grace::JPG  = 1;
int Grace::PNG  = 2;
int Grace::EPS1 = 3;
int Grace::PS   = 4;
int Grace::PDF  = 5;

int Grace::NORMAL_SCALE = 1;
int Grace::LOG_SCALE    = 2;

void Grace::printToFile(std::string &filename, int format)
{
    // set output format (DEVICE)
    std::stringstream ss1;
    //ss1 << "print to ";
    ss1  << "hardcopy device ";

    if(format == Grace::JPG)
	ss1 << "\"JPG\"";
    else if(format == Grace::PNG)
	ss1 << "\"PNG\"";
    else if(format == Grace::PS)
	ss1 << "\"PS\"";
    else if(format == Grace::EPS1)
	ss1 << "\"EPS1\"";
    else if(format == Grace::PDF)
	ss1 << "\"PDF\"";
    else throw std::runtime_error("Unknown format in Grace::printToFile()");

    sendCommand(ss1.str());

    std::stringstream ss2;
    ss2 << "print to \"" << filename << "\"";
    sendCommand(ss2.str());

    std::stringstream ss3;
    ss3 << "print";
    //ss3 << "hardcopy";
    sendCommand(ss3.str());
}




static void my_error_function(const char *msg)
{
  std::cout <<  "Grace library message: \"" << msg << "\"" << std::endl;
};



Grace::Grace(int xsize, int ysize, int Ngraph, bool show) : xmin_(0.0), xmax_(10.0), ymin_(0.0), 
			   ymax_(10.0), nMaxDataSet_(0), Ngraph_(Ngraph)
{
  if(!show) return;

  GraceRegisterErrorFunction(my_error_function);

  char st[32];
  sprintf(st,"xmgrace");

  std::string s0 = "-free";
  std::string s1 = "-nosafe";
  std::string s2 = "-geometry";
  std::stringstream arg2;
  arg2 << xsize << "x" << ysize;
  
  /* Start Grace with a buffer size of 2048 and open the pipe */
  if (GraceOpenVA(st,2048, s0.c_str(), s1.c_str(), s2.c_str(), arg2.str().c_str(),NULL) == -1)
    throw std::runtime_error("Can't run Grace. \n");
      
  /* Send some initialization commands to Grace */
      
  std::stringstream s;
  if(Ngraph > 1)
    {
      int rows = (Ngraph <= 2 ? 1 : 2);
      int cols = Ngraph/rows;
      cols += (Ngraph-rows*cols);

      s << "ARRANGE(" << rows << ", " << cols << ", 0.1, 0.15, 0.2)";
      sendCommand(s.str());
      s.str(std::string());
    }

  s << "world xmax " << xmax_;
  sendCommand(s.str());
  s.str(std::string());

  s << "world xmin " << xmin_;
  sendCommand(s.str());
  s.str(std::string());

  s << "world ymax "  << ymax_;
  sendCommand(s.str());
  s.str(std::string());
  
  s << "xaxis tick major 5";
  sendCommand(s.str());
  s.str(std::string());
  
  s << "xaxis tick minor 1";
  sendCommand(s.str());
  s.str(std::string());
  
  s << "yaxis tick major " << int(ymax_);
  sendCommand(s.str());
  s.str(std::string());

  s << "yaxis tick minor " << int(ymax_/2);
  sendCommand(s.str());
  s.str(std::string());
      
  s << "AUTOSCALE ONREAD XYAXES";
  sendCommand(s.str());
  s.str(std::string());
}

void Grace::setLimits(double xmin, double xmax, double ymin, double ymax, int Graph)
{
  std::stringstream s;

  if(Graph >= 0) 
    {
      s << "FOCUS G" << Graph;
      sendCommand(s.str());
      s.str(std::string());
    }


  if(xmin < xmax)
    {
      s << "world xmin " << xmin;
      sendCommand(s.str());
      s.str(std::string());
      
      s << "world xmax " << xmax;
      sendCommand(s.str());
      s.str(std::string());
    }
  
  if(ymin < ymax)
    {
      s << "world ymin " << ymin;
      sendCommand(s.str());
      s.str(std::string());
      
      s << "world ymax " << ymax;
      sendCommand(s.str());
      s.str(std::string());
    }
  redraw(0);
}

void Grace::setTicks(double dx, double dy, int Graph)
{
  std::stringstream s;

  if(Graph >= 0) 
    {
      s << "FOCUS G" << Graph;
      sendCommand(s.str());
      s.str(std::string());
    }

  if(dx > 0)
    {
      s << "xaxis tick major " << dx;
      sendCommand(s.str());
      s.str(std::string());
    }

  if(dy > 0)
    {
      s << "yaxis tick major " << dy;
      sendCommand(s.str());
      s.str(std::string());
    }
}

void Grace::setNumberXMinorTicks(int n, int Graph)
{
  std::stringstream s;

  if(Graph >= 0) 
    {
      s << "FOCUS G" << Graph;
      sendCommand(s.str());
      s.str(std::string());
    }

  s << "xaxis  tick minor ticks " << n;
  sendCommand(s.str());
  s.str(std::string());
}


void Grace::setNumberYMinorTicks(int n, int Graph)
{
  std::stringstream s;

  if(Graph >= 0) 
    {
      s << "FOCUS G" << Graph;
      sendCommand(s.str());
      s.str(std::string());
    }

  s << "yaxis  tick minor ticks " << n;
  sendCommand(s.str());
}






void Grace::setSymbol(int set, int symbol, int Graph)
{
  std::stringstream ss;

  if(Graph >= 0) 
    {
      ss << "FOCUS G" << Graph;
      sendCommand(ss.str());
      ss.str(std::string());
    }
  

  ss << "s" << set << " symbol " << symbol;
  sendCommand(ss.str());
}

void Grace::setLegend(const char *s, int set)
{
  std::stringstream ss;
  /*
  if(Graph >= 0) 
    {
      ss << "FOCUS G" << Graph;
      sendCommand(ss.str());
      ss.str(std::string());
    }
  */
  ss << "s" << set << " legend " << "\"" << s << "\"";
  sendCommand(ss.str());
}

void Grace::setXAxisLabel(const char *s, int Graph)
{ 
  std::stringstream ss;

  if(Graph >= 0) 
    {
      ss << "FOCUS G" << Graph;
      sendCommand(ss.str());
      ss.str(std::string());
    }

  ss << "XAXIS LABEL \"" << s << "\"";
  sendCommand(ss.str());
}

void Grace::setXAxisScale(const int scale, int Graph)
{ 
  std::stringstream ss;

  if(Graph >= 0) 
    {
      ss << "FOCUS G" << Graph;
      sendCommand(ss.str());
      ss.str(std::string());
    }

  ss << "XAXES SCALE ";
  if(scale == Grace::NORMAL_SCALE) ss << "Normal";
  else if(scale == Grace::LOG_SCALE) ss << "Logarithmic";
  else throw std::runtime_error("Unknown scale request in Grace::setXAxisScale");

  sendCommand(ss.str());
}

void Grace::setYAxisScale(const int scale, int Graph)
{ 
  std::stringstream ss;

  if(Graph >= 0) 
    {
      ss << "FOCUS G" << Graph;
      sendCommand(ss.str());
      ss.str(std::string());
    }

  ss << "YAXES SCALE "; 
  if(scale == Grace::NORMAL_SCALE) ss << "Normal";
  else if(scale == Grace::LOG_SCALE) ss << "Logarithmic";
  else throw std::runtime_error("Unknown scale request in Grace::setYAxisScale");

  sendCommand(ss.str());
}



void Grace::setCharSize(double s, int Graph)
{ 
  std::stringstream ss;

  if(Graph >= 0) 
    {
      ss << "FOCUS G" << Graph;
      sendCommand(ss.str());
      ss.str(std::string());
    }

  ss << "XAXIS LABEL CHAR SIZE " << s;
  sendCommand(ss.str());

  ss.str(std::string());
  ss << "YAXIS LABEL CHAR SIZE " << s;
  sendCommand(ss.str());

  ss.str(std::string());
  ss << "XAXIS TICKLABEL CHAR SIZE " << s;
  sendCommand(ss.str());

  ss.str(std::string());
  ss << "YAXIS TICKLABEL CHAR SIZE " << s;
  sendCommand(ss.str());      
}


void Grace::setYAxisLabel(const char *s, int Graph)
{ 
  std::stringstream ss;

  if(Graph >= 0) ss << "FOCUS G" << Graph;
  sendCommand(ss.str());
  ss.str(std::string());

  ss << "YAXIS LABEL \"" << s << "\"";
  sendCommand(ss.str());
}

void Grace::store() const 
{
  std::string s = "saveall \"sample.agr\"";
  int r = system("rm sample.agr");
  sendCommand(s);
}

void Grace::store(std::string &file) const
{
    std::stringstream s;
    s << "saveall \"" << file << "\"";
    sendCommand(s.str());
}


void Grace::deleteDataSet(int N, int G) const
{
  if(G < 0) G = 0;

  std::stringstream s;
  s << "kill g" << G << ".s" << N;
  sendCommand(s.str());
  setColor((N%10)+1,N);
}

void Grace::autotic(int Graph) const
{
  std::stringstream s;

  if(Graph >= 0) s << "FOCUS G" << Graph;
  sendCommand(s.str());
  s.str(std::string());

  sendCommand("AUTOTICKS");
}

void Grace::redraw(int autoScale, int Graph) const 
{  
  std::stringstream s;

  if(Graph >= 0) s << "FOCUS G" << Graph;
  sendCommand(s.str());
  s.str(std::string());



  if(autoScale)
    {
      s << "AUTOSCALE"; 
      sendCommand(s.str());
      s.str(std::string());

      s << "AUTOTICKS";
      sendCommand(s.str());
      s.str(std::string());

    }
  s << "redraw";
  sendCommand(s.str());
  s.str(std::string());
}


int Grace::addDataSet(std::vector<double> const &x, std::vector<double> const &y)
{
  int n = nMaxDataSet_++;

  for(unsigned int i=0;i<x.size();i++)
    addPoint(x[i],y[i],n);

  return n;
}

void Grace::replaceDataSet(std::vector<double> const &x, std::vector<double> const &y, int N) const
{
  deleteDataSet(N);
  for(unsigned int i=0;i<x.size();i++)
    addPoint(x[i],y[i],N);
}

void Grace::setColor(int color, int set, int G) const
{
  if(G < 0) G = 0;

  std::stringstream s;
  s << "g" << G << ".s" << set << " line color " << color;
  sendCommand(s.str());
  std::stringstream s1;
  s1 << "g" << G << ".s" << set << " symbol color " << color;
  sendCommand(s1.str());
}

void Grace::setSymbolFill(int set, int color, int G)
{
  if(G < 0) G = 0;

  std::stringstream s;
  s << "g" << G << ".s" << set << " symbol fill " << 1;
  sendCommand(s.str());

  std::stringstream s1;
  s1 << "g" << G << ".s" << set << " symbol fill color " << color;
  sendCommand(s1.str());
}

void Grace::addPoint(double x, double y, int N, int G) const
{
  if(G < 0) G = 0;
  std::stringstream s;
  s << "g" << G << ".s" << N << " point " << x << "," <<  y;
  sendCommand(s.str());
}

void Grace::setTitle(const char *s)
{
  std::stringstream ss;
  ss << "Title \"" << s << "\"";
  sendCommand(ss.str());

}

void Grace::setSubTitle(const char *s)
{
  std::stringstream ss;
  ss << "Subtitle \"" << s << "\"";
  sendCommand(ss.str());

}

void Grace::noLine(int dataSet, int Graph)
{
  if(Graph < 0) Graph = 0;
    std::stringstream ss;
    //    ss << "s" << dataSet << " line type 0";
    ss << "g" << Graph << ".s" << dataSet << " line type 0";
    sendCommand(ss.str());
}

void Grace::Line(int dataSet, int Graph)
{
  if(Graph < 0) Graph = 0;
    std::stringstream ss;
   ss << "g" << Graph << ".s" << dataSet << " line type 1";
    sendCommand(ss.str());
}

void Grace::circle(int dataSet, int Graph)
{
  if(Graph < 0) Graph = 0;
    std::stringstream ss;
    ss << "g" << Graph << ".s" << dataSet << " symbol 1";
    sendCommand(ss.str());
}

void Grace::square(int dataSet, int Graph)
{
  if(Graph < 0) Graph = 0;
    std::stringstream ss;
    ss << "g" << Graph << ".s" << dataSet << " symbol 2";
    sendCommand(ss.str());
}

void Grace::size(int dataSet, double size, int Graph)
{
  if(Graph < 0) Graph = 0;
    std::stringstream ss;
    if(Graph >= 0) 
    {
      ss << "FOCUS G" << Graph;
      sendCommand(ss.str());
      ss.str(std::string());
    }

    ss << "g" << Graph << ".s" << dataSet << " symbol size " << size;
    sendCommand(ss.str());
}

void Grace::format_phys_rev(bool overrite_sets, int Graph)
{
  std::stringstream ss;
  
  if(Graph < 0) Graph = 0;
  else 
  {
    ss << "FOCUS G" << Graph;
    sendCommand(ss.str());
    ss.str(std::string());
  }
  
  /*
  ss << "page size 794, 595";
  sendCommand(ss.str());
  ss.str(std::string());
  */
  
  ss << "map color 1 to (0, 0, 0),   \"myblack\"";  sendCommand(ss.str()); ss.str(std::string());
  ss << "map color 2 to (255, 20, 20), \"myred\"";    sendCommand(ss.str()); ss.str(std::string());
  ss << "map color 3 to (75, 75, 255), \"myblue\"";   sendCommand(ss.str()); ss.str(std::string());
  ss << "map color 4 to (0, 175, 0), \"mygreen\"";  sendCommand(ss.str()); ss.str(std::string());
  
  ss << "default linewidth 2.0";   sendCommand(ss.str()); ss.str(std::string());
  ss << "default linestyle 1";     sendCommand(ss.str()); ss.str(std::string());
  ss << "default color 1";         sendCommand(ss.str()); ss.str(std::string());
  ss << "default pattern 1";       sendCommand(ss.str()); ss.str(std::string());
  ss << "default font 0";          sendCommand(ss.str()); ss.str(std::string());
  ss << "default char size 1.5";   sendCommand(ss.str()); ss.str(std::string());
  ss << "default symbol size 1.0"; sendCommand(ss.str()); ss.str(std::string());
  
  ss << "title size 1.5";       sendCommand(ss.str()); ss.str(std::string());
  ss << "subtitle size 1.5";    sendCommand(ss.str()); ss.str(std::string());
  ss << "legend char size 1.5"; sendCommand(ss.str()); ss.str(std::string());
  
  ss << "xaxis label     char size 1.5"; sendCommand(ss.str()); ss.str(std::string());
  ss << "xaxis ticklabel char size 1.5"; sendCommand(ss.str()); ss.str(std::string());
  ss << "xaxis tick minor ticks 1";      sendCommand(ss.str()); ss.str(std::string());
  
  ss << "yaxis label     char size 1.5"; sendCommand(ss.str()); ss.str(std::string());
  ss << "yaxis ticklabel char size 1.5"; sendCommand(ss.str()); ss.str(std::string());
  ss << "yaxis tick minor ticks 1";      sendCommand(ss.str()); ss.str(std::string());
  
  if (overrite_sets)
  {
    //ss << "s0 line type 1";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 line linewidth 2.0";     sendCommand(ss.str()); ss.str(std::string());
    //ss << "s0 line linestyle 1";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 line color 1";           sendCommand(ss.str()); ss.str(std::string());
    //ss << "s0 fill type 2";            sendCommand(ss.str()); ss.str(std::string());
    //ss << "s0 fill rule 0";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 fill color 1";           sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 fill pattern 4";         sendCommand(ss.str()); ss.str(std::string());
    //ss << "s0 symbol 1";               sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 symbol size 1.0";        sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 symbol linewidth 2.0";   sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 symbol linestyle 1";     sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 symbol color 1";         sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 symbol pattern 1";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 symbol fill color 1";    sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 symbol fill pattern 1";  sendCommand(ss.str()); ss.str(std::string());
    ss << "s0 errorbar linewidth 2.0"; sendCommand(ss.str()); ss.str(std::string());
    
    //ss << "s1 line type 1";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 line linewidth 2.0";     sendCommand(ss.str()); ss.str(std::string());
    //ss << "s1 line linestyle 3";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 line color 2";           sendCommand(ss.str()); ss.str(std::string());
    //ss << "s1 fill type 2";            sendCommand(ss.str()); ss.str(std::string());
    //ss << "s1 fill rule 0";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 fill color 2";           sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 fill pattern 5";         sendCommand(ss.str()); ss.str(std::string());
    //ss << "s1 symbol 2";               sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 symbol size 1.0";        sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 symbol linewidth 2.0";   sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 symbol linestyle 1";     sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 symbol color 2";         sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 symbol pattern 1";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 symbol fill color 2";    sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 symbol fill pattern 1";  sendCommand(ss.str()); ss.str(std::string());
    ss << "s1 errorbar linewidth 2.0"; sendCommand(ss.str()); ss.str(std::string());
    
    //ss << "s2 line type 1";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 line linewidth 2.0";     sendCommand(ss.str()); ss.str(std::string());
    //ss << "s2 line linestyle 2";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 line color 3";           sendCommand(ss.str()); ss.str(std::string());
    //ss << "s2 fill type 2";            sendCommand(ss.str()); ss.str(std::string());
    //ss << "s2 fill rule 0";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 fill color 3";           sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 fill pattern 6";         sendCommand(ss.str()); ss.str(std::string());
    //ss << "s2 symbol 3";               sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 symbol size 1.0";        sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 symbol linewidth 2.0";   sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 symbol linestyle 1";     sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 symbol color 3";         sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 symbol pattern 1";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 symbol fill color 3";    sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 symbol fill pattern 1";  sendCommand(ss.str()); ss.str(std::string());
    ss << "s2 errorbar linewidth 2.0"; sendCommand(ss.str()); ss.str(std::string());
    
    //ss << "s3 line type 1";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 line linewidth 2.0";     sendCommand(ss.str()); ss.str(std::string());
    //ss << "s3 line linestyle 5";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 line color 4";           sendCommand(ss.str()); ss.str(std::string());
    //ss << "s3 fill type 2";            sendCommand(ss.str()); ss.str(std::string());
    //ss << "s3 fill rule 0";            sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 fill color 4";           sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 fill pattern 7";         sendCommand(ss.str()); ss.str(std::string());
    //ss << "s3 symbol 10";              sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 symbol size 1.0";        sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 symbol linewidth 2.0";   sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 symbol linestyle 1";     sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 symbol color 4";         sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 symbol pattern 1";       sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 symbol fill color 4";    sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 symbol fill pattern 1";  sendCommand(ss.str()); ss.str(std::string());
    ss << "s3 errorbar linewidth 2.0"; sendCommand(ss.str()); ss.str(std::string());
  }
  
  ss << "redraw";
  sendCommand(ss.str());
  ss.str(std::string());
}


void Grace::duplicate(int set, int Graph) const
{
  std::stringstream s;

  if(Graph >= 0) s << "FOCUS G" << Graph;
  sendCommand(s.str());
  s.str(std::string());

  std::stringstream s1;
  s1 << "DUPLICATE s" << set;
  
  sendCommand(s1.str());
}

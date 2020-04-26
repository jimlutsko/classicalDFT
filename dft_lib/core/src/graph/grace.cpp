#include "dft_lib/graph/grace.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

namespace dft_core {
  namespace grace_plot {

    void send_command(const std::string &cmd) {
      if (GraceIsOpen()) {
        GracePrintf(cmd.c_str());
      }
    }

    void error_function(const char *msg)
    {
      std::cout <<  "Grace library message: \"" << msg << "\"" << std::endl;
    };

    void register_grace_error_function()
    {
      GraceRegisterErrorFunction(error_function);
    }

    void start_grace_communication(const double& x_size, const double& y_size, int buffer_size)
    {
      char grace_name[32];
      sprintf(grace_name,"xmgrace");

      std::string free_opt = "-free";
      std::string no_safe_opt = "-nosafe";
      std::string geometry_opt = "-geometry";

      std::stringstream geometry_spec;
      geometry_spec << x_size << "x" << y_size;

      // Start Grace with a buffer size and open the pipe
      auto response = GraceOpenVA(grace_name, buffer_size, free_opt.c_str(), no_safe_opt.c_str(), geometry_opt.c_str(), geometry_spec.str().c_str(), NULL);
      if (-1 == response) {
        throw std::runtime_error("Can't run Grace. \n");
      }
    }

    // Send some initialization commands to Grace
    void setup_grace(const double& x_max, const double& x_min, const double& y_max, const int& number_of_graphs)
    {
      std::stringstream s;
      if(number_of_graphs > 1)
      {
        int rows = (number_of_graphs <= 2 ? 1 : 2);
        int cols = number_of_graphs / rows;
        cols += (number_of_graphs - rows*cols);

        s << "ARRANGE(" << rows << ", " << cols << ", 0.1, 0.15, 0.2)";
        send_command(s.str());
        s.str(std::string());
      }

      s << "world xmax " << x_max;
      send_command(s.str());
      s.str(std::string());

      s << "world xmin " << x_min;
      send_command(s.str());
      s.str(std::string());

      s << "world ymax "  << y_max;
      send_command(s.str());
      s.str(std::string());

      s << "xaxis tick major 5";
      send_command(s.str());
      s.str(std::string());

      s << "xaxis tick minor 1";
      send_command(s.str());
      s.str(std::string());

      s << "yaxis tick major " << int(y_max);
      send_command(s.str());
      s.str(std::string());

      s << "yaxis tick minor " << int(y_max/2);
      send_command(s.str());
      s.str(std::string());

      s << "AUTOSCALE ONREAD XYAXES";
      send_command(s.str());
      s.str(std::string());
    }

    Grace::Grace(int x_size, int y_size, int n_graph, bool show) :
        x_min_(min_axis_value), x_max_(max_axis_value),
        y_min_(min_axis_value), y_max_(min_axis_value),
        n_max_data_set_(default_dataset_number), n_graph_(n_graph),
        show_(show)
    {
      if (show) {
        register_grace_error_function();
        start_grace_communication(x_size, y_size);
        setup_grace(x_max_, x_min_, x_max_, n_graph_);
      }
    }
  }
}

/*
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
  s1 << "g" << G << ".s" << set << " symbol color " << color;
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
*/
// The order in which headers are imported follows Google's C++ style guide:
// link: https://google.github.io/styleguide/cppguide.html
#include "dft_lib/graph/grace.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

#include "dft_lib/exceptions/grace_exception.h"

namespace dft_core
{
  namespace grace_plot
  {
    void SendCommand(const std::string& cmd)
    {
      if (GraceIsOpen()) {
        try { GracePrintf(cmd.c_str()); }
        catch (...) { throw GraceCommunicationFailedException(); }
      }
      else { throw GraceNotOpenedException(); }
    }

    void ErrorParsingFunction(const char* msg)
    {
      std::cout << "Grace library message: \"" << msg << "\"" << std::endl;
    };

    void RegisterGraceErrorFunction()
    {
      try {
        GraceRegisterErrorFunction(ErrorParsingFunction);
      }
      catch (const std::exception& e) {
        throw GraceException("ErrorParsingFunction could not be registered by xmgrace");
      }
    }

    void StartGraceCommunication(const double& x_size, const double& y_size, int buffer_size)
    {
      if ((x_size <= 0) || (y_size <= 0)) { throw GraceException("The xmgrace-canvas dimensions cannot be negative!"); }
      if (buffer_size <= 0) { throw GraceException("The communication buffer size cannot be negative!"); }

      char grace_name[32];
      sprintf(grace_name,"xmgrace");

      std::string free_opt = "-free";
      std::string no_safe_opt = "-nosafe";
      std::string geometry_opt = "-geometry";

      std::stringstream geometry_spec;
      geometry_spec << x_size << "x" << y_size;

      // Start Grace with a buffer size and open the pipe
      auto response = GraceOpenVA(grace_name, buffer_size, free_opt.c_str(), no_safe_opt.c_str(), geometry_opt.c_str(), geometry_spec.str().c_str(), NULL);

      // Throwing exception in case of wrong communication, to keep track of what's happening
      if (-1 == response) { throw GraceCommunicationFailedException(); }
    }

    // TODO: Refactor this method and add testing
    void SetupGrace(const double& x_max, const double& x_min, const double& y_max, const int& number_of_graphs)
    {
      std::stringstream s;
      if(number_of_graphs > 1)
      {
        int rows = (number_of_graphs <= 2 ? 1 : 2);
        int cols = number_of_graphs / rows;
        cols += (number_of_graphs - rows*cols);

        s << "ARRANGE(" << rows << ", " << cols << ", 0.1, 0.15, 0.2)";
        SendCommand(s.str());
        s.str(std::string());
      }

      s << "world xmax " << x_max;
      SendCommand(s.str());
      s.str(std::string());

      s << "world xmin " << x_min;
      SendCommand(s.str());
      s.str(std::string());

      s << "world ymax "  << y_max;
      SendCommand(s.str());
      s.str(std::string());

      s << "xaxis tick major 5";
      SendCommand(s.str());
      s.str(std::string());

      s << "xaxis tick minor 1";
      SendCommand(s.str());
      s.str(std::string());

      s << "yaxis tick major " << int(y_max);
      SendCommand(s.str());
      s.str(std::string());

      s << "yaxis tick minor " << int(y_max/2);
      SendCommand(s.str());
      s.str(std::string());

      s << "AUTOSCALE ONREAD XYAXES";
      SendCommand(s.str());
      s.str(std::string());
    }

    Grace::Grace(int x_size, int y_size, int n_graph, bool show) :
        x_min_(min_axis_value), x_max_(max_axis_value),
        y_min_(min_axis_value), y_max_(min_axis_value),
        n_max_data_set_(default_dataset_number), n_graph_(n_graph),
        show_(show)
    {
      if (show) {
        RegisterGraceErrorFunction();
        StartGraceCommunication(x_size, y_size);
        SetupGrace(x_max_, x_min_, x_max_, n_graph_);
      }
    }
  }
}

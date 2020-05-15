// The order in which headers are imported follows Google's C++ style guide:
// link: https://google.github.io/styleguide/cppguide.html
#include "dft_lib/graph/grace.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dft_lib/exceptions/grace_exception.h"
#include "dft_lib/utils/console.h"

namespace dft_core
{
  namespace grace_plot
  {
    namespace command
    {
      std::string ArrangeCommand(const int& number_of_rows, const int& number_of_columns, const float& offset, const float& horizontal_gap, const float& vertical_gap)
      {
        std::string cmd = "ARRANGE(" + std::to_string(number_of_rows) + ", "
                          + std::to_string(number_of_columns) + ", "
                          + std::to_string(offset) + ", "
                          + std::to_string(horizontal_gap) + ", "
                          + std::to_string(vertical_gap) + ")";
        return cmd;
      }
      std::string SetXMinCommand(const double& x_min)
      {
        auto x_rounded = static_cast<float>(x_min);
        std::string cmd = "WORLD XMIN " + std::to_string(x_rounded);
        return cmd;
      }
      std::string SetXMaxCommand(const double& x_max)
      {
        auto x_rounded = static_cast<float>(x_max);
        std::string cmd = "WORLD XMAX " + std::to_string(x_rounded);
        return cmd;
      }
      std::string SetYMinCommand(const double& y_min)
      {
        auto y_rounded = static_cast<float>(y_min);
        std::string cmd = "WORLD YMIN " + std::to_string(y_rounded);
        return cmd;
      }
      std::string SetYMaxCommand(const double& y_max)
      {
        auto y_rounded = static_cast<float>(y_max);
        std::string cmd = "WORLD YMAX " + std::to_string(y_rounded);
        return cmd;
      }
      std::string AddPointCommand(const double& x, const double& y, const int& dataset_id, const int& graph_id)
      {
        std::string cmd = "G" + std::to_string(graph_id) + ".S" + std::to_string(dataset_id)
                          + " POINT " + std::to_string(x) + "," + std::to_string(y);
        return cmd;
      }
      std::string RedrawCommand()
      {
        std::string cmd = "REDRAW";
        return cmd;
      }
      std::string AutoScaleCommand()
      {
        std::string cmd = "AUTOSCALE";
        return cmd;
      }
      std::string AutoTicksCommand()
      {
        std::string cmd = "AUTOTICKS";
        return cmd;
      }
      std::string FocusCommand(const int& graph_id)
      {
        std::string cmd = "FOCUS G" + std::to_string(graph_id);
        return cmd;
      }
      std::string KillSetCommand(const int& dataset_id, const int& graph_id)
      {
        std::string cmd = "KILL G" + std::to_string(graph_id) + "." + "S" + std::to_string(dataset_id);
        return cmd;
      }
      std::string SetLegendCommand(const std::string& legend, const int& dataset_id, const int& graph_id)
      {
        std::string cmd = "G" + std::to_string(graph_id)
                          + ".S" + std::to_string(dataset_id)
                          + " LEGEND \"" + legend + "\"";
        return cmd;
      }
      std::string SetLineColorCommand(const grace_plot::Color& color, const int& dataset_id, const int& graph_id)
      {
        std::string cmd = "G" + std::to_string(graph_id)
                          + ".S" + std::to_string(dataset_id)
                          + " LINE COLOR " + std::to_string(static_cast<int>(color));
        return cmd;
      }
      std::string SetSymbolColorCommand(const grace_plot::Color& color, const int& dataset_id, const int& graph_id) {
        std::string cmd = "G" + std::to_string(graph_id)
                          + ".S" + std::to_string(dataset_id)
                          + " SYMBOL COLOR " + std::to_string(static_cast<int>(color));
        return cmd;
      }
      std::string SetAxisLabelCommand(const std::string& label, const grace_plot::Axis& axis)
      {
        std::string axis_string = (axis == Axis::X ? "X" : "Y");
        std::string cmd = axis_string + "AXIS LABEL \"" + label + "\"";
        return cmd;
      }
      std::string SetTitleCommand(const std::string& title)
      {
        std::string cmd = "TITLE \"" + title + "\"";
        return cmd;
      }
      std::string SetSubtitleCommand(const std::string& subtitle)
      {
        std::string cmd = "SUBTITLE \"" + subtitle + "\"";
        return cmd;
      }
    }
  }
}


namespace dft_core
{
  namespace grace_plot
  {
    //region static methods

    static void AddDatasetToGraceObject(const Grace* g, std::vector<double> const& x, std::vector<double> const& y, const int& dataset_id, const int& graph_id)
    {
      unsigned int length = x.size();
      for (uint k = 0; k < length; k++)
      {
        g->AddPoint(x[k], y[k], dataset_id, graph_id);
      }
    }

    static bool CheckEqualLength(std::vector<double> const& x, std::vector<double> const& y)
    {
      if (x.size() != y.size())
      {
        throw exception::GraceException("The dataset must have equal-size X and Y:"
                                        " x_size = " + std::to_string(x.size()) + "!=" +
                                        " y_size = " + std::to_string(y.size()));
      }
      return true;
    }

    static bool CheckGraphIdInBounds(const int& graph_id, const int& total_number_of_graphs)
    {
      if ((graph_id > (total_number_of_graphs - 1)) || (graph_id < 0))
      {
        throw exception::GraceException("The graph id is out of bounds: Min id = 0; Max id =" + std::to_string(total_number_of_graphs));
      }
      return true;
    }

    static bool CheckDatasetInBounds(const int& dataset_id, const int& last_dataset_id)
    {
      if (dataset_id < 0 || dataset_id > last_dataset_id)
      {
        throw exception::GraceException("The dataset id is out of bounds: Min id = 0; Last dataset id =" + std::to_string(last_dataset_id));
      }
      return true;
    }

    //endregion

    //region General methods

    void SendCommand(const std::string& cmd)
    {
      if (GraceIsOpen())
      {
        try { GracePrintf(cmd.c_str()); }
        catch (...) { throw dft_core::exception::GraceCommunicationFailedException(); }
      }
      else { throw dft_core::exception::GraceNotOpenedException(); }
    }

    void ErrorParsingFunction(const char* msg)
    {
      std::cout << "Grace library message: \"" << msg << "\"" << std::endl;
    }

    void RegisterGraceErrorFunction()
    {
      try
      {
        GraceRegisterErrorFunction(ErrorParsingFunction);
      }
      catch (const std::exception& e)
      {
        throw dft_core::exception::GraceException("ErrorParsingFunction could not be registered by xmgrace");
      }
    }

    void StartGraceCommunication(const int& x_size, const int& y_size, int buffer_size)
    {
      if ((x_size <= 0) || (y_size <= 0))
      {
        throw dft_core::exception::GraceException("The xmgrace-canvas dimensions cannot be negative!");
      }
      if (buffer_size <= 0)
      {
        throw dft_core::exception::GraceException("The communication buffer size cannot be negative!");
      }

      char grace_name[32];
      sprintf(grace_name,"xmgrace");

      // Start Grace with a buffer size and open the pipe
      std::string geometry_spec = std::to_string(x_size) + "x" + std::to_string(y_size);
      auto response = GraceOpenVA(grace_name, buffer_size, option::FREE.c_str(), option::NO_SAFE.c_str(), option::GEOMETRY.c_str(), geometry_spec.c_str(), NULL);

      // Throwing exception in case of wrong communication, to keep track of what's happening
      if (-1 == response) { throw dft_core::exception::GraceCommunicationFailedException(); }
    }

    int GetNumberOfRows(const int& number_of_graphs)
    {
      if (number_of_graphs < 1) { throw dft_core::exception::GraceException("Number of graphs cannot be lesser than one"); }
      int n_rows = (number_of_graphs <= 2 ? 1 : 2);
      return n_rows;
    }

    int GetNumberOfColumns(const int& number_of_graphs, const int& number_of_rows)
    {
      if (number_of_graphs < 1) { throw dft_core::exception::GraceException("Number of graphs cannot be lesser than one"); }
      else if (number_of_rows < 1) { throw dft_core::exception::GraceException("Number of rows cannot be lesser than one"); }

      int number_of_columns = number_of_graphs / number_of_rows;
      number_of_columns += (number_of_graphs - number_of_rows * number_of_columns);

      return number_of_columns;
    }
    //endregion

    //region Grace

    /// The configuration command which eventually sets up the graph
    void SetupGrace(
        const double& x_min, const double& x_max,
        const double& y_min, const double& y_max,
        const int& number_of_graphs, const float& offset,
        const float& hspace, float const& vspace)
    {
      if(number_of_graphs > 1)
      {
        int number_of_rows = GetNumberOfRows(number_of_graphs);
        int number_of_columns = GetNumberOfColumns(number_of_graphs, number_of_rows);
        SendCommand(command::ArrangeCommand(number_of_rows, number_of_columns, offset, hspace, vspace));
      }

      SendCommand(command::SetXMinCommand(x_min));
      SendCommand(command::SetXMaxCommand(x_max));

      // tick major must be set up befor minor to avoid glitch
      SendCommand( "XAXIS TICK MAJOR 5");
      SendCommand( "XAXIS TICK MINOR 1");

      SendCommand("WORLD YMIN " + std::to_string(y_min));
      SendCommand("WORLD YMAX " + std::to_string(y_max));

      // tick major must be set up befor minor to avoid glitch
      SendCommand( "YAXIS TICK MAJOR " + std::to_string(static_cast<int>(y_max)));
      SendCommand( "YAXIS TICK MINOR " + std::to_string(static_cast<int>(y_max/2)));

      SendCommand( "AUTOSCALE ONREAD XYAXES");
    }

    Grace::Grace(int x_size, int y_size, int n_graph, bool show) :
        x_min_(default_min_axis_value), x_max_(default_max_axis_value),
        y_min_(default_min_axis_value), y_max_(default_max_axis_value),
        offset_(default_offset),
        horizontal_space_(default_horizontal_space),
        vertical_space_(default_vertical_space),
        last_dataset_id_(default_dataset_id), number_of_graphs_(n_graph),
        show_(show)
    {
      if (show_) {
        RegisterGraceErrorFunction();
        StartGraceCommunication(x_size, y_size);
        SetupGrace(x_min_, x_max_, y_min_, y_max_, number_of_graphs_, offset_, horizontal_space_, vertical_space_);
      }
    }

    void Grace::SetXMin(const double& value)
    {
      x_min_ = value;
    }

    void Grace::SetXMax(const double& value)
    {
      x_max_ = value;
    }

    void Grace::SetYMin(const double& value)
    {
      y_min_ = value;
    }

    void Grace::SetYMax(const double& value)
    {
      y_max_ = value;
    }

    void Grace::SetXLimits(const double& x_min, const double& x_max)
    {
      if (x_min > x_max) { throw exception::GraceException("Lower limit cannot be greater than upper limit!"); }

      this->SetXMin(x_min);
      this->SetXMax(x_max);

      if (this->show_) {
        SendCommand(command::SetXMinCommand(this->x_min()));
        SendCommand(command::SetXMaxCommand(this->x_max()));
      }
    }

    void Grace::SetYLimits(const double& y_min, const double& y_max)
    {
      if (y_min > y_max) { throw exception::GraceException("Lower limit cannot be greater than upper limit!"); }

      this->SetYMin(y_min);
      this->SetYMax(y_max);

      if (this->show_)
      {
        SendCommand(command::SetYMinCommand(this->y_min()));
        SendCommand(command::SetYMaxCommand(this->y_max()));
      }
    }

    void Grace::SetLimits(const double &x_min, const double &x_max, const double &y_min, const double &y_max)
    {
      this->SetXLimits(x_min, x_max);
      this->SetYLimits(y_min, y_max);
    }

    void Grace::SetLimits(const std::vector<double>& x_limits, const std::vector<double>& y_limits)
    {
      this->SetLimits(x_limits.front(), x_limits.back(), y_limits.front(), y_limits.back());
    }

    void Grace::Close() const
    {
      if (this->show_)
      {
        GraceClose();
      }
    }

    void Grace::AddPoint(const double &x, const double &y, const int &dataset_id, const int &graph_id) const
    {
      CheckGraphIdInBounds(graph_id, this->number_of_graphs());

      if (this->show_)
      {
        SendCommand(command::AddPointCommand(x, y, dataset_id, graph_id));
      }
    }

    void Grace::IncreaseLastDatasetId()
    {
      this->last_dataset_id_ += 1;
    }

    void Grace::DecreaseLastDatasetId()
    {
      this->last_dataset_id_ -= 1;
    }

    int Grace::AddDataset(std::vector<double> const& x, std::vector<double> const& y, const int& graph_id)
    {
      CheckEqualLength(x, y);
      this->IncreaseLastDatasetId();
      AddDatasetToGraceObject(this, x, y, this->last_dataset_id(), graph_id);
      return last_dataset_id_;
    }

    void Grace::DeleteDataset(const int& dataset_id, const int& graph_id)
    {
      CheckDatasetInBounds(dataset_id, this->last_dataset_id());
      CheckGraphIdInBounds(graph_id, this->number_of_graphs());

      SendCommand(command::KillSetCommand(dataset_id, graph_id));
      //this->DecreaseLastDatasetId();

      this->SetColor(static_cast<Color>((this->last_dataset_id() % 10) + 1), dataset_id, graph_id);
    }

    void Grace::ReplaceDataset(std::vector<double> const& x, std::vector<double> const& y, const int& dataset_id, const int& graph_id)
    {
      CheckEqualLength(x, y);
      CheckDatasetInBounds(dataset_id, this->last_dataset_id());
      CheckGraphIdInBounds(graph_id, this->number_of_graphs());

      this->DeleteDataset(dataset_id, graph_id);
      AddDatasetToGraceObject(this, x, y, dataset_id, graph_id);
    }

    void Grace::Redraw(const bool& auto_scale, const bool& auto_ticks, const int& graph_id) const
    {
      if (this->show_)
      {
        if (CheckGraphIdInBounds(graph_id, this->number_of_graphs()))
        {
          SendCommand(command::FocusCommand(graph_id));
        }

        if (auto_scale) {
          SendCommand(command::AutoScaleCommand());
        }

        if (auto_ticks) {
          SendCommand(command::AutoTicksCommand());
        }

        SendCommand(command::RedrawCommand());
      }
    }

    void Grace::Wait() const
    {
      console::wait();
    }

    void Grace::RedrawAndWait(const bool& auto_scale, const bool& auto_ticks, const int& graph_id) const
    {
      this->Redraw(auto_scale, auto_ticks, graph_id);
      this->Wait();
    }

    void Grace::SetLegend(const std::string& legend, const int& dataset_id, const int& graph_id) const
    {
      CheckDatasetInBounds(dataset_id, this->last_dataset_id());
      CheckGraphIdInBounds(graph_id, this->number_of_graphs());

      SendCommand(command::SetLegendCommand(legend, dataset_id, graph_id));
    }

    void Grace::SetColor(const Color& color, const int& dataset_id, const int& graph_id) const
    {
      CheckDatasetInBounds(dataset_id, this->last_dataset_id());
      CheckGraphIdInBounds(graph_id, this->number_of_graphs());

      SendCommand(command::SetLineColorCommand(color, dataset_id, graph_id));
      SendCommand(command::SetSymbolColorCommand(color, dataset_id, graph_id));
    }

    void Grace::SetLabel(const std::string &label, const Axis &axis, const int& graph_id) const
    {
      CheckGraphIdInBounds(graph_id, this->number_of_graphs());

      if (graph_id > 0) { SendCommand(command::FocusCommand(graph_id)); }
      SendCommand(command::SetAxisLabelCommand(label, axis));
    }

    void Grace::SetTitle(const std::string &title) const
    {
      SendCommand(command::SetTitleCommand(title));
    }

    void Grace::SetSubtitle(const std::string &subtitle) const
    {
      SendCommand(command::SetSubtitleCommand(subtitle));
    }
    //endregion
  }
}

#ifndef __CLASSICALDFT_GRACE__
#define __CLASSICALDFT_GRACE__

#include <unistd.h>
#include <grace_np.h>

#include <string>
#include <vector>

namespace dft_core {

/**
  *  @brief  UTILITY: wrapper for xmgrace graphics program.
  *
  *  @detailed This class provides an interface for the xmgrace graphics program which can be called programatically to provide a cheap solution  for displaying dynamically updated line graphs etc.
  */
  namespace grace_plot {

    enum ExportFormat{
      JPG,
      PNG,
      EPS1,
      PS,
      PDF
    };

    enum Shape {
      CIRCLE,
      SQUARE
    };

    // Default configuration values
    const int default_x_size = 800;
    const int default_y_size = 600;
    const int default_dataset_number = 0;
    const int default_number_of_graphs = 1;
    const auto min_axis_value = 0.0;
    const auto max_axis_value = 10.0;

    void send_command(const std::string &cmd);
    void error_function(const char *msg);
    void register_grace_error_function();
    void start_grace_communication(const double& x_size, const double& y_size, int buffer_size = 2048);

    class Grace {
    public:
      explicit Grace(int x_size = default_x_size, int y_size = default_y_size, int n_graph = default_number_of_graphs, bool show = true);
      ~Grace() = default;

      bool is_initialised() { return show_; }
      double get_x_min() { return x_min_; }
      double get_x_max() { return x_max_; }
      double get_y_min() { return x_min_; }
      double get_y_max() { return x_max_; }

      /*
      void Close() { GraceClose(); }

      void Store() const;

      void Store(std::string &file) const;

      void Kill() const { sendCommand("Kill G0"); }

      void AddPoint(double x, double y, int N = 0, int G = 0) const;

      int AddDataSet(std::vector<double> const &x, std::vector<double> const &y);

      void ReplaceDataSet(std::vector<double> const &x, std::vector<double> const &y, int N = 0) const;

      void DeleteDataSet(int N = 0, int G = 0) const;

      void SetLegend(const char *s, int n);

      void SetColor(int color, int nm, int G = 0) const;

      void SetXAxisLabel(const char *s, int Graph = -1);

      void SetYAxisLabel(const char *s, int Graph = -1);

      void SetTitle(const char *s);

      void SetSubTitle(const char *s);

      void Redraw(int autoScale = 1, int Graph = -1) const;

      void SetLimits(double xmin, double xmax, double ymin, double ymax, int Graph = -1);

      void SetTicks(double dx, double dy, int Graph = -1);

      void SetSymbol(int set, int symbol, int Graph = -1);

      void SetSymbolFill(int set, int color, int Graph = -1);

      void Pause() const {
        char cc;
        std::cout << "Enter any character to continue: ";
        std::cin >> cc;
      }

      void NoLine(int dataSet, int Graph = -1);

      void Line(int dataSet, int Graph = -1);

      void Circle(int dataSet, int Graph = -1);

      void Square(int dataSet, int Graph = -1);

      void Size(int dataSet, double size, int Graph = -1);

      static int CIRCLE;
      static int SQUARE;

      void PrintToFile(std::string &file, int format);

      void AutoTic(int Graph = -1) const;
    */
    private:
      double x_min_ = min_axis_value;
      double x_max_ = max_axis_value;
      double y_min_ = min_axis_value;
      double y_max_ = max_axis_value;
      int n_max_data_set_ = default_dataset_number;
      int n_graph_ = default_number_of_graphs;
      bool show_ = false;

    };
  }
}

#endif

#ifndef CLASSICALDFT_GRACE_H
#define CLASSICALDFT_GRACE_H

#include <unistd.h>
#include <grace_np.h>

#include <string>
#include <vector>

namespace dft_core
{
  namespace grace_plot
  {
    /**
     * @brief Supported file formats when saving a graph
     *
     * These values are used in arguments to various functions in this package.
     */
    enum ExportFormat {
      JPG,
      PNG,
      EPS1,
      PS,
      PDF
    };

    /**
     * @brief Supported shapes for the data points
     */
    enum Shape {
      CIRCLE,
      SQUARE
    };

    namespace option
    {
      /// Option to specify free page layout
      const std::string FREE = "-free";
      /// Option to disable safe mode
      const std::string NO_SAFE = "-nosafe";
      /// Option to
      const std::string GEOMETRY = "-geometry";
    }

    /**
     * @brief Collection of wrappers for xmgrace commands (see https://plasma-gate.weizmann.ac.il/Grace/doc/UsersGuide.html)
     */
    namespace command
    {
      /**
       * @brief Returns the ARRANGE(nrows, ncols, offset, hgap, vgap) command as string
       *
       * @param number_of_rows the number of rows to be used
       * @param number_of_columns the number of columns to be used
       * @param offset the space left at each page edge with
       * @param horizontal_gap horizontal spacing
       * @param vertical_gap vertical spacing
       */
      std::string Arrange(const int& number_of_rows, const int& number_of_columns, const float& offset, const float& horizontal_gap, const float& vertical_gap);

      /**
       * @brief Returns the "WORLD XMIN x_min" command as string
       *
       * @param x_min the minimum value the X-axis will show
       */
      std::string SetXMinCommand(const double& x_min);

      /**
       * @brief Returns the "WORLD XMAX x_min" command as string
       *
       * @param x_max the max value the X-axis will show
       */
      std::string SetXMaxCommand(const double& x_max);

      /**
       * @brief Returns the "WORLD YMIN y_min" command as string
       *
       * @param x_min the minimum value the Y-axis will show
       */
      std::string SetYMinCommand(const double& y_min);

      /**
       * @brief Returns the "WORLD YMAX y_min" command as string
       *
       * @param x_max the max value the Y-axis will show
       */
      std::string SetYMaxCommand(const double& y_max);
    }

    /// The default X-size of the grace canvas
    const int default_x_size = 800;
    /// The default Y-size of the grace canvas
    const int default_y_size = 600;

    const int default_dataset_number = 0;
    /// The default number of graphs to be drawn on the same canvas
    const int default_number_of_graphs = 1;
    /// The value to be used as the default minimum of an axis
    const auto default_min_axis_value = 0.0;
    /// The value to be used as the default maximum of an axis
    const auto default_max_axis_value = 10.0;

    /// The default offset to be used when setting up the appearance of the graphs
    const auto default_offset = 0.1;
    /// The default relative horizontal space to be left
    const auto default_horizontal_space = 0.15;
    /// The default relative vertical space to be left
    const auto default_vertical_space = 0.2;

    /**
     * @brief Sends a string command to the Grace CLI
     *
     * Helper function which sends a command string to xmgrace by using `GracePrintf` from `grace_np.h`.
     * @param[in] cmd the command string to be run in xmgrace
     * @throw GraceNotOpenedException when the subprocess `xmgrace` is not running, i.e. Grace is not "opened".
     * @throw GraceCommunicationFailedException when the subprocess `xmgrace` is running but something happened inside
     */
    void SendCommand(const std::string &cmd);

    /**
     * @brief Parse and prints an error message from xmgrace.
     *
     * The library `grace_np.h` allows for registering a personalised function which will act as an error parser.
     * @param[in] msg the message (char*) received from `xmgrace` which will be printed out.
     * @remark Unfortunately, due to legacy reasons, such a function must deal with char* instead of modern strings.
     */
    void ErrorParsingFunction(const char *msg);

    /**
     * @brief Set up of the personalised error parsing function `ErrorParsingFunction`
     *
     * It uses `GraceRegisterErrorFunction` from `gnuplot_np.h` to register our personalised `ErrorParsingFunction`
     * as the error function `xmgrace` will communicate to when an error occurs.
     * @throw GraceException when there is an error registering the error function.
     */
    void RegisterGraceErrorFunction();

    /**
     * @brief Initialisation of the communication pipe with `xmgrace`
     *
     * This is the first mandatory step to start the communication pipeline with `xmgrace`. This method utilises the
     * `GraceOpenVA()` method from `grace_np.h`, which basically requires as first argument the name of the `xmgrace`
     * executable, a `buffer_size` and a variable number of options to be passed to `xmgrace`. This method will only
     * pass the following parameters:
     * @param x_size the x-size of the canvas, must be positive (default = `default_x_size`)
     * @param y_size the y-size of the canvas, must be positive (default = `default_y_size`)
     * @param buffer_size the buffer size, must be positive (default = 2048)
     * @throw GraceException when any of the parameters are not strictly positive
     * @throw GraceCommunicationFailedException  when something goes wrong in the communication with the communication
     */
    void StartGraceCommunication(const int& x_size, const int& y_size, int buffer_size = 2048);

    /**
     * @brief Returns the number of rows in which the graphs will be placed
     *
     * This simple method just returns 1 if the `number_of_graphs` = 1, and 2 if `number_of_graphs` > 1
     * @param number_of_graphs the number of graphs to be shown
     * @throw GraceException when any of the parameters are not strictly positive
     */
    int GetNumberOfRows(const int& number_of_graphs);

    /**
     * @brief Returns the number of columns in which the graphs will be placed
     *
     * @param number_of_graphs the number of graphs to be shown
     * @param number_of_rows the number of rows to use
     * @throw GraceException when any of the parameters are not strictly positive
     */
    int GetNumberOfColumns(const int& number_of_graphs, const int& number_of_rows);

    /**
     * @brief  Utility: Wrapper for xmgrace graphics program.
     * @detailed This class provides an interface for the xmgrace graphics program which can be called programmatically
     * to provide a cheap solution  for displaying dynamically updated line graphs etc.
    */
    class Grace
    {
      public:
        /// Explicit constructor, which avoids implicit conversions
        explicit Grace(int x_size = default_x_size, int y_size = default_y_size, int n_graph = default_number_of_graphs, bool show = true);
        /// Default destructor
        ~Grace() { GraceClose(); };

        /// Inspector of the property `show_` which governs the construction of a Grace object in case of being `false`
        const bool& is_initialised() const { return show_; }
        /// The minimum value for the X axis
        const double& x_min() const { return x_min_; }
        /// The maximum value for the X axis
        const double& x_max() const { return x_max_; }
        /// The minimum value for the Y axis
        const double& y_min() const { return x_min_; }
        /// The maximum value for the Y axis
        const double& y_max() const { return x_max_; }

        /// The offset
        const float& offset() const { return offset_; }
        /// The hspace to be used when setting up the graphs
        const float& horizontal_space() const { return horizontal_space_; }
        /// The vspace to be used when setting up the graphs
        const float& vertical_space() const { return vertical_space_; }

        /// Setters:
        void SetXLimits(const double& x_min, const double& x_max) const;
        void SetYLimits(const double& y_min, const double& y_max) const;

      private:
        double x_min_ = default_min_axis_value;
        double x_max_ = default_max_axis_value;
        double y_min_ = default_min_axis_value;
        double y_max_ = default_max_axis_value;

        int n_max_data_set_ = default_dataset_number;
        int n_graph_ = default_number_of_graphs;
        bool show_ = false;

        float offset_ = default_offset;
        float horizontal_space_ = default_horizontal_space;
        float vertical_space_ = default_vertical_space;
    };
  }
}

#endif

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
    enum Symbol {
      CIRCLE=1,
      SQUARE,
      DIAMOND,
      TRIANGLE_UP,
      TRIANGLE_LEFT,
      TRIANGLE_DOWN,
      TRIANGLE_RIGHT,
      PLUS,
      CROSS,
      STAR,
    };

    /**
     * @brief The possible colors we can use
     */
    enum Color {
      WHITE = 0,
      BLACK,
      RED,
      GREEN,
      BLUE,
      YELLOW,
      BROWN,
      GREY,
      VIOLET,
      CYAN,
      MAGENTA,
      ORANGE,
      INDIGO,
      MAROON,
      TURQUOISE,
      DARKGREEN
    };

    enum Axis
    {
      X,
      Y
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
       * @return std::string
       */
      std::string ArrangeCommand(const int& number_of_rows, const int& number_of_columns, const float& offset, const float& horizontal_gap, const float& vertical_gap);

      /**
       * @brief Returns the "WORLD XMIN x_min" command as string
       *
       * @param x_min the minimum value the X-axis will show
       * @return std::string
       */
      std::string SetXMinCommand(const double& x_min);

      /**
       * @brief Returns the "WORLD XMAX x_min" command as string
       *
       * @param x_max the max value the X-axis will show
       * @return std::string
       */
      std::string SetXMaxCommand(const double& x_max);

      /**
       * @brief Returns the "WORLD YMIN y_min" command as string
       *
       * @param x_min the minimum value the Y-axis will show
       * @return std::string
       */
      std::string SetYMinCommand(const double& y_min);

      /**
       * @brief Returns the "WORLD YMAX y_min" command as string
       *
       * @param x_max the max value the Y-axis will show
       * @return std::string
       */
      std::string SetYMaxCommand(const double& y_max);

      /**
       * @brief Returns the "G{N}.S{M} point {X},{Y} command as string
       *
       * @param x the x-coordinate of the point
       * @param y the y-coordinate of the point
       * @param dataset_id the integer number identifying the dataset the point will be associated top
       * @param graph_id the integer number identifying the graph the point will be represented on
       * @return std::string
       */
      std::string AddPointCommand(const double& x, const double& y, const int& dataset_id, const int& graph_id);

      /**
        * @brief Returns the "REDRAW" command as string
        * @return std::string
        */
      std::string RedrawCommand();

      /**
       * @brief Returns the "AUTOSCALE" command as string
       * @return std::string
       */
      std::string AutoScaleCommand();

      /**
       * @brief Returns the "AUTOTICKS" command as string
       * @return std::string
       */
      std::string AutoTicksCommand();

      /**
       * @brief Returns the "FOCUS G{N}" command as string
       * @param graph_id the graph int identifier to focus
       * @return std::string
       */
      std::string FocusCommand(const int& graph_id);

      /**
       * @brief Returns the "KILL G{N}.S{M}" command as string
       *
       * @param dataset_id the dataset int identifier to remove
       * @param graph_id the graph int identifier where the dataset lives in
       * @return std::string
       */
      std::string KillSetCommand(const int& dataset_id, const int& graph_id);

      /**
       * @brief Returns the "S LEGEND `legend`" command as string
       * @param legend the legend to be set
       * @return std::string
       */
      std::string SetLegendCommand(const std::string& legend, const int& dataset_id, const int& graph_id);

      /**
       * @brief Returns the "G{N}.S{M} LINE COLOR {ID}" command as string
       * @param color one of the possible grace_plot::Color values
       * @param dataset_id the integer number identifying the dataset the point will be associated top
       * @param graph_id the integer number identifying the graph the point will be represented on
       * @return std::string
       */
      std::string SetLineColorCommand(const grace_plot::Color& color, const int& dataset_id, const int& graph_id);

      /**
       * @brief Returns the "G{N}.S{M} SYMBOL COLOR {ID}" command as string
       * @param color one of the possible grace_plot::Color values
       * @param dataset_id the integer number identifying the dataset the point will be associated top
       * @param graph_id the integer number identifying the graph the point will be represented on
       * @return std::string
       */
      std::string SetSymbolColorCommand(const grace_plot::Color& color, const int& dataset_id, const int& graph_id);

      /**
       * @brief Returns the "{XY}AXIS LABEL {TEXT}" command as string
       * @param label the text to be set as label of the axis
       * @param axis the `grace_plot::Axis` where the text will serve as label, either X or Y
       * @return std::string
       */
      std::string SetAxisLabelCommand(const std::string& label, const grace_plot::Axis& axis);

      /**
       * @brief Returns the "TITLE '{TEXT}'" command as string
       * @param title the text to be set as the title of the graph
       * @return std::string
       */
      std::string SetTitleCommand(const std::string& title);

      /**
       * @brief Returns the "SUBTITLE '{TEXT}'" command as string
       * @param title the text to be set as the subtitle of the graph
       * @return std::string
       */
      std::string SetSubtitleCommand(const std::string& subtitle);

      /**
       * @brief Returns the "{XY}AXIS TICK MAJOR {TICK_SIZE}" command as string
       * @param tick_sep the space separation between the different ticks
       * @return std::string
       */
      std::string SetTicksCommand(const double& tick_sep, const Axis& axis);


      /**
       * @brief Returns the "G{N}.S{M} SYMBOL {ID}" command as string
       * @param symbol_id one of the possible symbol values
       * @param dataset_id the integer number identifying the dataset the point will be associated top
       * @param graph_id the integer number identifying the graph the point will be represented on
       * @return std::string
       */
      std::string SetSymbolCommand(const Symbol& symbol, const int& dataset_id, const int& graph_id);
    }

    /// The default X-size of the grace canvas
    const int default_x_size = 800;
    /// The default Y-size of the grace canvas
    const int default_y_size = 600;
    /// The default dataset id
    const int default_dataset_id = 0;
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
      ~Grace() { this->Close(); };

      /// Inspector of the property `show_` which governs the construction of a Grace object in case of being `false`
      const bool& is_initialised() const { return show_; }
      /// The minimum value for the X axis
      const double& x_min() const { return x_min_; }
      /// The maximum value for the X axis
      const double& x_max() const { return x_max_; }
      /// The minimum value for the Y axis
      const double& y_min() const { return y_min_; }
      /// The maximum value for the Y axis
      const double& y_max() const { return y_max_; }
      /// The total number of graphs
      const int& number_of_graphs() const { return number_of_graphs_; }
      /// The integer id of the last dataset
      const int& last_dataset_id() const { return last_dataset_id_; }

      /// The offset
      const float& offset() const { return offset_; }
      /// The hspace to be used when setting up the graphs
      const float& horizontal_space() const { return horizontal_space_; }
      /// The vspace to be used when setting up the graphs
      const float& vertical_space() const { return vertical_space_; }

      //region Setters

      void SetXMin(const double& value);
      void SetXMax(const double& value);
      void SetYMin(const double& value);
      void SetYMax(const double& value);

      //region SetLimits

      void SetXLimits(const double& x_min, const double& x_max);
      void SetYLimits(const double& y_min, const double& y_max);
      void SetLimits(const double& x_min, const double& x_max, const double& y_min, const double& y_max);
      void SetLimits(const std::vector<double>& x_limits, const std::vector<double>& y_limits);

      //endregion

      //endregion


      //region Methods

      /**
       * @brief Adds a point to a certain dataset (`dataset_id`) for given graph (`graph_id`)
       *
       * @param x the x-coordinate of the point
       * @param y the y-coordinate of the point
       * @param dataset_id the integer number identifying the dataset the point will be associated top
       * @param graph_id the integer number identifying the graph the point will be represented on (default=0)
       * @throw GraceException when the graph_id given is out of bounds
       */
      void AddPoint(const double& x, const double& y, const int& dataset_id=0, const int& graph_id=0) const;

      /**
       * @brief Adds a tuple (X, Y) of `std::vector` objects which represent a dataset
       *
       * @param x a std::vector of X double-like values
       * @param y a std:vector of Y double-like values
       * @param graph_id the integer number identifying the graph the dataset will be represented on (default=0)
       * @returns the integer identifying the dataset within the graph
       * @throw GraceException when the dataset size is not well-balanced, i.e. x.size() != y.size()
       * @throw GraceException when the graph_id given is out of bounds
       */
      int AddDataset(std::vector<double> const &x, std::vector<double> const &y, const int& graph_id = 0);

      /**
       * @brief Replace an existing dataset with a tuple (X, Y) of `std::vector` objects which represent a new dataset
       *
       * @param x a std::vector of X double-like values
       * @param y a std:vector of Y double-like values
       * @param dataset_id the integer identifier of the dataset to be replaced with new data
       * @param graph_id the integer number identifying the graph the dataset will be represented on (default=0)
       * @returns the integer identifying the dataset within the graph
       * @throw GraceException when the dataset size is not well-balanced, i.e. x.size() != y.size()
       * @throw GraceException when the graph_id given is out of bounds
       */
      void ReplaceDataset(std::vector<double> const &x, std::vector<double> const &y, const int& dataset_id, const int& graph_id = 0);

      /**
       * @brief Deletes a dataset identified by `dataset_id` on the graph identified by `graph_id`
       *
       * @param dataset_id the integer identifier of the dataset to be deleted
       * @param graph_id the integer number identifying the graph where the dataset lives (default=0)
       * @throws GraceException in case the dataset_id given is out of bounds
       * @throws GraceException in case the graph_id given is out of bounds
       */
      void DeleteDataset(const int& dataset_id, const int& graph_id = 0);

      /**
       * @brief Closes the opened session (if it was opened with the constructor, via `show=true`
       */
      void Close() const;

      /**
       * @brief Pauses the terminal awaiting for the user interaction
       */
      void Wait() const;

      /**
       * @brief Updates the graph by using the `redraw` command
       *
       * @param auto_scale if true will autoscale the graph when redrawing
       * @param auto_ticks if true will set automatically the best-fit ticks for X and Y axes
       * @param graph_id the integer number identifying the graph the point will be represented on (default=0)
       * @throws GraceException in case the graph_id given is out of bounds
       */
      void Redraw(const bool& auto_scale = false, const bool& auto_ticks = true, const int& graph_id = 0) const;

      /**
       * @brief Wrapper of Redraw and Wait functionality to play together for convenience
       *
       * @param auto_scale if true will autoscale the graph when redrawing
       * @param auto_ticks if true will set automatically the best-fit ticks for X and Y axes
       * @param graph_id the integer number identifying the graph the point will be represented on (default=0)
       * @throws GraceException in case the graph_id given is out of bounds
       */
      void RedrawAndWait(const bool& auto_scale = false, const bool& auto_ticks = true, const int& graph_id = 0) const;

      /**
       * @brief Sets the legend of the graph
       *
       * @param legend the string to be set as legend of the graph
       * @param dataset_id the integer identifier of the dataset to be deleted
       * @param graph_id the integer number identifying the graph where the dataset lives (default=0)
       * @throws GraceException in case the dataset_id given is out of bounds
       * @throws GraceException in case the graph_id given is out of bounds
       */
      void SetLegend(const std::string& legend, const int& dataset_id, const int& graph_id = 0) const;

      /**
       * @brief Sets a given `dataset_id` of a given `graph_id` with the specified `color_id`
       *
       * @param color the integer identifying the color
       * @param dataset_id the integer identifier of the dataset to be deleted
       * @param graph_id the integer number identifying the graph where the dataset lives (default=0)
       * @throws GraceException in case the dataset_id given is out of bounds
       * @throws GraceException in case the graph_id given is out of bounds
       */
      void SetColor(const Color& color, const int& dataset_id, const int& graph_id = 0) const;

      /**
       * @brief Sets the Axis label of a given `graph_id`
       *
       * @param label the text which will be set as label of the axis
       * @param axis either Axis::X or Axis::Y
       * @param graph_id the integer number identifying the graph where the dataset lives (default=0)
       * @throws GraceException in case the graph_id given is out of bounds
       */
      void SetLabel(const std::string& label, const Axis& axis, const int& graph_id = 0) const;

      /**
       * @brief Sets the graph's Title
       *
       * @param title the text which will be set as the title of the graph
       */
      void SetTitle(const std::string& title) const;

      /**
       * @brief Sets the graph's Subitle
       *
       * @param title the text which will be set as the subtitle of the graph
       */
      void SetSubtitle(const std::string& subtitle) const;

      /**
       * @brief Sets the tick distancing on the X and Y axis
       * @param dx the spacing between X-axis ticks
       * @param dy the spacing between Y-axis ticks
       * @param graph_id
       */
      void SetTicks(const double& dx, const double& dy, const int& graph_id = 0) const;

      /**
       * @brief Sets the symbol shape of given `dataset_id` of a given `graph_id`
       *
       * @param symbol_id one of the possible values of enum::Symbol
       * @param dataset_id the integer identifier of the dataset to be deleted
       * @param graph_id the integer number identifying the graph where the dataset lives (default=0)
       * @throws GraceException in case the dataset_id given is out of bounds
       * @throws GraceException in case the graph_id given is out of bounds
       */
      void SetSymbol(const Symbol& symbol, const int& dataset_id, const int& graph_id = 0) const;

      //endregion

    private:
      double x_min_ = default_min_axis_value;
      double x_max_ = default_max_axis_value;
      double y_min_ = default_min_axis_value;
      double y_max_ = default_max_axis_value;

      int last_dataset_id_ = default_dataset_id;
      int number_of_graphs_ = default_number_of_graphs;
      bool show_ = false;

      float offset_ = default_offset;
      float horizontal_space_ = default_horizontal_space;
      float vertical_space_ = default_vertical_space;

      /**
       * @brief Adds the unit to the last dataset id registered to keep track of how many are active
       */
      void IncreaseLastDatasetId();
      /**
       * @brief Removes the unit to the last dataset id registered to keep track of how many are active
       */
      void DecreaseLastDatasetId();
    };
  }
}

#endif

#ifndef CLASSICALDFT_PARAMETER_EXCEPTION_H
#define CLASSICALDFT_PARAMETER_EXCEPTION_H

#include <string>
#include <utility>

namespace dft_core
{
  namespace exception
  {

    /**
     * @brief Exceptions related with input parameters of the library
     *
     * This class will serve as the base class for the exceptions we want to register regarding undesired or uncontrolled
     * behaviours when dealing with wrong input parameters in methods or constructors
     */
    class WrongParameterException: public std::exception
    {
      private:
        /// The error message which will be deliver when the exception is thrown
        std::string error_message_;

      protected:
        /// A setter of the error message for internal use in constructors of derived classes
        void set_error_message(const std::string& message) { error_message_ = message; }

      public:
        /// Explicit constructor of the exception, to avoid unintentional implicit conversions
        explicit WrongParameterException(std::string msg): error_message_(std::move(msg)) {}

        /// Standard default destructor to be used, it should not throw an exception and if it does the program just crash
        ~WrongParameterException() noexcept override = default;

        /// Overrides the `what` method which is used when catching the message of a `std::exception`
        virtual const char* what() const noexcept override { return error_message_.c_str(); }

        /// Gets the error message as a constant reference to the private field `error_message_`
        const std::string& error_message() { return error_message_; }
    };

    /**
     * @brief Exception to be thrown when detecting a negative parameter is given but positive is the only option
     *
     * This class is just a derived class from the general `WrongParameterException` which particularises the message to
     * specify that the input parameter is negative when it shouldn't be.
     */
    class NegativeParameterException: public WrongParameterException
    {
      public:
        explicit NegativeParameterException(std::string msg);
    };
  }
}

#endif

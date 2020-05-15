#include "classical_dft"

int main(int argc, char **argv)
{
  console::info("Hello world");
  console::warning("This is a warning");
  console::error("This is an error!!");
  console::debug("This is a debugging message");
  console::write_line(console::format::bold("Bold text"));
  console::write_line(console::format::blink("Blinking"));
  console::wait();
}
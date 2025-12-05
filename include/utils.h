#include <iostream>
#include <string>

void
removeSuffix(std::string & str, const std::string & suffix)
{
  if (str.length() >= suffix.length() && str.substr(str.length() - suffix.length()) == suffix)
  {
    str.erase(str.length() - suffix.length());
  }
}
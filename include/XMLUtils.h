#pragma once
#include "tinyxml2.h"
#include <stdexcept>
#include <string>
#include <type_traits>

template <typename T>
T
getAttributeOrThrow(const tinyxml2::XMLElement * elem, const char * name)
{
  auto * attr = elem->FindAttribute(name);
  if (!attr)
    throw std::runtime_error(std::string("Missing attribute: ") + name);

  if constexpr (std::is_integral<T>::value)
  {
    if constexpr (std::is_unsigned<T>::value)
      return attr->UnsignedValue();
    else
      return attr->IntValue();
  }
  else if constexpr (std::is_floating_point<T>::value)
  {
    return attr->DoubleValue();
  }
  else
  {
    static_assert(!sizeof(T), "Unsupported type for getAttributeOrThrow");
  }
};
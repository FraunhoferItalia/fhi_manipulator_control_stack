/*
 * Copyright 2022-2024 Fraunhofer Italia Research
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <fhi_manipulator_control/parameters.hpp>

namespace fmc
{
namespace parameters
{

Parameter::Parameter(const Description & parameter_description)
: name{parameter_description.name},
  type{parameter_description.type},
  lower_limit{parameter_description.lower_limit},
  upper_limit{parameter_description.upper_limit}
{
  switch (type) {
    case Types::BOOL:
      set_value(parameter_description.bool_value);
      break;
    case Types::DOUBLE:
      set_value(parameter_description.double_value);
      break;
    default:
      break;
  }
}

double Parameter::get_double_value() { return double_value; }
bool Parameter::get_bool_value() { return bool_value; }
bool Parameter::set_value(const double & value)
{
  if (value < lower_limit || value > upper_limit) return false;
  double_value = value;
  return true;
}
bool Parameter::set_value(const bool & value)
{
  bool_value = value;
  return true;
}
}  // namespace parameters
}  // namespace fmc

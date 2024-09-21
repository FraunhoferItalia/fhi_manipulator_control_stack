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
#pragma once
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

namespace fmc
{
namespace parameters
{

enum Types { BOOL, DOUBLE };
struct Description
{
  Description(
    const std::string & param_name, const double & param_initial_value,
    const double & param_lower_limit = std::numeric_limits<double>::min(),
    const double & param_upper_limit = std::numeric_limits<double>::max())
  : name{param_name},
    type{Types::DOUBLE},
    double_value{param_initial_value},
    lower_limit{param_lower_limit},
    upper_limit{param_upper_limit} {};

  Description(const std::string & param_name, const bool & param_initial_value)
  : name{param_name},
    type{Types::BOOL},
    bool_value{param_initial_value},
    lower_limit{0.0},
    upper_limit{1.0} {};

  const std::string name;
  const Types type;
  const bool bool_value{false};
  const double double_value{0.0};
  const double lower_limit;
  const double upper_limit;
};

class Parameter
{
public:
  using SharedPtr = std::shared_ptr<Parameter>;

  Parameter(const Description & parameter_description);

  double get_double_value();
  bool get_bool_value();
  bool set_value(const double & value);
  bool set_value(const bool & value);

  const std::string name;
  const Types type;
  const double lower_limit;
  const double upper_limit;

protected:
  bool bool_value{false};
  double double_value{0.0};
};

template <
  typename T, typename = std::enable_if_t<std::is_same_v<T, bool> || std::is_same_v<T, double>>>
T get_parameter(Parameter::SharedPtr const & parameter)
{
  switch (parameter->type) {
    case Types::BOOL:
      return parameter->get_bool_value();
    case Types::DOUBLE:
      return parameter->get_double_value();
    default:
      throw std::invalid_argument("Invalid parameter type");
  }
}

}  // namespace parameters
}  // namespace fmc

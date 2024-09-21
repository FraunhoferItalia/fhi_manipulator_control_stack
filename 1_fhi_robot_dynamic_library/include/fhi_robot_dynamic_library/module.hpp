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
#include <eigen3/Eigen/Dense>
#include <iostream>

#include "yaml-cpp/yaml.h"

namespace frdl
{

enum ModuleType { REVOLUTE, PRISMATIC, FIXED };

std::map<std::string, ModuleType> MODULE_TYPE_MAP = {
  {"revolute", ModuleType::REVOLUTE},
  {"prismatic", ModuleType::PRISMATIC},
  {"fixed", ModuleType::FIXED}};

struct Module
{
  std::string name;
  ModuleType type;
  double mass;
  Eigen::Vector3d CoM;
  Eigen::Matrix3d inertia_matrix;
  double motor_inertia;
  double static_friction;
  double dynamic_friction;
  double gear_ratio;
  Eigen::Vector3d axis;
  Eigen::Affine3d origin;
  Eigen::Vector2d position_limits;
  double speed_limit;
  double effort_limit;

  void print()
  {
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Name: " << name << std::endl;
    std::cout << "Type: " << type << std::endl;
    std::cout << "Mass: " << mass << std::endl;
    std::cout << "CoM: " << CoM << std::endl;
    std::cout << "Inertia Matrix: " << inertia_matrix << std::endl;
    std::cout << "Motor Inertia: " << motor_inertia << std::endl;
    std::cout << "Static Friction: " << static_friction << std::endl;
    std::cout << "Dynamic Friction: " << dynamic_friction << std::endl;
    std::cout << "Gear Ratio: " << gear_ratio << std::endl;
    std::cout << "Axis: " << axis << std::endl;
    std::cout << "Origin: " << origin.matrix() << std::endl;
    std::cout << "Position limits: " << position_limits << std::endl;
    std::cout << "Speed limit: " << speed_limit << std::endl;
    std::cout << "Effort limit: " << effort_limit << std::endl;
    std::cout << std::endl;
  }
};

Eigen::Affine3d vector_to_affine3d(const std::vector<double> & vec)
{
  if (vec.size() != 6) {
    throw std::invalid_argument("Invalid vector size. Expected size: 6.");
  }
  return Eigen::Translation3d(vec[0], vec[1], vec[2]) *
         Eigen::AngleAxisd(vec[5], Eigen::Vector3d::UnitZ()) *
         Eigen::AngleAxisd(vec[4], Eigen::Vector3d::UnitY()) *
         Eigen::AngleAxisd(vec[3], Eigen::Vector3d::UnitX());
}

class ChainModules
{
public:
  using SharedPtr = std::shared_ptr<ChainModules>;
  using UniquePtr = std::unique_ptr<ChainModules>;

  ChainModules(const std::string & modules_yaml_description)
  {
    YAML::Node config = YAML::Load(modules_yaml_description);

    if (!config.IsMap()) {
      std::cout << "Unable to load modules" << std::endl;
      return;
    }

    for (auto it = config["modules"].begin(); it != config["modules"].end(); ++it) {
      Module m;
      m.name = it->first.as<std::string>();
      for (auto it_par = it->second.begin(); it_par != it->second.end(); ++it_par) {
        if (it_par->first.as<std::string>() == "type") {
          m.type = MODULE_TYPE_MAP.at(it_par->second.as<std::string>());
        } else if (it_par->first.as<std::string>() == "mass") {
          m.mass = it_par->second.as<double>();
        } else if (it_par->first.as<std::string>() == "CoM") {
          m.CoM = Eigen::Map<Eigen::Vector3d>(it_par->second.as<std::vector<double>>().data());
        } else if (it_par->first.as<std::string>() == "inertia_matrix") {
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              m.inertia_matrix(i, j) = it_par->second[i][j].as<double>();
            }
          }
        } else if (it_par->first.as<std::string>() == "motor_inertia") {
          m.motor_inertia = it_par->second.as<double>();
        } else if (it_par->first.as<std::string>() == "static_friction") {
          m.static_friction = it_par->second.as<double>();
        } else if (it_par->first.as<std::string>() == "dynamic_friction") {
          m.dynamic_friction = it_par->second.as<double>();
        } else if (it_par->first.as<std::string>() == "gear_ratio") {
          m.gear_ratio = it_par->second.as<double>();
        } else if (it_par->first.as<std::string>() == "axis") {
          m.axis = Eigen::Map<Eigen::Vector3d>(it_par->second.as<std::vector<double>>().data());
        } else if (it_par->first.as<std::string>() == "origin") {
          m.origin = vector_to_affine3d(it_par->second.as<std::vector<double>>());
        } else if (it_par->first.as<std::string>() == "position_limits") {
          m.position_limits =
            Eigen::Map<Eigen::Vector2d>(it_par->second.as<std::vector<double>>().data());
        } else if (it_par->first.as<std::string>() == "speed_limit") {
          m.speed_limit = it_par->second.as<double>();
        } else if (it_par->first.as<std::string>() == "effort_limit") {
          m.effort_limit = it_par->second.as<double>();
        } else {
          std::cout << "Unrecognized parameter " << it_par->first.as<std::string>() << std::endl;
        }
      }
      modules.push_back(m);
      switch (m.type) {
        case ModuleType::PRISMATIC:
        case ModuleType::REVOLUTE:
          active_modules.push_back(m);
          break;
        case ModuleType::FIXED:
          fixed_modules.push_back(m);
          break;
        default:
          break;
      }
    }
  };
  size_t size() const { return active_modules.size(); }

  std::vector<Module> active_modules;
  std::vector<Module> fixed_modules;
  std::vector<Module> modules;
};

}  // namespace frdl
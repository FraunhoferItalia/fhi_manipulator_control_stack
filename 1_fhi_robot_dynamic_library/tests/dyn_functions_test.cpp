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
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <fhi_robot_dynamic_library/chain.hpp>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

void printProgressBar(int progress, int total, int barWidth = 50)
{
  float progressRatio = static_cast<float>(progress) / total;
  int filledWidth = static_cast<int>(progressRatio * barWidth);

  std::cout << "[";

  for (int i = 0; i < barWidth; ++i) {
    if (i < filledWidth)
      std::cout << "=";
    else
      std::cout << " ";
  }

  std::cout << "] " << int(progressRatio * 100.0) << "%\r";
  std::cout.flush();
}

double computeAverage(const std::vector<double> & numbers)
{
  double sum = 0.0;
  for (const auto & num : numbers) {
    sum += num;
  }
  return sum / numbers.size();
}

double computeStandardDeviation(const std::vector<double> & numbers)
{
  double mean = computeAverage(numbers);
  double squaredDifferencesSum = 0.0;
  for (const auto & num : numbers) {
    double difference = num - mean;
    squaredDifferencesSum += difference * difference;
  }
  double meanSquaredDifferences = squaredDifferencesSum / numbers.size();
  return std::sqrt(meanSquaredDifferences);
}

int random_index(int min, int max)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(min, max);
  return dis(gen);
}

// Read vector of Eigen::VectorXd from CSV
std::vector<Eigen::VectorXd> vectorXd_from_CSV(const std::string & filename)
{
  std::vector<Eigen::VectorXd> data;
  std::ifstream file(filename);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      Eigen::VectorXd row;
      std::stringstream ss(line);
      std::string value;
      while (std::getline(ss, value, ',')) {
        row.conservativeResize(row.size() + 1);
        row(row.size() - 1) = std::stod(value);
      }
      data.push_back(row);
    }
    file.close();
  } else {
    std::cerr << "Failed to open file " << filename << " for reading." << std::endl;
  }
  return data;
}

Eigen::Affine3d affine3d_from_vectorXd(const Eigen::VectorXd & vec)
{
  Eigen::Affine3d affine;
  affine.matrix() = Eigen::Map<const Eigen::Matrix<double, 4, 4, Eigen::RowMajor>>(vec.data());
  return affine;
}

Eigen::MatrixXd matrixXd_from_vectorXd(
  const Eigen::VectorXd & vec, const int & rows, const int & cols)
{
  Eigen::MatrixXd matrix = Eigen::Map<const Eigen::MatrixXd>(vec.data(), cols, rows).transpose();
  return matrix;
}

std::vector<Eigen::Affine3d> affine3d_from_vectorXd(const std::vector<Eigen::VectorXd> & vec)
{
  std::vector<Eigen::Affine3d> affine;
  for (auto & v : vec) {
    affine.push_back(affine3d_from_vectorXd(v));
  }
  return affine;
}

std::vector<Eigen::MatrixXd> square_matrixXd_from_vectorXd(const std::vector<Eigen::VectorXd> & vec)
{
  int size = vec[0].size();
  int rows = std::sqrt(size);
  assert(rows * rows == size && "Error: Vector size is not a perfect square.");
  std::vector<Eigen::MatrixXd> matrix;
  for (auto & v : vec) {
    matrix.push_back(matrixXd_from_vectorXd(v, rows, rows));
  }
  return matrix;
}

std::vector<Eigen::MatrixXd> matrixXd_from_vectorXd(
  const std::vector<Eigen::VectorXd> & vec, const int & rows, const int & cols)
{
  std::vector<Eigen::MatrixXd> matrix;
  for (auto & v : vec) {
    matrix.push_back(matrixXd_from_vectorXd(v, rows, cols));
  }
  return matrix;
}

frdl::ChainModules get_modules(std::string const & data_dir)
{
  std::ifstream fin(data_dir + "/chain_config.yaml");
  if (!fin) {
    std::cerr << "Failed to open file: " << data_dir << "/chain_config.yaml" << std::endl;
  }
  std::stringstream buffer;
  buffer << fin.rdbuf();
  std::string yamlContent = buffer.str();
  frdl::ChainModules chain_modules = frdl::ChainModules(yamlContent);
  for (auto m : chain_modules.modules) m.print();
  return chain_modules;
}

using namespace std::literals;

class TestResult
{
public:
  TestResult(const std::string & name, size_t test_size) : name_{name}, size{test_size} {};
  void print_stats()
  {
    std::cout << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "-----    " << name_ << " RESULT    ----- " << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    std::cout << "Passed percentage: " << equal_count / size * 100 << std::endl;
    auto maxElement = std::max_element(times.begin(), times.end());
    maxElement = std::min_element(times.begin(), times.end());

    double average = computeAverage(times);
    std::cout << "Execution time average: " << average << "microseconds" << std::endl;

    double standardDeviation = computeStandardDeviation(times);
    std::cout << "Execution time standard deviation: " << standardDeviation << "microseconds"
              << std::endl;

    std::cout << std::endl;
  }
  void start()
  {
    stopped_ = false;
    start_time = std::chrono::high_resolution_clock::now();
  }

  void stop()
  {
    if (!stopped_) {
      times.push_back(std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start_time)
                        .count());
      stopped_ = true;
    }
  }

  template <typename T>
  void validate(const T & result, const T & gt, double precision = 1e-6)
  {
    stop();
    if (result.isApprox(gt, precision)) {
      equal_count++;
    }
  }
  std::vector<double> times{};
  double equal_count{0.0};

private:
  std::string name_;
  size_t size;
  std::chrono::_V2::system_clock::time_point start_time;
  bool stopped_{false};
};

Eigen::MatrixXd get_dJ_numeric(
  const frdl::ChainModules & modules, const Eigen::VectorXd & q, const Eigen::VectorXd & dq,
  double dt = 0.001)
{
  Eigen::VectorXd q_next = q + dq * dt;
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, modules.size());
  Eigen::MatrixXd J_next = Eigen::MatrixXd::Zero(6, modules.size());

  frdl::Kinematics::get_jacobian(modules, q, J);
  frdl::Kinematics::get_jacobian(modules, q_next, J_next);
  return (J_next - J) / dt;
}

int main(int argc, char ** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <path_to_data_dir>" << std::endl;
    return 1;
  }
  std::string data_dir = argv[1];
  frdl::ChainModules modules = get_modules(data_dir);
  std::vector<Eigen::VectorXd> q = vectorXd_from_CSV(data_dir + "/input/q.csv");
  std::vector<Eigen::VectorXd> dq = vectorXd_from_CSV(data_dir + "/input/dq.csv");
  std::vector<Eigen::VectorXd> ddq = vectorXd_from_CSV(data_dir + "/input/ddq.csv");

  std::vector<Eigen::Affine3d> x_gt =
    affine3d_from_vectorXd(vectorXd_from_CSV(data_dir + "/gt/x.csv"));
  std::vector<Eigen::MatrixXd> J_gt =
    matrixXd_from_vectorXd(vectorXd_from_CSV(data_dir + "/gt/J.csv"), 6, modules.size());
  std::vector<Eigen::MatrixXd> dJ_gt =
    matrixXd_from_vectorXd(vectorXd_from_CSV(data_dir + "/gt/dJ.csv"), 6, modules.size());
  std::vector<Eigen::VectorXd> tau_gt = vectorXd_from_CSV(data_dir + "/gt/tau.csv");
  std::vector<Eigen::VectorXd> g_gt = vectorXd_from_CSV(data_dir + "/gt/g.csv");
  std::vector<Eigen::MatrixXd> M_gt =
    square_matrixXd_from_vectorXd(vectorXd_from_CSV(data_dir + "/gt/M.csv"));
  std::vector<Eigen::MatrixXd> C_gt =
    square_matrixXd_from_vectorXd(vectorXd_from_CSV(data_dir + "/gt/C.csv"));

  size_t size = q.size();
  std::cout << "Testing on " << size << " data" << std::endl;
  assert(size == dq.size());
  assert(size == ddq.size());
  assert(size == x_gt.size());
  assert(size == J_gt.size());
  assert(size == dJ_gt.size());
  assert(size == tau_gt.size());
  assert(size == g_gt.size());
  assert(size == M_gt.size());
  assert(size == C_gt.size());

  Eigen::Affine3d x = Eigen::Affine3d::Identity();
  Eigen::VectorXd ik = Eigen::VectorXd::Zero(modules.size());
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(modules.size());
  Eigen::VectorXd g = Eigen::VectorXd::Zero(modules.size());
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, modules.size());
  Eigen::MatrixXd dJ = Eigen::MatrixXd::Zero(6, modules.size());
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(modules.size(), modules.size());
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(modules.size(), modules.size());

  // size = 1000;
  TestResult fk_result("FK", size);
  TestResult ik_result("IK", size);
  TestResult J_result("JACOBIAN", size);
  TestResult dJ_result("dJACOBIAN", size);
  TestResult rne_result("RNE", size);
  TestResult g_result("gravity", size);
  TestResult M_result("M", size);
  TestResult C_result("C", size);
  for (size_t i = 0; i < size; i++) {
    printProgressBar(i, size);
    fk_result.start();
    frdl::Kinematics::get_forward_kinematic(modules, q[i], x);
    fk_result.validate(x, x_gt[i]);
    // std::cout << "x" << std::endl;
    // std::cout << x.matrix() << std::endl;
    // std::cout << "x_gt" << std::endl;
    // std::cout << x_gt[i].matrix() << std::endl;

    // frdl::Kinematics::IKParameters ik_parameters(q[random_index(0, size - 1)]);
    // ik_result.start();
    // if (frdl::Kinematics::get_inverse_kinematic(modules, x_gt[i], ik_parameters, ik)) {
    //   frdl::Kinematics::get_forward_kinematic(modules, ik, x);
    //   ik_result.validate(x, x_gt[i]);
    // }

    // J_result.start();
    // frdl::Kinematics::get_jacobian(modules, q[i], J);
    // J_result.validate(J, J_gt[i]);
    // std::cout << "J" << std::endl;
    // std::cout << J << std::endl;
    // std::cout << "J_gt" << std::endl;
    // std::cout << J_gt[i] << std::endl;

    dJ_result.start();
    J_result.start();
    frdl::Kinematics::get_jacobian_dot(modules, q[i], dq[i], J, dJ);
    J_result.stop();
    dJ_result.stop();
    J_result.validate(J, J_gt[i]);
    dJ_gt[i] = get_dJ_numeric(modules, q[i], dq[i]);
    dJ_result.validate(dJ, dJ_gt[i], 1e-2);
    // std::cout << "dJ" << std::endl;
    // std::cout << dJ << std::endl;
    // std::cout << "dJ_gt" << std::endl;
    // std::cout << dJ_gt[i] << std::endl;

    rne_result.start();
    frdl::Dynamics::get_rne(modules, q[i], dq[i], dq[i], ddq[i], tau);
    // std::cout << "tau" << std::endl;
    // std::cout << tau << std::endl;
    // std::cout << "tau_gt" << std::endl;
    // std::cout << tau_gt[i] << std::endl;
    rne_result.validate(tau, tau_gt[i], 1e-3);

    g_result.start();
    frdl::Dynamics::get_torque_gravity(modules, q[i], g);
    // std::cout << "g" << std::endl;
    // std::cout << g << std::endl;
    // std::cout << "g_gt" << std::endl;
    // std::cout << g_gt[i] << std::endl;
    g_result.validate(g, g_gt[i], 1e-3);

    M_result.start();
    frdl::Dynamics::get_M(modules, q[i], M);
    // std::cout << "M" << std::endl;
    // std::cout << M.matrix() << std::endl;
    // std::cout << "M_gt" << std::endl;
    // std::cout << M_gt[i].matrix() << std::endl;
    M_result.validate(M, M_gt[i], 1e-3);

    C_result.start();
    frdl::Dynamics::get_C(modules, q[i], dq[i], C);
    C_result.stop();
    frdl::Dynamics::get_rne(
      modules, q[i], dq[i], dq[i], Eigen::VectorXd::Zero(modules.size()), tau,
      frdl::Dynamics::RNEOptions::zero_gravity());
    Eigen::VectorXd tauC = C * dq[i];
    C_result.validate(tauC, tau, 1e-2);
  }

  fk_result.print_stats();
  ik_result.print_stats();
  J_result.print_stats();
  dJ_result.print_stats();
  rne_result.print_stats();
  g_result.print_stats();
  M_result.print_stats();
  C_result.print_stats();
  return 0;
}

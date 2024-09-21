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

// #include <atomic> // used for double thread

// Include Files
#include "rtwtypes.h"
#include <atomic> // used for multiple threads
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <vector>
#include <iostream>

// Include Files for parameter tuning and robust controller
#include "matlab_interface.h"
#include "coder_bounded_array.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

// ----------------- Print Function for debugging --------------

inline void print_coder(coder::array<double, 1U> &var)
{
    for (int i = 0; i < var.size(0); i++)
    {
        std::cout << var.at(i) << " ";
    }
    std::cout << "\n";
}

inline void print_coder(coder::array<double, 2U> &var)
{

    for (int i = 0; i < var.size(0); i++)
    {
        for (int j = 0; j < var.size(1); j++)
        {
            std::cout << var.at(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

inline void print_coder(coder::array<double, 3U> &var)
{

    for (int k = 0; k < var.size(2); k++)
    {
        std::cout << "--- " << (k + 1) << " matrix ---\n";
        for (int i = 0; i < var.size(0); i++)
        {
            for (int j = 0; j < var.size(1); j++)
            {
                std::cout << var.at(i, j, k) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

// ----------------- Json interface functions ------------------

template <typename T>
bool read_scalar(T &var, json const &d, std::string const &name)
{
    if (std::is_same<T, double>::value || std::is_same<T, int>::value)
    {
        if (d.find(name) != d.end())
        {
            var = d.at(name);
            std::cout << name << "=" << var << "\n";
        }
        else
        {
            std::cout << "\nkey " << name << " not found \n";
            return false;
        }
        return true;
    }
    else
    {
        return false;
    }
}

inline void read_double_array(double *var, json d, const char *name, int imax, int jmax = 0, int kmax = 0)
{
    int err = 0;

    if (d.find(name) == d.end())
    {
        // std::cout << "\nkey " << name << " not found \n";
        return;
    }

    for (int i = 0; i < imax; i++)
    {
        if (jmax == 0)
        {
            var[i] = d.at(name)[i];
            //          var[i]=d[name][i];
        }
        else
        {
            if (kmax == 0)
            {
                for (int j = 0; j < jmax; j++)
                {
                    var[i * jmax + j] = d.at(name)[i][j];
                    //                  std::cout << name << " = " << d.at(name)[i][j];
                    //                  var[i*jmax+j]=d[name][i][j];
                }
            }
            else
            {
                for (int j = 0; j < jmax; j++)
                {
                    for (int k = 0; k < kmax; k++)
                    {
                        var[i * jmax * kmax + j * kmax + k] = d.at(name)[i][j][k];
                    }
                }
            }
        }
    }
}

inline bool read_vect(coder::array<double, 1U> &var, json const &d, const std::string &name, int imax) // to coder::array
{
    if (d.find(name) == d.end())
    {
        std::cout << "\nkey " << name << " not found \n";
        return false;
    }

    var.set_size(imax);
    if (imax == 1)
    {
        var[0] = d.at(name);
    }
    else
    {
        for (int i = 0; i < imax; i++)
        {
            var[i] = d.at(name)[i];
        }
    }
    return true;
}

inline bool read_vect(coder::array<double, 2U> &var, json const &d, const std::string &name, int imax, int jmax) // to coder::array
{
    if (d.find(name) == d.end())
    {
        std::cout << "\nkey " << name << " not found \n";
        return false;
    }

    var.set_size(imax, jmax);

    if (imax == 1)
    {
        for (int j = 0; j < jmax; j++)
        {
            var.at(0, j) = d.at(name)[j];
        }
    }
    else
    {
        for (int i = 0; i < imax; i++)
        {
            for (int j = 0; j < jmax; j++)
            {
                var.at(i, j) = d.at(name)[i][j];
            }
        }
    }

    return true;
}

inline bool read_vect(coder::array<double, 3U> &var, json const &d, const std::string &name, int imax, int jmax, int kmax) // to coder::array
{
    if (d.find(name) == d.end())
    {
        std::cout << "\nkey " << name << " not found \n";
        return false;
    }

    var.set_size(imax, jmax, kmax);

    for (int i = 0; i < imax; i++)
    {
        for (int j = 0; j < jmax; j++)
        {
            for (int k = 0; k < kmax; k++)
            {
                var.at(i, j, k) = d.at(name)[i][j][k];
            }
        }
    }

    return true;
}

inline bool read_struct(coder::array<double, 1U> &var, json const &d, const char *key1, const char *key2, int imax)
{
    json d_sub;

    d_sub = d[key1][key2];

    var.set_size(imax);

    if (imax == 1)
    {
        var[0] = d_sub;
    }
    else
    {
        for (int i = 0; i < imax; i++)
        {
            var[i] = d_sub[i];
        }
    }

    return true;
}

inline bool read_struct(coder::array<double, 2U> &var, json const &d, const char *key1, const char *key2, int imax, int jmax)
{
    json d_sub;

    d_sub = d[key1][key2];

    var.set_size(imax, jmax);

    if (imax == 1)
    {
        for (int j = 0; j < jmax; j++)
        {
            var.at(0, j) = d_sub[j];
        }
    }
    else
    {
        for (int i = 0; i < imax; i++)
        {
            for (int j = 0; j < jmax; j++)
            {
                var.at(i, j) = d_sub[i][j];
            }
        }
    }

    return true;
}

inline bool read_struct(coder::array<double, 3U> &var, json const &d, const char *key1, const char *key2, int imax, int jmax, int kmax)
{
    json d_sub;

    d_sub = d[key1][key2];

    var.set_size(imax, jmax, kmax);

    for (int i = 0; i < imax; i++)
    {
        for (int j = 0; j < jmax; j++)
        {
            for (int k = 0; k < kmax; k++)
            {
                var.at(i, j, k) = d_sub[i][j][k];
            }
        }
    }

    return true;
}

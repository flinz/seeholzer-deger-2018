/*
 *  mymodule.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef MYMODULE_H
#define MYMODULE_H

#include "slimodule.h"
#include "slifunction.h"
#include "topologymodule.h"
#include "parameter.h"
#include "randomgen.h"
#include "normal_randomdev.h"

// Put your stuff into your own namespace.
namespace mynest
{

/**
 * Class defining your model.
 * @note For each model, you must define one such class, with a unique name.
 */
class MyModule : public SLIModule
{
public:
  // Interface functions ------------------------------------------

  /**
   * @note The constructor registers the module with the dynamic loader.
   *       Initialization proper is performed by the init() method.
   */
  MyModule();

  /**
   * @note The destructor does not do much in modules.
   */
  ~MyModule();

  /**
   * Initialize module by registering models with the network.
   * @param SLIInterpreter* SLI interpreter
   */
  void init( SLIInterpreter* );

  /**
   * Return the name of your model.
   */
  const std::string name( void ) const;

  /**
   * Return the name of a sli file to execute when mymodule is loaded.
   * This mechanism can be used to define SLI commands associated with your
   * module, in particular, set up type tries for functions you have defined.
   */
  const std::string commandstring( void ) const;

public:
  // Classes implementing your functions -----------------------------
  
  /**
   * Defines a Gaussian weight distribution with additive (normally 
   * distributed) noise on weights controlled by sigma_noise.
   */
  class GaussianNoisyParameter : public nest::Parameter
  {
    public:
      /**
       * Parameters:
       * c        - constant offset
       * p_center - value at center of gaussian
       * mean     - distance to center
       * sigma    - width of gaussian
       * sigma_noise  - standard deviation of additive normal noise 
       */
      GaussianNoisyParameter(const DictionaryDatum& d):
        Parameter(d),
        c_(0.0),
        p_center_(1.0),
        mean_(0.0),
        sigma_(1.0),
        sigma_noise_(0.0),
        rdev()
        {
          updateValue<double_t>(d, nest::names::c, c_);
          updateValue<double_t>(d, nest::names::p_center, p_center_);
          updateValue<double_t>(d, nest::names::mean, mean_);
          updateValue<double_t>(d, nest::names::sigma, sigma_);
          updateValue<double_t>(d, "sigma_noise", sigma_noise_);
        }


      double_t raw_value(const nest::Position<2> &p, librandom::RngPtr& rng) const
        {
          return std::max(
              c_ + p_center_*
            std::exp(-std::pow(p.length() - mean_,2)/(2*std::pow(sigma_,2)))
              + rdev(rng) * sigma_noise_
            ,0.);
        }

      Parameter * clone() const
        { return new GaussianNoisyParameter(*this); }

    private:
      double_t c_, p_center_, mean_, sigma_, sigma_noise_;
      librandom::NormalRandomDev rdev;

  };
};
} // namespace mynest

#endif

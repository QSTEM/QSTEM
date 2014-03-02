/*
  QSTEM - image simulation for TEM/STEM/CBED
  Copyright (C) 2000-2010  Christoph Koch
  Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RANDOM_H
#define RANDOM_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>

namespace QSTEM
{

static boost::random::mt11213b _rng;         // produces randomness out of thin air.  This is the algorithm.  To get numbers,
//  use ran1(_rng) for uniform values, or 
static boost::random::uniform_01<float_tt> _ran1;    // This returns uniform random values between 0 and 1
inline float_tt ran1(){return _ran1(_rng);}
static boost::random::normal_distribution<float_tt> _gasdev(0,1); // This returns uniform random values with a Gaussian normal distribution
inline float_tt gasdev(){return _gasdev(_rng);}

} // end namespace QSTEM

#endif

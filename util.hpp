/*
 * util.hpp
 *
 *  Created on: Sep 5, 2022
 *      Author: valeh
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include "params.hpp"

using P = params::poly_t;
using V = typename P::value_type;

P sample_c(uint8_t mode, uint8_t* hash);

void sample_alpha(P *out, uint8_t* hash, int j);

#endif /* UTIL_HPP_ */

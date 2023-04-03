/*
 * params.hpp
 *
 *  Created on: Sep 5, 2022
 *      Author: valeh
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include "nfl.hpp"

namespace params{
	using poly_t = nfl::poly_from_modulus<uint64_t, 4096, nfl::params<uint64_t>::kModulusBitsize>;
	poly_t::value_type plaintextModulus = 10000009;  //
	long plaintextModulusBitsize = 24;
	const size_t MU=1, LAMBDA=1;
}

#define NV 10   // number of voters
#define NC 2   // number of candidates
#define LC 1   // maximum vote that can be given to a single candidate
#define LS 1   // total amount of votes a single voter can give - vote budget
#define delta_1 (1L<<32)  // rejection sampling parameter
#define beta_1 4096        // rejection sampling and MSIS parameter
#define k 2            // soundness boost and automorphism order parameter

#define SYMBYTES 32       // output of hash functions in bytes

#endif /* PARAMS_HPP_ */

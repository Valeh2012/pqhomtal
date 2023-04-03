#pragma once

#include <cstddef>

#include <algorithm>
#include <iostream>
#include <nfl.hpp>


namespace BDLOP{
template <std::size_t KAPPA>
class comm_rnd_t;
template <std::size_t KAPPA>
class comm_pk_t;
template <std::size_t KAPPA>
class comm_t;
}

namespace BDLOP{
template < std::size_t KAPPA>
class comm_rnd_t{
	using P = params::poly_t;

public:
	std::array<P, params::MU> rs;
	std::array<P, KAPPA> re;
	std::array<P, params::LAMBDA> rem;


	comm_rnd_t(){

		for(std::size_t i=0;i<params::MU;i++){
			rs[i] = nfl::non_uniform(2);
			rs[i].ntt_pow_phi();
		}

		for(std::size_t i=0;i<KAPPA;i++){
			re[i] = nfl::non_uniform(2);
			re[i].ntt_pow_phi();
		}

		for(std::size_t i=0;i<params::LAMBDA;i++){
			rem[i] = nfl::non_uniform(2);
			rem[i].ntt_pow_phi();
		}
	}

//	~comm_rnd_t(){}
};
}

namespace BDLOP{
template <std::size_t KAPPA>
class comm_pk_t{
	using P = params::poly_t;

public:
	std::array<std::array<P, params::LAMBDA>, params::MU> B0l;
	std::array<std::array<P, KAPPA>, params::MU> B0k;
	std::array<std::array<P, params::LAMBDA>, KAPPA> b;


	comm_pk_t(){

		for(std::size_t i=0;i<params::MU;i++){
			for(std::size_t j=0;j<params::LAMBDA;j++){
				B0l[i][j] = nfl::uniform();
			}
		}

		for(std::size_t i=0;i<params::MU;i++){
			for(std::size_t j=0;j<KAPPA;j++){
				B0k[i][j] = nfl::uniform();
			}
		}
		for(std::size_t i=0;i<KAPPA;i++){
			for(std::size_t j=0;j<params::LAMBDA;j++){
				b[i][j] = nfl::uniform();
			}
		}
	}
//	~comm_pk_t(){}
};


}

namespace BDLOP{
template <std::size_t KAPPA>
class comm_t{
	using P = params::poly_t;
public:

	std::array<P, params::MU> t0;
	std::array<P, KAPPA> t;

	// constructor
	comm_t(){};

	template <class RND, class PK>
	comm_t(const PK& pk, const RND& randomness){

		for(std::size_t i=0; i<params::MU;i++){
			t0[i] = std::inner_product(pk.B0l[i].begin(), pk.B0l[i].end(), randomness.rem.begin(), P{0});
		}
		std::array<P, params::MU> tag;
		for(std::size_t i=0; i<params::MU;i++){
			tag[i] = std::inner_product(pk.B0k[i].begin(), pk.B0k[i].end(), randomness.re.begin(), P{0});
		}
		for(std::size_t i=0; i<KAPPA;i++){
			t[i] = std::inner_product(pk.b[i].begin(), pk.b[i].end(), randomness.rem.begin(), P{0});
		}

		for(std::size_t i=0; i<params::MU; i++){
			t0[i] = t0[i] + tag[i] +randomness.rs[i];
		}

		for(std::size_t i=0; i<KAPPA;i++){
			t[i] = t[i] + randomness.re[i];
		}

	}
//	~comm_t(){}

	// Assign operators
};
}

namespace BDLOP{
template <class COMM>
void commit_polyvec(COMM &comm, std::vector<params::poly_t> &polyvec, std::size_t pos, bool isNTT){

	std::size_t MSGLEN = polyvec.size();
	assert(pos + MSGLEN <= comm.t.size());

	if(!isNTT){
		for(auto &v : polyvec)
			v.ntt_pow_phi();
	}

	for(std::size_t i=0;i<MSGLEN; i++){
		comm.t[pos+i] = comm.t[pos+i] + polyvec[i];
	}
}

template <class COMM>
void commit_poly(COMM &comm, params::poly_t &poly, std::size_t pos, bool isNTT){

	assert(pos < comm.t.size());

	if(!isNTT){
		poly.ntt_pow_phi();
	}

	comm.t[pos] = comm.t[pos] + poly;
}
}

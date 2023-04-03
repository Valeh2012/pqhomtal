#include "util.hpp"

extern "C" {
	#include "fips202.h"
}

P sample_ternary(uint8_t mode, uint8_t* hash){
	P p;
	const V pm = P::get_modulus(0) - 1u;
	uint8_t rnd[P::degree];
	V tmp[P::degree];
	if(hash == nullptr) nfl::fastrandombytes(rnd, P::degree); // TODO: remove randomness for null pointer
	else shake256(rnd,sizeof(rnd),hash,32);  // TODO: remove XOF after hashing
	for(std::size_t i=0;i<P::degree;i++){
		tmp[i] = rnd[i] <= mode ? pm + (rnd[i] & 2)  : 0u;
	}
	p.set(tmp, tmp+P::degree, true);
	return p;
}

void sample_uniform(P *out, uint8_t* hash, int j){

	size_t outlen = j*sizeof(P);
	uint8_t *rnd = (uint8_t*) malloc(outlen);
//	uint8_t rnd[(j+1)*param_N*2];
	memset(rnd, 0, outlen); //TODO: change to incremental api.
	shake256(rnd,outlen,hash, 32);
	for(int i=0;i<j;i++){
		V * ptr = (V*) (rnd+i*sizeof(P));
		out[i].set(ptr, ptr+ P::degree, true);
	}
	free(rnd);
}

#include <iostream>
#include <string>
#include <memory>
#include <chrono>
#include <vector>
#include <nfl.hpp>

using namespace std;

#include "params.hpp"
#include <BGV.hpp>


typedef BGV::ciphertext_t ctx_t;
typedef BGV::pk_t pk_t;
typedef BGV::sk_t sk_t;

#define NV 100000
#define NC 100
#define LC 1


int main(){

	srand(0);

	sk_t secret_key;
	pk_t public_key(secret_key);

	using P = params::poly_t;

	ctx_t tally;
	P zero = P{0};
	ctx_t tmp;
	for(int j=0;j<NC;j++){
		int c = rand() % (LC+1);
		zero.data()[j] = c;
	}
	zero.ntt_pow_phi();
	BGV::encrypt_poly(tmp, public_key, zero);

	auto start = chrono::steady_clock::now();
	for(int i=0;i<NV;i++){
		tally += tmp;
	}
	auto end = chrono::steady_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();


	P result;
	BGV::decrypt_poly(result, secret_key, public_key, tally);
	cout << duration <<"ms\n";

	//cout << result << endl;
//	P::value_type sum = 0;
//	for(auto &v : result){
//		sum += v;
//	}
//
//	cout << sum << endl;
	return 0;
}

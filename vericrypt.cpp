#include <iostream>
#include <string>
#include <memory>
#include <random>
#include <chrono>
#include <vector>
#include <numeric>


extern "C" {
#include "fips202.h"
}

using namespace std;


template <class T, size_t Align, class... Args>
T* alloc_aligned(size_t n, Args&& ... args)
{
	T* ret;
	if (posix_memalign((void**) &ret, Align, sizeof(T)*n) != 0) {
		throw std::bad_alloc();
	}
	for (size_t i = 0; i < n; i++) {
		new (&ret[i]) T(std::forward<Args>(args)...);
	}
	return ret;
}


#include "params.hpp"
#include <BGV.hpp>
#include <BDLOP.hpp>

typedef BGV::ciphertext_t ctx_t;
typedef BGV::pk_t pk_t;
typedef BGV::sk_t sk_t;



#define kappa (LC+8)

#define nmod(a, m) ((a) < 0 ? (m - (-1*(a))%m) : (a)%m)


using P = params::poly_t;
using V = typename P::value_type;

array<V, NC> ballot_generator(){
	bool flag = false;
	int lim;
	array<V, NC> tmp;

	while(!flag){
		fill(tmp.begin(), tmp.end(), 0);
		lim = LS;
		for(int i=0; i<NC;i++){
			if(lim >0){
				auto b = rand() % (LC+1);
				tmp[i] = (V) b;
				lim -= b;
			}
		}
		auto acc = accumulate(tmp.begin(), tmp.end(), 0);
		flag = (acc == (V) LS);
	}
	return tmp;
}

void random_ballot(P& out, array<V, NC> ballot){
	shuffle(ballot.begin(), ballot.end(), default_random_engine(chrono::system_clock::now().time_since_epoch().count()));
	memcpy(out.begin(), ballot.begin(), NC*sizeof(V));
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

void sample_uniform_scalar_ntt(P *out, uint8_t* hash,int j){
	size_t outlen = j*sizeof(V);
	uint8_t *rnd = (uint8_t*) malloc(outlen);
//	uint8_t rnd[(j+1)*param_N*2];
	memset(rnd, 0, outlen); //TODO: change to incremental api.
	shake256(rnd,outlen,hash, 32);
	for(int i=0;i<j;i++){
		V* ptr = (V*) (rnd + i*sizeof(V));
		fill(out[i].begin(), out[i].end(), *ptr);
	}
	free(rnd);
}

void poly_transpose(P &in){
	V tmp;
	V q = P::get_modulus(0);
	auto l = P::degree;
	in.data()[l/2] = q - in.data()[l/2];
	for(size_t i=1;i< l / 2;i++){
		tmp = in.data()[i];
		in.data()[i] = q - in.data()[l-i];
		in.data()[l-i] = q - tmp;
	}
}


void poly_inner_product_Zq(mpz_t &res, P &a, P &b){
	mpz_t tmp1, tmp2;
	mpz_init(tmp1);
	mpz_init(tmp2);
	for(size_t i=0;i<P::degree;i++){
		mpz_set_ui(tmp1, a.data()[i]);
		mpz_set_ui(tmp2, b.data()[i]);
		mpz_addmul(res, tmp1, tmp2);
	}
}


void galois_transform(P &out, P &p, size_t i){
	// re-think about implementation

	// check gcd(i,2*degree) == 1

	size_t l  = P::degree;

	if( gcd(l, i) != 1 ){
		cerr << "no such transformation" << endl;
		exit(-1);
		return;
	}

	V modulus = P::get_modulus(0);
	P res{0}, tmp;
	int q,r;
	for(size_t j=0;j<l;j++){
		if(p.data()[j] == 0) continue;
		tmp.set(0);
		q = j*i / l;
		r = (j*i)%l;
		if(q & 1) tmp.data()[r] = modulus - p.data()[j];
		else tmp.data()[r] = p.data()[j];
		res = res + tmp;
	}


	out = res;
}

template<typename _Mn>
_Mn modInverse(_Mn a, _Mn m){
	_Mn m0 = m;
	_Mn y = 0, x = 1;

	if (m == 1)
		return 0;

	while (a > 1)
	{
		// q is quotient
		_Mn q = a / m;
		_Mn t = m;

		// m is remainder now, process same as
		// Euclid's algo
		m = a % m, a = t;
		t = y;

		// Update y and x
		y = x - q * y;
		x = t;
	}

	// Make x positive
	if (x < 0)
		x += m0;

	return x;
}

void inv_galois_transform(P& out, P &p, int i){
	// re-think about implementation

	int l = P::degree;

	if( !gcd(l, i) ){
		cerr << "no such transformation" << endl;
		exit(-1);
		return;
	}

//	cout << i << endl;
	// find multiplicative inverse mod 2l
	i = modInverse(i, 2*l);
//	cout << i << endl;
//	poly_type res;
	p.invntt_pow_invphi();
	galois_transform(out,p,i);
	out.ntt_pow_phi();
//	memcpy(&out, &res, sizeof(out));
}


V _abs(V a){
	P::signed_value_type tmp1;
	BGV::util::center(tmp1, a, P::get_modulus(0), P::get_modulus(0) >> 1);
	return (V) abs(tmp1);
}

template<std::size_t LEN>
V inf_norm(std::array<P, LEN> &polyvec){
	V max_el = 0;
	V tmp = 0;
	for(auto &it : polyvec){
		for(auto &el : it){
			tmp =  _abs(el);
			max_el = max(max_el, tmp);
		}
	}

	return max_el;
}


template<size_t KAPPA>
struct Proof{
	BDLOP::comm_t<KAPPA> t;
	array<BDLOP::comm_rnd_t<KAPPA>, k > z;
	uint8_t chash[SYMBYTES];
	ctx_t ct;
	P h;
};


using proof_t = Proof<kappa>;
typedef BDLOP::comm_pk_t<kappa> comm_key_t;

void vericrypt(proof_t &p, const pk_t &pubkey,const comm_key_t& commkey, const P &message){

	uint8_t symbuf[4*SYMBYTES];
	memset(symbuf, 0, 4*SYMBYTES);
	uint8_t *thash = symbuf; // to store hash of public variables before while loop
	uint8_t *whash = symbuf+2*SYMBYTES; // to store hash of variables withing while loop i.e depends on w
	uint8_t *chash = symbuf+3*SYMBYTES;

	shake128incctx state;

	P poly_m = message;
	poly_m.ntt_pow_phi();

	P t = {params::plaintextModulus};
	t.ntt_pow_phi();

	P r_enc = nfl::non_uniform(2);
	P e_u = nfl::non_uniform(2);
	P e_v = nfl::non_uniform(2);

	r_enc.ntt_pow_phi();
	p.ct.r = r_enc;

	p.ct.c0 = e_u;
	p.ct.c0.ntt_pow_phi();
	p.ct.e1 = p.ct.c0;
	p.ct.c0 = p.ct.c0 * t + r_enc * pubkey.a;

	p.ct.c1 = e_v;
	p.ct.c1.ntt_pow_phi();
	p.ct.e2 = p.ct.c1;
	p.ct.c1 = p.ct.c1 * t + r_enc * pubkey.b + poly_m;

	p.ct.isnull = false;

	vector<P> shat;
	shat.resize(0);
	r_enc.invntt_pow_invphi();
	shat.push_back(r_enc);
	shat.push_back(e_u);
	shat.push_back(e_v);
	shat.push_back(message);

	P tmp = message;
	P tmp2{0}, tmp3{0};
	for(int i=1; i<= LC;i++){
		fill(tmp2.begin(), tmp2.end(), (V) i);
		tmp2 = message - tmp2;
		tmp = tmp * tmp2;
		shat.push_back(tmp);
	}

	BDLOP::comm_rnd_t<kappa> r_comm;
	p.t = BDLOP::comm_t<kappa>(commkey, r_comm);

	BDLOP::commit_polyvec(p.t, shat, 0, true);

	P g = nfl::uniform();
	fill(g.begin(), g.begin() + k, 0);
	g.ntt_pow_phi();

	BDLOP::commit_poly(p.t, g, shat.size(), true);

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, (uint8_t *)p.t.t0.begin(), params::MU*sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)p.t.t.begin(),  (shat.size() + 1)*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(thash,SYMBYTES, &state);

	array<BDLOP::comm_rnd_t<kappa>, k> y;

	// start rejection sampling
	P* alphas = alloc_aligned<P, 32>(4*k, P{0});
	P* alphas_prime = alloc_aligned<P, 32>(LC*k, P{0});
	P* gammas = alloc_aligned<P, 32>(2*k, P{0});
	P* gammas_prime = alloc_aligned<P, 32>(k, P{0});

	bool rej = false;
	//auto counter = 1;

	p.ct.c0.invntt_pow_invphi();
	p.ct.c1.invntt_pow_invphi();

	P t5,t6,t7;

	while(!rej){
		//cout << "trial " << counter++ << endl;
		for(auto &el : y){
			for(auto &yi : el.re){
				yi = nfl::non_uniform(delta_1);
				yi.ntt_pow_phi();
			}
			for(auto &yi : el.rs){
				yi = nfl::non_uniform(delta_1);
				yi.ntt_pow_phi();
			}
			for(auto &yi : el.rem){
				yi = nfl::non_uniform(delta_1);
				yi.ntt_pow_phi();
			}
		}


		array<array<P, params::MU>, k> w;
		for(int i=0; i<k;i++){
			for(size_t j =0; j<params::MU;j++){
				w[i][j] = inner_product(commkey.B0l[j].begin(), commkey.B0l[j].end(), y[i].rem.begin(), P{0});
				tmp = inner_product(commkey.B0k[j].begin(), commkey.B0k[j].end(), y[i].re.begin(), P{0});
				w[i][j] = w[i][j] + tmp + y[i].rs[j];
			}
		}

		shake128_inc_init(&state);
		shake128_inc_absorb(&state, thash,  SYMBYTES);
		shake128_inc_absorb(&state, (uint8_t *)w.begin(), k*params::MU*sizeof(P));
		shake128_inc_finalize(&state);
		shake128_inc_squeeze(whash,2*SYMBYTES, &state);

		P one{1};fill(one.begin(), one.end(), (P::value_type) 1);

		sample_uniform(alphas, whash, 4*k);
		sample_uniform(alphas_prime, chash, LC*k);

		shake128_inc_init(&state);
		shake128_inc_absorb(&state, thash,  SYMBYTES);
		shake128_inc_absorb(&state, (uint8_t *)alphas, 4*k*sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *)alphas_prime, LC*k*sizeof(P));
		shake128_inc_finalize(&state);
		shake128_inc_squeeze(whash,2*SYMBYTES, &state);

		sample_uniform(gammas,whash, 2*k);
		sample_uniform_scalar_ntt(gammas_prime, chash, k);


		tmp = inner_product(commkey.b[LC+6].begin(), commkey.b[LC+6].end(), y[0].rem.begin(), P{0});
		tmp = tmp + y[0].re[LC+6];

		P v = inner_product(commkey.b[LC+5].begin(), commkey.b[LC+5].end(), y[0].rem.begin(), P{0});
		v = v + y[0].re[LC+5];

		P three, bjyi;
		fill(three.begin(), three.end(), (P::value_type) 3);
		tmp3 = {0};
		for(int i=0;i<k;i++){
			for(size_t j=0; j<3; j++){
				// <b_j, y_i>
				bjyi = inner_product(commkey.b[j].begin(), commkey.b[j].end(), y[i].rem.begin(), P{0});
				bjyi = bjyi + y[i].re[j];

				// m_nd+1 = <b_nd+2, y_0> - sum_ij alphas_ij * sigma^(-i)((3s_j)<b_j, y_i>^2))
				tmp2 = bjyi * bjyi;
				tmp2 = three * shat[j] * tmp2;
				for(int ii=0; ii<i; ii++) inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas[3*i+j] * tmp2;
				tmp  = tmp - tmp2;

				// m_nd+2 = sum_ij alphas_ij * sigma^(-i)((3s_j^2 - 1)<b_j, y_i>))
				tmp2 = three * shat[j] * shat[j];
				tmp2 = tmp2 - one;
				tmp2 = tmp2 * bjyi;
				for(int ii=0; ii<i; ii++) inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas[3*i+j] * tmp2;
				tmp3  = tmp3 +tmp2;

				// v = <b_nd+2, y_0> + sum_ij alphas_ij * sigma^(-i)(<b_j, y_i>^3))
				tmp2 = bjyi * bjyi * bjyi;
				for(int ii=0; ii<i; ii++) inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas[3*i+j] * tmp2;
				v  = v + tmp2;
			}
		}

		t5 = p.t.t[LC+5] + tmp; //BDLOP::commit_poly(p.t, tmp, LC+5, true);
		t6 = p.t.t[LC+6] + tmp3; //BDLOP::commit_poly(p.t, tmp3, LC+6, true);

		P v2 = inner_product(commkey.b[LC+7].begin(), commkey.b[LC+7].end(), y[0].rem.begin(), P{0});
		v2 = v2 + y[0].re[LC+7];

		tmp = {0};
		for(int i=0;i<k;i++){
			P b3yi = inner_product(commkey.b[3].begin(), commkey.b[3].end(), y[i].rem.begin(), P{0});
			b3yi = b3yi + y[i].re[3];
			for(size_t j=0;j<LC;j++){
				tmp2 = inner_product(commkey.b[3+j].begin(), commkey.b[3+j].end(), y[i].rem.begin(), P{0});
				tmp2 = tmp2 + y[i].re[3+j];

				tmp3 = b3yi * tmp2;
				for(int ii=0; ii<i; ii++) inv_galois_transform(tmp3, tmp3, 2*P::degree/k + 1);
				tmp3 = tmp3 * alphas_prime[i*LC+j];
				v2 = v2 + tmp3;

				fill(tmp3.begin(), tmp3.end(), (P::value_type) (j+1));
				tmp3 = shat[3] - tmp3;
				tmp2 = tmp2 * tmp3;
				tmp2 = tmp2 + b3yi * shat[3+j];
				tmp3 = inner_product(commkey.b[3+j+1].begin(), commkey.b[3+j+1].end(), y[i].rem.begin(), P{0});
				tmp3 = tmp3 + y[i].re[3+j+1];
				tmp2 = tmp3 - tmp2;
				for(int ii=0; ii<i; ii++) inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas_prime[i*LC+j] * tmp2;
				tmp  = tmp + tmp2;
			}
		}
		t7 = p.t.t[LC+7] + tmp; //BDLOP::commit_poly(p.t, tmp, LC+7, true);

		array<P, 4> pis;
		P pka = pubkey.a;
		pka.invntt_pow_invphi();
		poly_transpose(pka);
		pka.ntt_pow_phi();

		P pkb = pubkey.b;
		pkb.invntt_pow_invphi();
		poly_transpose(pkb);
		pkb.ntt_pow_phi();

		V kinv = (V) modInverse((P::signed_value_type) k, (P::signed_value_type) P::get_modulus(0)); // k * kinv = 1 mod q

		array<P, k> vs;
		P vprime;
		for(size_t i=0;i<k;i++){
			vs[i] = inner_product(commkey.b[LC+4].begin(), commkey.b[LC+4].end(), y[i].rem.begin(), P{0});
			vs[i] = vs[i] + y[i].re[LC+4];


			tmp = inner_product(commkey.b[LC+3].begin(), commkey.b[LC+3].end(), y[i].rem.begin(), P{0});
			tmp = tmp + y[i].re[LC+3];

			for(size_t ii=0;ii<i;ii++) inv_galois_transform(tmp, tmp, 2*P::degree/k + 1);
			tmp = alphas[3*k+i]*tmp;
			vprime = vprime + tmp;
		}

		P const_d, ls;
		fill(const_d.begin(), const_d.end(), (V) P::degree);
		fill(ls.begin(), ls.end(), (V)LS);

		p.h = g;
		for(int _mu = 0; _mu<k;_mu++){
			P Xmu;
			Xmu.data()[_mu] = kinv;
			Xmu.ntt_pow_phi(); // find faster way of getting X^mu

			pis[0] = pka * gammas[2*_mu] + pkb * gammas[2*_mu+1];
			pis[1] = t * gammas[2*_mu];
			pis[2] = t * gammas[2*_mu+1];
			pis[3] = gammas[2*_mu+1];

			for(auto &pi : pis){
				pi.invntt_pow_invphi();
			}

			mpz_t ug;
			mpz_init(ug);

			gammas[2*_mu].invntt_pow_invphi();
			gammas[2*_mu+1].invntt_pow_invphi();

			poly_inner_product_Zq(ug, p.ct.c0, gammas[2*_mu]);
			poly_inner_product_Zq(ug, p.ct.c1, gammas[2*_mu+1]);
			mpz_mod_ui(ug, ug, P::get_modulus(0));


			tmp = inner_product(pis.begin(), pis.end(), shat.begin(), P{0});
			tmp = const_d * tmp;
			tmp2 = gammas_prime[_mu]*(const_d * message - ls);
			tmp = tmp2 + tmp;
			tmp.invntt_pow_invphi();

			tmp = tmp - P{(V)mpz_get_ui(ug)};
			tmp2 = tmp;
			tmp3 = {0};
			for(int i=0;i<k;i++){
				tmp = tmp2;
				for(int j=0;j<i;j++) galois_transform(tmp, tmp, 2*P::degree/k + 1);
				tmp3 = tmp3 + tmp;
			}

			tmp3.ntt_pow_phi();
			tmp3 = tmp3 * Xmu;
			p.h = p.h + tmp3;


			for(auto i=0;i<k;i++){
				tmp3 = {0};
				for(auto _nu=0;_nu<k;_nu++){
					tmp2 = {0};
					for(auto j=0; j<4;j++){
						tmp = inner_product( commkey.b[j].begin(), commkey.b[j].end(), y[nmod(i-_nu, k)].rem.begin(), P{0});
						tmp = tmp + y[nmod(i-_nu, k)].re[j];
						tmp = tmp * pis[j];
						tmp = tmp * const_d;
						tmp2 = tmp2 + tmp;
					}
					tmp = inner_product( commkey.b[3].begin(), commkey.b[3].end(), y[nmod(i-_nu, k)].rem.begin(), P{0});
					tmp = tmp + y[nmod(i-_nu, k)].re[3];
					tmp = const_d*tmp*gammas_prime[_mu];
					tmp2 = tmp2 + tmp;

					tmp2.invntt_pow_invphi();
					tmp = tmp2;
					for(int j=0;j<_nu;j++) galois_transform(tmp, tmp, 2*P::degree/k + 1);
					tmp3 = tmp3 + tmp;
				}
				tmp3.ntt_pow_phi();
				tmp3 = tmp3 * Xmu;
				vs[i] = vs[i] + tmp3;
			}
		}
		p.h.invntt_pow_invphi();

		shake128_inc_init(&state);
		shake128_inc_absorb(&state, whash,  2*SYMBYTES);
		shake128_inc_absorb(&state, (uint8_t *) t5.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) t6.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) t7.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) p.h.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) v.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) vprime.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) v2.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) vs.data(), k*sizeof(P));
		shake128_inc_finalize(&state);
		shake128_inc_squeeze(chash,SYMBYTES, &state);

		P challenge = sample_ternary(0x80, chash);

		for(int i=0;i<k; i++){
			tmp = challenge;
			for(int j=0;j<i;j++) galois_transform(tmp, tmp, 2*P::degree/k + 1);
			tmp.ntt_pow_phi();
			for(size_t j=0;j<params::LAMBDA;j++) p.z[i].rem[j] = y[i].rem[j] + tmp * r_comm.rem[j];
			for(size_t j=0;j<params::MU;j++) p.z[i].rs[j] = y[i].rs[j] + tmp * r_comm.rs[j];
			for(size_t j=0;j<kappa;j++) p.z[i].re[j] = y[i].re[j] + tmp * r_comm.re[j];
		}

		memcpy(p.chash, chash, SYMBYTES);

		for(int i=0; i<k; i++){
			for(size_t j=0;j<params::LAMBDA; j++) p.z[i].rem[j].invntt_pow_invphi();
			for(size_t j=0;j<params::MU; j++) p.z[i].rs[j].invntt_pow_invphi();
			for(size_t j=0;j<kappa; j++) p.z[i].re[j].invntt_pow_invphi();
		}

		rej = true;
		for(int i=0; i<k; i++){
			auto nrm = inf_norm(p.z[i].rem);
			if( nrm >=  delta_1 - beta_1 ){
				rej &= false;
			}
			nrm = inf_norm(p.z[i].rs);
			if( nrm >= delta_1 - beta_1){
				rej &= false;
			}
			nrm = inf_norm(p.z[i].re);
			if( nrm >= delta_1 - beta_1){
				rej &= false;
			}
		}
	}

	p.t.t[LC+5] = t5;
	p.t.t[LC+6] = t6;
	p.t.t[LC+7] = t7;

	p.ct.c0.ntt_pow_phi();
	p.ct.c1.ntt_pow_phi();

//	for(uint8_t *c = chash; c < chash+SYMBYTES; c++){
//		printf("%02x", *c);
//	}printf("\n");

	free(alphas);
	free(alphas_prime);
	free(gammas);
	free(gammas_prime);

}

bool verify(proof_t &p, const comm_key_t &commkey, const pk_t &pubkey){


	bool  flag = true;
	for(int i=0; i<k; i++){
		auto nrm = inf_norm(p.z[i].rem);

		if( nrm >=  delta_1 - beta_1 ){
			flag &= false;
		}
		nrm = inf_norm(p.z[i].rs);

		if( nrm >= delta_1 - beta_1){
			flag &= false;
		}
		nrm = inf_norm(p.z[i].re);

		if( nrm >= delta_1 - beta_1){
			flag &= false;
		}
	}

	for(int i=0; i<k; i++){
		for(size_t j=0;j<params::LAMBDA; j++) p.z[i].rem[j].ntt_pow_phi();
		for(size_t j=0;j<params::MU; j++) p.z[i].rs[j].ntt_pow_phi();
		for(size_t j=0;j<kappa; j++) p.z[i].re[j].ntt_pow_phi();
	}

	if(!flag) return false;
//	cout << "check  1 ✅\n";

	uint8_t symbuf[4*SYMBYTES];
	memset(symbuf, 0, 4*SYMBYTES);
	uint8_t *thash = symbuf; // to store hash of public variables before while loop
	uint8_t *whash = symbuf+2*SYMBYTES; // to store hash of variables withing while loop i.e depends on w
	uint8_t *chash = symbuf+3*SYMBYTES;

	shake128incctx state;

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, (uint8_t *)p.t.t0.data(), params::MU*sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)p.t.t.data(),  (4+LC+1)*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(thash,SYMBYTES, &state);

	array<P, k> sigma_c;
	P tmp;

	P c = sample_ternary(0x80, p.chash);

//	p.c.invntt_pow_invphi();
	for(int i=0;i<k;i++){
		tmp = c;
		for(int j=0;j<i;j++) galois_transform(tmp, tmp, 2*P::degree/k + 1);
		tmp.ntt_pow_phi();
		sigma_c[i] = tmp;
	}
	c.ntt_pow_phi();

	array<array<P, params::MU>, k> w;
	for(int i=0; i<k;i++){
		for(size_t j=0; j<params::MU;j++){
			w[i][j] = inner_product(commkey.B0l[j].begin(), commkey.B0l[j].end(), p.z[i].rem.begin(), P{0});
			tmp = inner_product(commkey.B0k[j].begin(), commkey.B0k[j].end(), p.z[i].re.begin(), P{0});
			w[i][j] = w[i][j] + tmp + p.z[i].rs[j];
			w[i][j] = w[i][j] - sigma_c[i]* p.t.t0[j];
		}
	}

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, thash,  SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *)w.begin(), k*params::MU*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(whash,2*SYMBYTES, &state);

	P* alphas = alloc_aligned<P, 32>(4*k, P{0});
	P* alphas_prime = alloc_aligned<P, 32>(LC*k, P{0});
	P* gammas = alloc_aligned<P, 32>(2*k, P{0});
	P* gammas_prime = alloc_aligned<P, 32>(k, P{0});

	sample_uniform(alphas, whash, 4*k);
	sample_uniform(alphas_prime, chash, LC*k);

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, thash,  SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *)alphas, 4*k*sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)alphas_prime, LC*k*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(whash,2*SYMBYTES, &state);

	sample_uniform(gammas, whash, 2*k);
	sample_uniform_scalar_ntt(gammas_prime, chash, k);



	P f2, f3, f4, fji;
	f2 = inner_product(commkey.b[LC+5].begin(), commkey.b[LC+5].end(), p.z[0].rem.begin(), P{0});
	f2 = f2 + p.z[0].re[LC+5];
	f2 = f2 - c * p.t.t[LC+5];

	f3 = inner_product(commkey.b[LC+6].begin(), commkey.b[LC+6].end(), p.z[0].rem.begin(), P{0});
	f3 = f3 + p.z[0].re[LC+6];
	f3 = f3 - c * p.t.t[LC+6];

	P tmp2, tmp3;
	tmp = {0};
	for(int i=0; i<k;i++){
		for(size_t j=0;j<3;j++){
			fji = inner_product(commkey.b[j].begin(), commkey.b[j].end(), p.z[i].rem.begin(), P{0});
			fji = fji + p.z[i].re[j];
			fji = fji - sigma_c[i] * p.t.t[j];
			tmp2 = fji * (fji - sigma_c[i]) * (fji + sigma_c[i]);
			for(int ii=0; ii<i; ii++) inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
			tmp2 = tmp2 * alphas[i*3+j];
			tmp  = tmp + tmp2;
		}
	}
	P v = tmp + f2 + c*f3;

	f4 = inner_product(commkey.b[LC+7].begin(), commkey.b[LC+7].end(), p.z[0].rem.begin(), P{0});
	f4 = f4 + p.z[0].re[LC+7];
	f4 = f4 - c * p.t.t[LC+7];

	tmp = {0};
	P f3i;
	for(int i=0; i<k;i++){
		f3i = inner_product(commkey.b[3].begin(), commkey.b[3].end(), p.z[i].rem.begin(), P{0});
		f3i = f3i + p.z[i].re[3];
		f3i = f3i - sigma_c[i] * p.t.t[3];
		for(size_t j=0;j<LC;j++){
			fji = inner_product(commkey.b[3+j].begin(), commkey.b[3+j].end(), p.z[i].rem.begin(), P{0});
			fji = fji + p.z[i].re[3+j];
			fji = fji - sigma_c[i] * p.t.t[3+j];
			fill(tmp2.begin(), tmp2.end(), (P::value_type) (j+1));
			tmp2 = sigma_c[i] * tmp2;
			tmp2 = f3i + tmp2;
			tmp2 = fji * tmp2;
			fji = inner_product(commkey.b[3+j+1].begin(), commkey.b[3+j+1].end(), p.z[i].rem.begin(), P{0});
			fji = fji + p.z[i].re[3+j+1];
			fji = fji - sigma_c[i] * p.t.t[3+j+1];
			tmp2 = tmp2 + sigma_c[i] * fji;
			for(int ii=0; ii<i; ii++) inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
			tmp2 = tmp2 * alphas_prime[i*LC + j];
			tmp  = tmp + tmp2;
		}
	}
	P v2 = tmp + f4;


	for(int i=0;i<k;i++){
		if(p.h.data()[i] != 0) return false;
	}
//	cout << "check  2 ✅\n";


	array<P, k> vs;
	P vprime;
	for(size_t i=0; i<k;i++){
		vs[i] = inner_product(commkey.b[LC+4].begin(), commkey.b[LC+4].end(), p.z[i].rem.begin(), P{0});
		vs[i] = vs[i] + p.z[i].re[LC+4];

		tmp = inner_product(commkey.b[LC+3].begin(), commkey.b[LC+3].end(), p.z[i].rem.begin(), P{0});
		tmp = tmp + p.z[i].re[LC+3];
		tmp = tmp - sigma_c[i] * p.t.t[LC+3];
		tmp2 = tmp;
		for(size_t j=0;j<i;j++) inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
		tmp2 = alphas[3*k+i]*tmp2;
		vprime = vprime + tmp2;
	}

	array<P, 4> pis;
	P pka = pubkey.a;
	pka.invntt_pow_invphi();
	poly_transpose(pka);
	pka.ntt_pow_phi();

	P pkb = pubkey.b;
	pkb.invntt_pow_invphi();
	poly_transpose(pkb);
	pkb.ntt_pow_phi();

	//P Xmu, X{0,1};
	//X.ntt_pow_phi();
	//Xmu = one;
	V kinv = (V) modInverse((P::signed_value_type) k, (P::signed_value_type) P::get_modulus(0)); // k * kinv = 1 mod q
	p.ct.c0.invntt_pow_invphi();
	p.ct.c1.invntt_pow_invphi();

	P const_d, ls;
	fill(const_d.begin(), const_d.end(), (V) P::degree);
	fill(ls.begin(), ls.end(), (V) LS);

	P t = {params::plaintextModulus};
	t.ntt_pow_phi();
	P tau;
	for(int _mu = 0; _mu<k;_mu++){
		P Xmu;
		Xmu.data()[_mu] = kinv;
		Xmu.ntt_pow_phi(); // find faster way of getting X^mu

		pis[0] = pka * gammas[2*_mu] + pkb * gammas[2*_mu+1];
		pis[1] = t * gammas[2*_mu];
		pis[2] = t * gammas[2*_mu+1];
		pis[3] = gammas[2*_mu+1];

		for(auto &pi : pis){
			pi.invntt_pow_invphi();
		}

		mpz_t ug;
		mpz_init(ug);

		gammas[2*_mu].invntt_pow_invphi();
		gammas[2*_mu+1].invntt_pow_invphi();

		poly_inner_product_Zq(ug, p.ct.c0, gammas[2*_mu]);
		poly_inner_product_Zq(ug, p.ct.c1, gammas[2*_mu+1]);
		mpz_mod_ui(ug, ug, (unsigned long int) P::get_modulus(0));

		tmp = inner_product(pis.begin(), pis.end(), p.t.t.begin(), P{0});
		tmp = const_d * tmp;
		tmp2 = gammas_prime[_mu]*(const_d*p.t.t[3] - ls);
		tmp = tmp + tmp2;
		tmp.invntt_pow_invphi();

		tmp = tmp - P{(V)mpz_get_ui(ug)};
		tmp2 = tmp;
		tmp3 = {0};
		for(int i=0;i<k;i++){
			tmp = tmp2;
			for(int j=0;j<i;j++) galois_transform(tmp, tmp, 2*P::degree/k + 1);
			tmp3 = tmp3 + tmp;
		}
		tmp3.ntt_pow_phi();
		tmp3 = tmp3 * Xmu;
		tau = tau + tmp3;

		for(auto i=0;i<k;i++){
			tmp3 = {0};
			for(auto _nu=0;_nu<k;_nu++){
				tmp2 = {0};
				for(auto j=0; j<4;j++){
					tmp = inner_product( commkey.b[j].begin(), commkey.b[j].end(), p.z[nmod(i-_nu, k)].rem.begin(), P{0});
					tmp = tmp + p.z[nmod(i-_nu, k)].re[j];
					tmp = tmp * pis[j];
					tmp = tmp * const_d;
					tmp2 = tmp2 + tmp;
				}
				tmp = inner_product( commkey.b[3].begin(), commkey.b[3].end(), p.z[nmod(i-_nu, k)].rem.begin(), P{0});
				tmp = tmp + p.z[nmod(i-_nu, k)].re[3];
				tmp = const_d*tmp*gammas_prime[_mu];
				tmp2 = tmp2 + tmp;

				tmp2.invntt_pow_invphi();
				tmp = tmp2;
				for(int j=0;j<_nu;j++) galois_transform(tmp, tmp, 2*P::degree/k + 1);
				tmp3 = tmp3 + tmp;
			}
			tmp3.ntt_pow_phi();
			tmp3 = tmp3 * Xmu;
			vs[i] = vs[i] + tmp3;
		}
	}

	p.h.ntt_pow_phi();
	tmp2 = tau + p.t.t[LC+4] - p.h;
	for(auto i=0; i<k;i++){
		vs[i] = vs[i] - sigma_c[i] * tmp2;
	}
	p.h.invntt_pow_invphi();

	p.ct.c0.ntt_pow_phi();
	p.ct.c1.ntt_pow_phi();




	shake128_inc_init(&state);
	shake128_inc_absorb(&state, whash,  2*SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *) p.t.t[LC+5].data(), 3*sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *) p.h.data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *) v.data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *) vprime.data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *) v2.data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *) vs.data(), k*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(chash,SYMBYTES, &state);

//	for(uint8_t *c = chash; c < chash+SYMBYTES; c++){
//		printf("%02x", *c);
//	}printf("\n");


	P challenge = sample_ternary(0x80, chash);
	challenge.ntt_pow_phi();
	if(challenge != c) return false;
//	cout << "check  3 ✅\n";

	free(alphas);
	free(alphas_prime);
	free(gammas);
	free(gammas_prime);
	return true;
}

int main(){

	cout << "sizeof proof:" << sizeof(proof_t) << "B\n";

	srand(time(NULL));

	sk_t secret_key;
	pk_t public_key(secret_key);

	comm_key_t commkey = comm_key_t();

	auto total_duration_voter=0, total_duration_verifier=0, total_duration_tally=0;
	int rejected_ballots = 1;

	array<V, NC> valid_ballot = ballot_generator();

	BGV::ciphertext_t tally;
	P real_tally = {0};
	int i = 0;
	for(i=0;i<NV;i++){
		P m = {0};
		random_ballot(m, valid_ballot);

		proof_t proof;
		auto start = chrono::steady_clock::now();
		vericrypt(proof, public_key, commkey, m);
		auto end = chrono::steady_clock::now();
		auto duration_voter = chrono::duration_cast<chrono::milliseconds>(end-start).count();
		//cout << "verifiable encryption time: " << duration << "(ms)" << endl;

		start = chrono::steady_clock::now();
		bool res = verify(proof, commkey, public_key);
		end = chrono::steady_clock::now();
		auto duration_verifier = chrono::duration_cast<chrono::milliseconds>(end-start).count();

		if(!res) {
			cout << "rejected ballot #" << (rejected_ballots++) << endl;
			continue;
		}

		total_duration_voter += duration_voter;
		total_duration_verifier += duration_verifier;
		start = chrono::steady_clock::now();
		tally += proof.ct;
		end = chrono::steady_clock::now();
		auto duration_tally = chrono::duration_cast<chrono::microseconds>(end-start).count();
		total_duration_tally += duration_tally;
		real_tally = real_tally + m;
		//cout << "single encryption noise: " << BGV::noise(m, secret_key, public_key, proof.ct) << endl;
		//cout << "proof size:" << sizeof(proof) << "B\n";
	}

	cout << "total encryption time: " << total_duration_voter << "(ms)" << endl;
	cout << "total verifier time: " << total_duration_verifier << "(ms)" << endl;
	cout << "total tallying time: " << total_duration_tally << "(μs)" << endl;
	cout << "rejected ballots: " << rejected_ballots-1 << endl;

	P result;
	auto start = chrono::steady_clock::now();
	BGV::decrypt_poly(result, secret_key, public_key, tally);
	auto end = chrono::steady_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(end-start).count();
	cout << "Decryption time:" << duration << "(μs)" << endl;
	cout << "Decrypted tally is " << (result == real_tally ? "" : "not") << "equal to actual tally" << endl;
	//cout << result << endl;

	auto noise = BGV::noise(result, secret_key, public_key, tally);
	cout << "log2 noise:" << noise << endl;
	cout << "noise max: " << public_key.noise_max << endl;

	return 0;
}



/*
 * PARI code
 * \p 200
default(parisize, 12000000000);
pola = ((2/9)*X^-1 + 5/9 + (2/9)*X^1)^4096;
polb = (1/3)*X^-1 + 1/3 + (1/3)*X^1;
polc = pola*pola*polb;

s=0; for(n=-768,768, s=s+polcoeff(polc, n));
ans = log(1-s)/log(2) // -121.140540387017417895473837357...

pold = polc^11000;
 *
 *
 * */

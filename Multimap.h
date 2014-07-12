#ifndef __MULTIMAP_H
#define __MULTIMAP_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <gmpxx.h>
#include <assert.h>
#include <string>
#include <sys/time.h>

#define INSTANTIATION 1 		// 1- Small, 2- Medium, 3- Large, 4- Extra

#define kappa 6					// Maximal level
#define hBits 80 				// size of h_i's in v = BETA
#define theta  15 				// number of non-zero elements in subset sum during rerand 

#define sessionKeyBits 160 		// Bitsize of session key to derive
#define bound sessionKeyBits 	// bound to decide if it is zero or not
								// bound must be >= sessionKeyBits
#define alpha 80				// size of g_i's and elements of A

#if INSTANTIATION == 1 		// Small

	#define N  540  			// number of p_i's
	#define delta 23 			// sqrt(N)
	#define eta 1838			// size of p_i's
	#define etp 460				// Size of primes in p_i's; it is better to have eta%etp=0
	#define rho 41				// size of r_i's in xp_i's and y

#elif INSTANTIATION == 2 		// Medium

	#define N  2085  			// number of p_i's
	#define delta 45 			// sqrt(N)
	#define eta 2043			// size of p_i's
	#define etp 409				// Size of primes in p_i's; it is better to have eta%etp=0
	#define rho 56				// size of r_i's in xp_i's and y

#elif INSTANTIATION == 3 		// Large

	#define N  8250  			// number of p_i's
	#define delta 90 			// sqrt(N)
	#define eta 2261			// size of p_i's
	#define etp 453				// Size of primes in p_i's; it is better to have eta%etp=0
	#define rho 72				// size of r_i's in xp_i's and y

#elif INSTANTIATION == 4 		// Extra

	#define N  26115  			// number of p_i's
	#define delta 161 			// sqrt(N)
	#define eta 2438			// size of p_i's
	#define etp 407				// Size of primes in p_i's; it is better to have eta%etp=0
	#define rho 85				// size of r_i's in xp_i's and y

#endif



double currentTime();

class MMKey;

/* 
Class Ciphertext

Contain: 
- ciphertext value (large integer) `cval'
- ciphertext degree `degree'
- pointer to the key associated with the ciphertext `key'
*/
class Ciphertext {
private:
	mpz_class cval;
	long degree;
	MMKey* key;

public:
	Ciphertext();
	Ciphertext(const Ciphertext& c);
	Ciphertext(MMKey* mmkey, mpz_class c, long deg);

	void Decrypt_with_sk(mpz_class* m);
	long get_noise();

	long get_degree() const {return degree;};
	mpz_class get_cval() const {return cval;};

	mpz_class deriveSessionKey();

	Ciphertext& operator=(const Ciphertext&);
	Ciphertext& operator+=(const Ciphertext&);
	Ciphertext& operator+=(const mpz_class&);
	Ciphertext& operator*=(const Ciphertext&);
	Ciphertext& operator-=(const Ciphertext&);
	Ciphertext& operator-=(const mpz_class&);
	
	Ciphertext operator+(const Ciphertext& c) const {
	    Ciphertext c2(*this);
	    return (c2 += c);
	}
	Ciphertext operator+(const mpz_class a) const {
	    Ciphertext c2(*this);
	    return (c2 += a);
	}
	friend Ciphertext operator+(const mpz_class, const Ciphertext&);

	Ciphertext operator-(const Ciphertext& c) const {
	    Ciphertext c2(*this);
	    return (c2 -= c);
	}
	Ciphertext operator-(const mpz_class a) const {
	    Ciphertext c2(*this);
	    return (c2 -= a);
	}
	friend Ciphertext operator-(const mpz_class, const Ciphertext&);

	Ciphertext operator*(Ciphertext& c) const {
	    Ciphertext c2(*this);
	    return (c2 *= c);
	}
};

/* 
Class MMKey (Multilinear-Map Key)

Contain: 
- pointer to the gmp pseudorandom generator `rng'
- pointer to secret primes `p'
- public key value `x0' (=prod(p))
- private value `z' and `zkappa'=z^kappa
- private value `zinv' = z^(-1) mod x0
- pointer to private elements `g'
- public value `y'
- public zero-tester `v'
*/
class MMKey {
private:
	gmp_randclass* 	rng;
	mpz_class* 		p; 				//	[N];
	mpz_class 		x0, z, zkappa;
	mpz_class 		zinv; 			//	[N];
	mpz_class* 		crtCoeff; 		//	[N];
	mpz_class* 		g; 				//	[N];
	mpz_class 		y;
	mpz_class 		v;

public:
	MMKey(gmp_randclass* random);
	~MMKey();
	Ciphertext Sample(unsigned long k);
	mpz_class Encrypt_with_sk(mpz_class* m, long nbBits, long degree);
	mpz_class Encrypt_with_sk(unsigned long m, long nbBits, long degree);
	void Decrypt_with_sk(mpz_class* m, const Ciphertext& c);
	mpz_class reduce(const mpz_class &c);
	Ciphertext get_y();
	long get_noise(const mpz_class& c, long degree);
	mpz_class zero_test(const mpz_class &c, long degree);
	long nbBits(const mpz_class &v);
	bool is_zero(const Ciphertext &c);
	Ciphertext Rerand(const Ciphertext & c);
	mpz_class& get_x0() { return x0; };
};

#endif 
// #ifndef __MULTIMAP_H

#include <iostream>
#include <cstdlib>
#include <gmpxx.h>
#include <assert.h>
#include <sys/time.h>
#include "Multimap.h"

#define LINUX false			// Timings using CLOCK_MONOTONIC
#define VERBOSE false			// Display additional information
#define DISPLAY_MESSAGES false		// Decrypt & display the messages (only if VERBOSE==true)

#define USERS kappa+1			// Number of Users: make sure (USERS-1) <= kappa

std::ostream& operator<<(std::ostream& os, Ciphertext& c) {
    os << "<Ciphertext of degree=" << c.get_degree();

    #if VERBOSE
    os << " and with bit-noise=" << c.get_noise() << "/" << (eta-alpha);
	
	#if DISPLAY_MESSAGES
	unsigned i;
	mpz_class* m;
	m = new mpz_class[N];
  	c.Decrypt_with_sk(m);
	os << std::endl;
	os << "          m=(";
    for (i=0; i<N; i++)
    	os << m[i] << ((i==(N-1))?") ":" ");
    #endif

    #endif 
    // #if VERBOSE
    os << ">";
    return os;
}



int main()
{
	long i;

#if LINUX
	timespec start, finish;
#else
	double startTime;
#endif

	// PRNG
	gmp_randclass* random = new gmp_randclass(gmp_randinit_default);

	// KeyGen
#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &start);
#else
	startTime=currentTime();
#endif

	MMKey key(random);
	Ciphertext y = key.get_y();

#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &finish);
	std::cout << "Keygen: " << (double)((finish.tv_nsec-start.tv_nsec)/1000000000.0+(finish.tv_sec-start.tv_sec)) << "s" << std::endl;
#else
	std::cout << "Keygen: " << (double)(currentTime()-startTime) << "s" << std::endl;
#endif

	std::cout << "y=" << y << std::endl;
	std::cout << std::flush;
	
	// USERS-multipartite Diffie Hellman
	Ciphertext* l0[USERS];
	Ciphertext* l1[USERS];
	Ciphertext* keys[USERS];
	
	long j;
	bool b[USERS][ell];
	
#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &start);
#else
	startTime=currentTime();
#endif

	for (j=0; j<USERS; j++)
	{
		for (i=0; i<ell; i++)
			b[j][i] = rand()%2;
		l0[j] = new Ciphertext(key.Encrypt(b[j])); // Secret value of each user

		std::cout << "User #" << j << " (level 0): " << *l0[j] << std::endl; 
	}
	std::cout << std::endl;

	for (j=0; j<USERS; j++)
	{
		Ciphertext c = (*l0[j])*y;
		l1[j] = new Ciphertext(key.Rerand(c)); // Public value of each user

		std::cout << "User #" << j << " (level 1): " << *l1[j] << std::endl; 
	}
	std::cout << std::endl;
#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &finish);
	std::cout << "(Samp+enc+rerand): " << (double)((finish.tv_nsec-start.tv_nsec)/1000000000.0+(finish.tv_sec-start.tv_sec)) << "s" << std::endl;
#else
	std::cout << "(Samp+enc+rerand): " << (double)(currentTime()-startTime) << "s" << std::endl;
#endif
	std::cout << std::flush;


#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &start);
#else
	startTime=currentTime();
#endif

	for (j=0; j<USERS; j++)
	{
		keys[j] = new Ciphertext((*l0[j]));
		for (i=0; i<USERS; i++)
			if (i != j) *keys[j] *= (*l1[i]); 
			// Each user compute the (USERS-1) products using public keys of the other users

		std::cout << "User #" << j << " (level " << (USERS-1) << "): " << *keys[j] << std::endl; 
	}
	std::cout << std::endl;

	std::cout << "Session keys of " << bound << " bits" << std::endl;
	for (j=0; j<USERS; j++)
		std::cout << "User #" << j << ": " << (*keys[j]).deriveSessionKey() << std::endl;

#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &finish);
	std::cout << "(Product+Extract): " << (double)((finish.tv_nsec-start.tv_nsec)/1000000000.0+(finish.tv_sec-start.tv_sec)) << "s" << std::endl;
#else
	std::cout << "(Product+Extract): " << (double)(currentTime()-startTime) << "s" << std::endl;
#endif

std::cout << std::flush;


// How long could have been the sessionKey?
#if VERBOSE
	std::cout << "Tests if keys[i]-keys[j] \"=\" 0: " << std::endl;
	long maxBits = 0, keyBits;
	#pragma omp parallel for private(j, keyBits)
	for (i=0; i<USERS; i++)
		for (j=0; j<USERS; j++)
		{
			if (i != j) 
			{
				Ciphertext diff(*keys[i]);
				diff-=(*keys[j]);
				keyBits = key.nbBits(key.zero_test(diff.get_cval(), diff.get_degree()));
				#pragma omp critical
				{
					if (keyBits > maxBits) maxBits = keyBits;
				}
				std::cout << ((keyBits < key.nbBits(key.get_x0())-bound) ? "true" : "false") << " ";
				std::cout << std::flush;
			}
		}
	std::cout << std::endl;
	std::cout << "Maximum size of session key " << (key.nbBits(key.get_x0()) - maxBits) << " bits" << std::endl;
#endif


	return 0;
}

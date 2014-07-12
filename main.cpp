#include <iostream>
#include <cstdlib>
#include <gmpxx.h>
#include <assert.h>
#include <sys/time.h>
#include "Multimap.h"

#define LINUX false			// Timings using CLOCK_MONOTONIC
#define VERBOSE false			// Display additional information
#define DISPLAY_MESSAGES false		// Decrypt & display the messages (only if VERBOSE==true)


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
	
	// kappa level-1 encodings to generate and multiply
	Ciphertext* encodings[kappa];
	Ciphertext* result;
	
	long j;
	
#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &start);
#else
	startTime=currentTime();
#endif

	// generate kappa random level-1 encodings
	//TODO
	for (j=0; j<kappa; j++)
	{
		// generate random level-1 encoding
		encodings[j] = new Ciphertext(key.Sample(1)); // Secret value of each user

	}

#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &finish);
	std::cout << "(Encoding generation): " << (double)((finish.tv_nsec-start.tv_nsec)/1000000000.0+(finish.tv_sec-start.tv_sec)) << "s" << std::endl;
#else
	std::cout << "(Encoding generation): " << (double)(currentTime()-startTime) << "s" << std::endl;
#endif
	std::cout << std::flush;


#if LINUX
	clock_gettime(CLOCK_MONOTONIC, &start);
#else
	startTime=currentTime();
#endif

	// multiply all encodings

	result = new Ciphertext((*encodings[0]));
	for (j=1; j<kappa; j++)
	{
		*result *= (*encodings[j]);

	}

	// perform extraction

	std::cout << "Extracted value:" << (*result).deriveSessionKey() << std::endl;

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

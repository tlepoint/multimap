#include "Multimap.h"

/*
	Return current time in sec.
*/
double currentTime()
{
	timeval t;
	gettimeofday(&t,NULL);
	return (double) (t.tv_sec+(double)(t.tv_usec/1000000.0));
}

/*
	Generate a (centered) random of nbBits using rng
*/
inline mpz_class generateRandom(long nbBits, gmp_randclass* rng)
{
	if (nbBits <= 1) return rand()%2;
	else return rng->get_z_bits(nbBits)-(mpz_class(1)<<(nbBits-1));
}

/*
	Compute a mod b, a>=0
*/
inline mpz_class mod(const mpz_class &a, const mpz_class &b)
{
	mpz_class res;
	mpz_mod (res.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
	return res;
}

/*
	Compute a mod b, -b/2 < a <= b/2
*/
inline mpz_class modNear(const mpz_class &a, const mpz_class &b)
{
	mpz_class res;
	mpz_mod (res.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
	if (res > (b>>1)) res -= b;
	return res;
}

/*
	Compute nearest value of a/b
*/
inline mpz_class quotNear(const mpz_class &a, const mpz_class &b)
{
	return (a-modNear(a,b))/b;
}



/*
	Initialization of a ciphertext to 0
*/
Ciphertext::Ciphertext()
{
	cval = 0;
	degree = 0;
}

/*
	Initialization of a ciphertext from another ciphertext
*/
Ciphertext::Ciphertext(const Ciphertext& c)
{
	key = c.key;
	cval = c.cval;
	degree = c.degree;
}

/*
	Initialization of a ciphertext given a pointer to a MMKey, a value c and a degree deg
*/
Ciphertext::Ciphertext(MMKey* mmkey, mpz_class c, long deg)
{
	key = mmkey;
	cval = c;
	degree = deg;
}

/*
	Decrypt ciphertext c, called from c itself
*/
void Ciphertext::Decrypt_with_sk(mpz_class* m)
{
	key->Decrypt_with_sk(m, *this);
}

/*
	Get c noise from itself, calling the secret key
*/
long Ciphertext::get_noise()
{
	return key->get_noise(cval, degree);
}

/*
	Derive the session Key from c
*/
mpz_class Ciphertext::deriveSessionKey()
{
	return ((key->zero_test(cval, degree))>>(key->nbBits(key->get_x0())-sessionKeyBits));
}

/*
	Operations on the ciphertexts: +, -, *
*/
Ciphertext& Ciphertext::operator+=(const Ciphertext& c) {
	assert(degree == c.degree);
	cval += c.cval;
	cval = key->reduce(cval);
	return *this;
}

Ciphertext& Ciphertext::operator+=(const mpz_class& c) {
	cval += c;
	cval = key->reduce(cval);
	return *this;
}

Ciphertext& Ciphertext::operator-=(const Ciphertext& c) {
	assert(degree == c.degree);
	cval -= c.cval;
	cval = key->reduce(cval);
	return *this;
}

Ciphertext& Ciphertext::operator-=(const mpz_class& c) {
	cval -= c;
	cval = key->reduce(cval);
	return *this;
}

Ciphertext& Ciphertext::operator*=(const Ciphertext& c) {
	if (degree+c.degree>0)
		degree += c.degree;
	cval *= c.cval;
	cval = key->reduce(cval);
	return *this;
}

/*
	Generation of a MMKey given a pointer to the random generator
*/
MMKey::MMKey(gmp_randclass* random) {

	long i, j;
	double startTime;

	// Define PRNG
	rng = random;
	
	// Generate the p_i's
	// /!\ Generating primes of eta bits can be very long, 
	// so eventually we generate p as product of primes of etp bits
	x0 = 1;
	long niter=(unsigned)eta/etp;
	long psize;

	startTime=currentTime();
	std::cout << "Generate the p_i's and x0: " << std::flush;
	mpz_class p_tmp, p_unif;
	p = new mpz_class[N];
	#pragma omp parallel for private(p_tmp, p_unif, j)
	for (i=0; i<N; i++)
	{
		p[i] = 1;
		for (j=0; j<niter; j++)
		{
		  if (j<(niter-1)) psize=etp;
		  else psize=eta-etp*(niter-1);
		  
		  p_unif = rng->get_z_bits(psize);
		  
		  mpz_nextprime(p_tmp.get_mpz_t(), p_unif.get_mpz_t());
		  p[i] *= p_tmp;
		}
		#pragma omp critical
		{
			x0 *= p[i];
		}
	}
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;

	// Generate the CRT Coefficients
	// /!\ This requires a lot of RAM (>100 GB for extra instantiation)
	startTime=currentTime();
	std::cout << "Generate the crtCoeff_i's: " << std::flush;
	mpz_class Q;
	crtCoeff = new mpz_class[N];
	#pragma omp parallel for private(Q)
	for (i=0; i<N; i++)
	{
		Q = x0/p[i];
		mpz_invert(crtCoeff[i].get_mpz_t(), Q.get_mpz_t(), p[i].get_mpz_t());
		crtCoeff[i] *= Q;
	}
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;

	// Generate the g_i's
	startTime=currentTime();
	std::cout << "Generate the g_i's: " << std::flush;
	mpz_class g_tmp;
	g = new mpz_class[N];
	#pragma omp parallel for private(g_tmp)
	for (i=0; i<N; i++)
	{
		g_tmp = rng->get_z_bits(alpha);
		mpz_nextprime(g[i].get_mpz_t(), g_tmp.get_mpz_t());
	}
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;

	// Generate z
	startTime=currentTime();
	std::cout << "Generate z and zinv: " << std::flush;
	do {
		z = rng->get_z_range(x0);
		mpz_invert(zinv.get_mpz_t(), z.get_mpz_t(), x0.get_mpz_t());
	}
	while (zinv == 0);
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;
	
	// Generate A
	std::cout << "Generate A: " << std::flush;
	startTime=currentTime();
	A = new mpz_class[ell*N];
	for (i=0; i<ell; i++)
		for (j=0; j<N; j++)
			A[i*N+j] = generateRandom(alpha, rng);
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;

	// Generate the xp_i's
	std::cout << "Generate the xp_i's: " << std::flush;
	startTime=currentTime();
	xp = new mpz_class[ell];
	#pragma omp parallel for
	for (i=0; i<ell; i++)
		xp[i] = Encrypt_with_sk(A+i*N, rho, 0);
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;

	// Generate the varpi_i's
	// Elements used for rerandomization, `x' in the article
	std::cout << "Generate the varpi[j]_i's: " << std::flush;
	startTime=currentTime();
	varpi = new mpz_class[2*delta];

	#pragma omp parallel for private(i)
	for (i=0; i<delta; i++)
	{
		varpi[i] = Encrypt_with_sk((unsigned long) 0, rho, 0);
		varpi[delta+i] = Encrypt_with_sk((unsigned long) alpha, rho, 1);
	}	
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;

	// Generate y
	std::cout << "Generate y: " << std::flush;
	startTime=currentTime();
	y = Encrypt_with_sk((unsigned long) 1, rho, 1);
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;

	// Generate zero-tester v
	std::cout << "Generate the zero-tester v: " << std::flush;
	startTime=currentTime();
	mpz_class input;
	zkappa=1;
	for (i=0; i<kappa; i++)
		zkappa = mod(zkappa*z, x0);
	v=0;
	#pragma omp parallel for private(i, input)
	for (i=0; i<N; i++)
	{
		mpz_invert(input.get_mpz_t(), g[i].get_mpz_t(), p[i].get_mpz_t());
		input = mod(input*zkappa, p[i])*generateRandom(hBits, rng)*(x0/p[i]); 
		#pragma omp critical
		{
			v += input;
		}
	}
	v = mod(v,x0);
	std::cout << (double)(currentTime()-startTime) << "s" << std::endl;
}

/*
	Encrypt bit vector b using public key (subset sum)
*/
Ciphertext MMKey::Encrypt(bool b[ell])
{
	mpz_class c=0;

	for (long i=0; i<ell; i++)
		if (b[i]) c += xp[i];

	return Ciphertext(this, mod(c, x0), 0);
}

/*
	Encrypt input ARRAY `m' with `nbBits' random and initial degree `degree' with secret key
*/
mpz_class MMKey::Encrypt_with_sk(mpz_class* m, long nbBits, long degree)
{
	std::cout << ". " << std::flush;

	long i, j;

	mpz_class res=0;

	for (i=0; i<N; i++)
	{
		res += (m[i] + g[i]*generateRandom(nbBits, rng)) * crtCoeff[i];
	}

	res = mod(res, x0);
	for (j=degree; j>0; j--)
		res = mod(res*zinv, x0);

	std::cout << "* " << std::flush;

	return res;
}

/*
	Destruction of MMKey
*/
MMKey::~MMKey()
{
	//
}

/*
	Encrypt input VALUE `m' with `nbBits' random and initial degree `degree' with secret key
*/
mpz_class MMKey::Encrypt_with_sk(unsigned long m, long nbBits, long degree)
{
	std::cout << ". " << std::flush;

	long i, j;

	mpz_class input, res=0;
	for (i=0; i<N; i++)
	{
		if (m <= 1)
			input = m + g[i]*generateRandom(nbBits, rng);
		else
			input = generateRandom(m, rng) + g[i]*generateRandom(nbBits, rng);

		res += input*crtCoeff[i];
	}

	res = mod(res, x0);
	for (j=degree; j>0; j--)
		res = mod(res*zinv, x0);

	std::cout << "* " << std::flush;

	return res;
}

/*
	Decrypt ciphertext c using the secret key
*/
void MMKey::Decrypt_with_sk(mpz_class* m, const Ciphertext& c)
{
	long i;
	mpz_class value = c.get_cval();
	for (i=c.get_degree(); i>0; i--)
		value = mod(value*z, x0);
	#pragma omp parallel for
	for (i=0; i<N; i++)
	{
		m[i] = modNear(modNear(value, p[i]), g[i]);
	}
}

/*
	Compute c mod x0
*/
mpz_class MMKey::reduce(const mpz_class &c)
{
	return mod(c, x0);
}

/*
	Return public value y converted into a ciphertext
*/
Ciphertext MMKey::get_y()
{
	return Ciphertext(this, y, 1);
}

/*
	Get noise in a ciphertext value c of degree `degree'
*/
long MMKey::get_noise(const mpz_class& c, long degree)
{
	long i;
	mpz_class value = c;
	long max = 0, nbBits;
	mpz_class noise;

	for (i=degree; i>0; i--)
		value *= z;

	#pragma omp parallel for private(noise, nbBits)
	for (i=0; i<N; i++)
	{
		noise = quotNear(modNear(value, p[i]), g[i]);
		nbBits = mpz_sizeinbase(noise.get_mpz_t(),2);
		#pragma omp critical
		{
		if (nbBits>max) max = nbBits;
		}
	}
	return max;
}

/*
	Return w = (c*z^(kappa-degree))*v mod x0
*/
mpz_class MMKey::zero_test(const mpz_class &c, long degree)
{
	assert(degree<=kappa);
	mpz_class value = modNear(c*v,x0);

	for (long i=kappa-degree; i>0; i--)
		value = modNear(value*y, x0);

	return value;
}

/*
	Return number of bits of v
*/
long MMKey::nbBits(const mpz_class &v)
{
	return mpz_sizeinbase(v.get_mpz_t(), 2);
}

/*
	Check whether c is an encoding of 0
*/
bool MMKey::is_zero(const Ciphertext &c)
{
	mpz_class value = zero_test(c.get_cval(), c.get_degree());
	return (nbBits(value)<(nbBits(x0)-bound)) ? 1 : 0;
}

/*
	Rerandomize a ciphertext c
*/
Ciphertext MMKey::Rerand(const Ciphertext & c)
{
	assert(c.get_degree() == 1);
	mpz_class cval = c.get_cval();
	long i, j;

	long indicesvarpi[theta];
	for (i=0; i<theta; i++)
	{
		indicesvarpi[i] = rand()%(delta*delta);
		for (j=0; j<i; j++)
			if (indicesvarpi[j] == indicesvarpi[i]) i--; // restart this indicesvarpi[i].
	}

	for (i=0; i<theta; i++)
		cval += varpi[indicesvarpi[i]%delta] * varpi[delta+indicesvarpi[i]/delta];
	
	return Ciphertext(this, mod(cval, x0), 1);
}

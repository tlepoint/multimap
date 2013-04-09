An implementation of multilinear maps over the integers
=======================================================

This is an implementation of cryptographic multilinear maps, used to perform a
one-round n-party (unauthenticated) Diffie-Hellman key exchange. This
implementation is described in the following article:

[1] J.-S. Coron, T. Lepoint and M. Tibouchi, "Practical Multilinear Maps over
the Integers". Available at http://eprint.iacr.org/2013/183.

This proof-of-concept implementation is done in C++ using the GMP library
http://gmplib.org/


WHAT ARE MULTILINEAR MAPS?
--------------------------

In cryptography, a multilinear map is a mapping between useful cryptographic
groups, which allows the construction of new cryptographic schemes based on
the reduction of one problem in one group to a different problem in the target
group.

A first candidate for multilinear maps has been described by Garg, Gentry and
Halevi in

[2] S. Garg, C. Gentry and S. Halevi, "Candidate Multilinear Maps from Ideal
Lattices and Applications". Proceedings of Eurocrypt 2013, available at
http://eprint.iacr.org/2012/610.

The main difference with bilinear pairings is that the encoding a_i * g of an
element a_i is randomized instead of deterministic; only the final multilinear
map e(a_1*g, ..., a_k*g) is a deterministic function of the  $a_i$'s only.


MULTILINEAR MAPS OVER THE INTEGERS
----------------------------------

This scheme (as in [2]) requires a trusted generator to compute the public
values. This entity knows a master secret that allows her to decode all
encodings (it is similar to a "trapdoor discrete logarithm" system).

Define n secret eta-bit primes p_i of product x0=p_1*...*p_n, n alpha-bit
primes g_i and a secret large integer z invertible mod x_0 (master secret) and
publish the product x0.

A level-k encoding of a vector m=(m_1, ..., m_n) if an integer defined modulo
x_0 such that

c = (r_i * g_i + m_i) / z^k    mod p_i

where the r_i's are small integers.

One can add two encodings at the same level k and obtain and encoding at level
k of the component-wise sum of the m vectors.

On can also multiply two encodings at levels k_1 and k_2 and obtain an
encoding of the component-wise product of the m vectors at level k_1+k_2.

However, the noises r_i grow in the encodings after additions or
multiplications. Therefore the parameters must be chosen so that a target
number of operations can be achieved.

For level kappa encodings, the trusted authority publishes a zero-testing
parameter p_zt which allows to test whether an encoding c at level kappa (only)
is an encoding of 0 or not. The element is given by the formula

p_zt = 	h_1*(z^kappa * g_1^(-1) mod p_1)*p_2*...*p_n + ... 
		+ h_n*(z^kappa * g_n^(-1) mod p_n)*p_1*...*p_(n-1) mod x0

where the h_i's are small integers. Therefore when we multiply a level kappa
encoding c with p_zt we get

p_zt*c = 	h_1*(r_1 + m_1*(g_1^(-1) mod p_1))*p_2*...*p_n + ... 
			+ h_n*(r_n + m_n*(g_n^(-1) mod p_1))*p_1*...*p_(n-1) mod x0

When all the m_i=0 for all i, since the r_i's and h_i's are small, p_zt*c is a
lot smaller than x0. Therefore the MSB of p_zt*c depends only on the m_i's and
not on the r_i's.

Notice that this zero-testing element does not allow to decrypt, just to test
whether the m_i are all equal to 0 or not.


WHAT IS IMPLEMENTED?
--------------------

We provide a proof-of-concept implementation of the multilinear maps scheme
running on a N-multipartite Diffie-Hellman key exchange protocol.

The implementation is done using C++ with the GMP library, available at
http://gmplib.org

Modify the Makefile accordingly and type:

```
$ make
$ ./multimap
```

The program demonstrates how a setup phase is performed, and simulates the
view of N users. The parameters can me modified in Multimap.h and main.cpp.
Notice the following parameters:

in main.cpp

```#define LINUX```
(if you run on Linux, uses CLOCK_MONOTONIC)

```#define VERBOSE```
(if you want to verbose the code)

in Multimap.h

```#define INSTANTIATION```

According to the value of INSTANTIATION, different parameters set are used.

/!\ WARNING: due to the design of the code, the CRT coefficients are stored
during the program execution. Therefore the parameters 3 and 4 require a large
amount of RAM (around 100GB for INSTANTIATION 4). The program is easy to
modify to recompute these CRT coefficients each time but this induces a
significant computational overhead.

NB: Instead of takings the p_i's as prime numbers the program generate them as
products of (not so small) primes; therefore one should check whether the
security claimed from the parameters with the p_i's primes still hold in that
case.

This simplification is due to the function mpz_nextprime() used in the program
(which is clearly not optimal). Besides, in a similar approach in a real
application, one should ensure that x0 is hard to factorize.

NB 2: Note that the Setup party is only to be run once by a (trusted) third
party. One might consider it to be accessible as a common (public) knowledge
by all the parties; therefore this step can be forgotten when considering the
performances of the multilinear maps.

LICENSE
-------

This code is licensed under Creative Commons, Attribution, NonCommercial, 
NoDerivs http://creativecommons.org/licenses/by-nc-nd/3.0/



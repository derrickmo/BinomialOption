#ifndef BinomialMethod_hpp
#define BinomialMethod_hpp

#include "lattice.cpp"
#include "BinomialLatticeStrategy.hpp"
#include <cmath>

#include <iostream>
using namespace std;

class BinomialMethod
{
private:
	// Underlying data structure
	Lattice<double, int, 2> lattice; // Magic number == 2 means binomial
	BinomialLatticeStrategy *str;	 // Pointer to an algorithm

	double disc;

public:
	// Default constructor
	BinomialMethod();

	// Constructor taking discount factor, strategy (e.g. CRR) and number of steps
	BinomialMethod(double discounting, BinomialLatticeStrategy &strategy, int N);

	// Initialise lattice data structure
	void buildLattice(int N);

	// Initialise lattice node values (Forward Induction)
	void modifyLattice(double U);

	// Calculate derivative price (Backward Induction)
	double getPrice(const Vector<double, int> &RHS);

	// Handy function to give us the size at expiry date
	Vector<double, int> BasePyramidVector() const;

	// Underlying lattice
	const Lattice<double, int, 2> &getLattice() const;
};

#endif

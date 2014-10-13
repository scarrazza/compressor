// $Id: randomgenerator.h 636 2013-03-25 16:58:08Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

/**
 * @brief The RandomGenerator class
 * This class manages the random number generators
 */
class RandomGenerator
{
private:
  gsl_rng *fR;

public:
  RandomGenerator(int,unsigned long int);
  ~RandomGenerator();

  unsigned long int GetRandomInt();
  unsigned long int GetRandomUniform(unsigned long int);
  double GetRandomUniform();
  double GetRandomUniform(double, double);
  double GetRandomUniformPos();
  double GetRandomGausDev(const double);

  // Set Methods
  void SetSeed(unsigned long int);
};

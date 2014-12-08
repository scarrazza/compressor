#pragma once

#include <vector>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
using std::vector;

#define VERSION "1.0.0"

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

/**
 * @brief randomize
 * @param rg
 * @param index
 */
void randomize(int max, RandomGenerator *rg, vector<int> &index);

double ERF(int f, int nx, double **x, double **xavg);
double ERFS(int f, int nx, int reg, double ***x, double ***xavg);
double ERFC(int size, double *x, double *xavg);

void ComputeCV(vector<double> x, double& cv, double &md,
	       double &dn50,double &up50,
	       double &dn68,double &up68,
	       double &dn90,double &up90);


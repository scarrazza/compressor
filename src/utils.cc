#include <iostream>
#include <cmath>
#include "utils.hh"
#include <algorithm>
using namespace std;

// See GSL documentation
const gsl_rng_type* rngs[] = { gsl_rng_ranlux,
                               gsl_rng_cmrg,
                               gsl_rng_mrg,
                               gsl_rng_mt19937,
                               gsl_rng_gfsr4,
                               gsl_rng_ran0,
                               gsl_rng_ran1,
                               gsl_rng_ran2,
                               gsl_rng_ran3,
                               gsl_rng_rand,
                               gsl_rng_ranlxd1,
                               gsl_rng_ranlxd2,
                               gsl_rng_ranlxs0,
                               gsl_rng_ranlxs1,
                               gsl_rng_ranlxs2,
                               gsl_rng_taus,
                               gsl_rng_taus2
                             };

/**
 * @brief RandomGenerator::RandomGenerator
 */
RandomGenerator::RandomGenerator(int method, unsigned long int seed)
{
  // just export GSL_RNG_TYPE=mrg to change the generator
  //gsl_rng_env_setup();
  const gsl_rng_type *T = rngs[method];

  fR = gsl_rng_alloc(T);
  cout << "* Random Generator allocated: " << gsl_rng_name(fR) << endl;
  gsl_rng_set(fR, seed);
}

/**
 * @brief RandomGenerator::~RandomGenerator
 */
RandomGenerator::~RandomGenerator()
{
  gsl_rng_free(fR);
}

/**
 * @brief RandomGenerator::SetSeed
 * @param s
 */
void RandomGenerator::SetSeed(unsigned long int s)
{
  gsl_rng_set(fR, s);
}

/**
 * @brief RandomGenerator::GetRandomInt
 * @return a random integer from the generator.
 * Min max are equally likely and depends on the
 * algorithm used.
 */
unsigned long int RandomGenerator::GetRandomInt()
{
  return gsl_rng_get(fR);
}

/**
 * @brief RandomGenerator::GetRandomUniform
 * @param n
 * @return an integer from 0 to n-1
 */
unsigned long int RandomGenerator::GetRandomUniform(unsigned long int n)
{
  return gsl_rng_uniform_int(fR, n);
}

/**
 * @brief RandomGenerator::GetRandomUniform
 * @return a double from [0,1)
 */
double RandomGenerator::GetRandomUniform()
{
  return gsl_rng_uniform(fR);
}

/**
 * @brief RandomGenerator::GetRandomUniform
 * @return a double from [a,b)
 */
double RandomGenerator::GetRandomUniform(double a, double b)
{
  return (b-a)*gsl_rng_uniform(fR) + a;
}

/**
 * @brief RandomGenerator::GetRandomUniformPos
 * @return a double from (0,1)
 */
double RandomGenerator::GetRandomUniformPos()
{
  return gsl_rng_uniform_pos(fR);
}

/**
 * @brief Gaussian Random number generator
 * using Box-Muller algorithm
 * @param x
 * @return
 */
double RandomGenerator::GetRandomGausDev(const double sigma)
{
  return gsl_ran_gaussian(fR, sigma);
}

void randomize(int max, RandomGenerator *rg, vector<int> &index)
{
  const size_t rep = index.size();
  for (size_t i = 0; i < rep; i++)
    {
      bool pass = false; int id;
      while (pass == false)
        {
          pass = true;
          id = rg->GetRandomUniform(max);
          for (size_t j = 0; j < rep; j++)
            if (index[j] == id) pass = false;
        }
      index[i] = id;
    }
}

double ERF(int f, int nx, double **x, double **xavg)
{
  double res = 0;
  for (int r = 0; r < f; r++)
    for (int ix = 0; ix < nx; ix++)
      {
        if (xavg[r][ix] != 0)
          res += pow( (x[r][ix]-xavg[r][ix])/xavg[r][ix], 2.0);
      }
  return res;
}

double ERFS(int f, int nx, int reg, double ***x, double ***xavg)
{
  double res = 0;
  for (int r = 0; r < f; r++)
    for (int ix = 0; ix < nx; ix++)
      for (int l = 0; l < reg; l++)
        {
          if (xavg[r][ix][l] != 0)
            res += pow( (x[r][ix][l]-xavg[r][ix][l])/xavg[r][ix][l], 2.0);
        }
  return res;
}

void ComputeCV(vector<double> x, double& cv, double &md,
	       double &dn50,double &up50,
	       double &dn68,double &up68,
	       double &dn90,double &up90)
{
  size_t rep = x.size();
  double sum = 0;
  for (size_t i = 0; i < rep; i++) sum += x[i];

  cv = sum / (double) rep;

  std::sort(x.begin(),x.end());

  int esc = (int) (rep*(1-0.50)/2);
  up50 = x[rep-1-esc];
  dn50 = x[esc];

  esc = (int) (rep*(1-0.68)/2);
  up68 = x[rep-1-esc];
  dn68 = x[esc];

  esc = (int) (rep*(1-0.90)/2);
  up90 = x[rep-1-esc];
  dn90 = x[esc];

  if (rep % 2 == 0)
    md = (x[rep / 2 - 1] + x[rep / 2]) / 2.0;
  else
    md = x[rep/ 2];
}


#pragma once

#include <vector>
#include "LHAPDF/LHAPDF.h"
using std::vector;

class EstimatorsM;
class EstimatorsS;
class EstimatorsC;
class Grid;
class RandomGenerator;

class Minimizer
{
public:
  Minimizer(vector<LHAPDF::PDF*> const& pdf, Grid* const& x, double const& Q);
  ~Minimizer();
  vector<EstimatorsM*> GetMomentEstimators()  const  { return _estM; }
  vector<EstimatorsS*> GetStatEstimators()  const  { return _estS;   }
  vector<EstimatorsC*> GetCorrEstimators()  const  { return _estC;   }

  vector<double**>  GetPriorMomentEstValues() const { return _estMval; }
  vector<double***> GetPriorStatEstValues()   const { return _estSval; }
  vector<double*>   GetPriorCorrEstValues()   const { return _estCval; }
  vector<int> GetIDS() const { return _ids; }
  double iterate();
  void setupminimizer(int rep, vector<double> N, RandomGenerator *rg);
  vector<int> getIndex() const { return _index; }


private:
  int _nf;  
  double _Q;
  Grid *_x;
  vector<int> _ids;
  vector<int> _index;
  vector<double> _N;
  vector<EstimatorsM*> _estM;
  vector<EstimatorsS*> _estS;
  vector<EstimatorsC*> _estC;
  vector<double**> _estMval;
  double** _iteMval;
  vector<double***> _estSval;
  double*** _iteSval;
  vector<double*> _estCval;  
  double* _iteCval;
  vector<LHAPDF::PDF*> _pdf;
  int _NX;
  int _rep;
  int _nmut;
  vector< vector<int> > _mut;
  RandomGenerator *_rg;
};

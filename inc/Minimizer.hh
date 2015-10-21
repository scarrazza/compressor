// Compressor - January 2015
// Author:  Stefano Carrazza
// Contact: stefano.carrazza@mi.infn.it
#pragma once

#include <vector>
#include "LHAPDF/LHAPDF.h"
#include "TMatrixD.h"
using std::vector;

class EstimatorsM;
class EstimatorsS;
class EstimatorsC;
class Grid;
class RandomGenerator;
class LocalPDF;

class Minimizer
{
public:
  Minimizer(LocalPDF* const& pdf, Grid* const& x, int const& nf);
  ~Minimizer();
  vector<EstimatorsM*> GetMomentEstimators()const  { return _estM; }
  vector<EstimatorsS*> GetStatEstimators()  const  { return _estS;   }
  vector<EstimatorsC*> GetCorrEstimators()  const  { return _estC;   }

  vector<double**>  GetPriorMomentEstValues() const { return _estMval; }
  vector<double***> GetPriorStatEstValues()   const { return _estSval; }
  vector<double*>   GetPriorCorrEstValues()   const { return _estCval; }
  TMatrixD          GetPriorInvMatrix()       const { return _invmatrix; }
  vector<int> GetIDS() const { return _ids; }
  double iterate();
  double iterate_w();
  void setupminimizer(int rep, vector<double> N, RandomGenerator *rg);
  vector<int> getIndex() const { return _index; }
  vector<double> getW() const { return _w; }


private:
  int _nf;  
  Grid *_x;
  vector<int> _ids;
  vector<int> _index;
  vector<double> _N;
  vector<EstimatorsM*> _estM;
  vector<EstimatorsS*> _estS;
  vector<EstimatorsC*> _estC;
  TMatrixD _invmatrix;
  vector<double**> _estMval;
  double** _iteMval;
  vector<double***> _estSval;
  double*** _iteSval;
  vector<double*> _estCval;  
  double* _iteCval;
  LocalPDF *_pdf;
  int _rep;
  int _nmut;
  int _Nx;
  vector< vector<int> > _mut;
  RandomGenerator *_rg;
  vector<double> _w;
  vector< vector<double> > _mutw;
};

// Compressor - January 2015
// Author:  Stefano Carrazza
// Contact: stefano.carrazza@mi.infn.it

#include "Minimizer.hh"
#include "Estimators.hh"
#include "Grid.hh"
#include "utils.hh"
#include "TMatrixD.h"
#include <iostream>
using namespace std;

Minimizer::Minimizer(LocalPDF* const& pdf,
                     Grid* const& x, int const& nf):
  _nf(nf),
  _x(x),
  _pdf(pdf),
  _Nx(x->size())
{
  // computing estimators
  _estM.push_back(new CentralValue());
  _estM.push_back(new StdDeviation());
  _estM.push_back(new Skewness());
  _estM.push_back(new Kurtosis());
  //_estM.push_back(new moment5th());
  //_estM.push_back(new moment6th());
  _estS.push_back(new Kolmogorov());
  _estC.push_back(new Correlation(2*_nf+1));  

  _w.resize(_pdf->size(),1.0);

  // extracting active flavors
  _ids.resize(2*_nf+1);
  for (int i = -_nf; i <= _nf; i++) _ids[i+_nf] = i;

  // fill evaluation vector
  for (size_t es = 0; es < _estM.size(); es++)
    {
      _estMval.push_back(new double*[_ids.size()]);
      for (size_t fl = 0; fl < _ids.size(); fl++)
        _estMval[es][fl] = new double[_x->size()];
    }

  // default index distribution
  _index.resize(_pdf->size()-1);
  for (size_t r = 0; r < _pdf->size()-1; r++)
    _index[r] = r+1;

  // computing estimators
  for (size_t fl = 0; fl < _ids.size(); fl++)
    for (int ix = 0; ix < _x->size(); ix++)
      for (size_t es = 0; es < _estM.size(); es++)
        _estMval[es][fl][ix] = _estM[es]->Evaluate(_pdf,_ids[fl],_index,_w,ix);

  // computing statistical tests
  for (size_t es = 0; es < _estS.size(); es++)
    {
      _estSval.push_back(new double**[_ids.size()]);
      for (size_t fl = 0; fl < _ids.size(); fl++)
        {
          _estSval[es][fl] = new double*[_x->size()];
          for (int ix = 0; ix < x->size(); ix++)
            _estSval[es][fl][ix] = new double[_estS[es]->getRegions()];
        }
    }

  // filling estimatorS
  for (size_t fl = 0; fl < _ids.size(); fl++)
    for (int ix = 0; ix < _x->size(); ix++)
      for (size_t es = 0; es < _estS.size(); es++)
        {
          vector<double> res = _estS[es]->Evaluate(_pdf,_ids[fl],_index,_w,ix);
          for (int l = 0; l < _estS[es]->getRegions(); l++) _estSval[es][fl][ix][l] = res[l];
        }  

  // filling estimatorsC
  for (size_t i = 0; i < _estC.size(); i++)
    {
      _estCval.push_back(new double[_estC[i]->getSize()]);
      TMatrixD m = _estC[i]->Evaluate(_pdf,_ids,_index,_w,_x);
      _invmatrix.ResizeTo(m); _invmatrix = m; _invmatrix.Invert();

      TMatrixD r = m*_invmatrix;
      for (int j = 0; j < _estC[i]->getSize(); j++) _estCval[i][j] = 0;
      for (int j = 0; j < _estC[i]->getSize(); j++) _estCval[i][0] += r(j,j);
    }

  // Preparing iteration containers
  _iteMval = new double*[_ids.size()];
  for (size_t fl = 0; fl < _ids.size(); fl++)
    _iteMval[fl] = new double[_x->size()];

  _iteSval = new double**[_ids.size()];
  for (size_t fl = 0; fl < _ids.size(); fl++)
    {
      _iteSval[fl] = new double*[_x->size()];
      for (int ix = 0; ix < x->size(); ix++)
        _iteSval[fl][ix] = new double[_estS[0]->getRegions()];
    }

  _iteCval = new double[_estC[0]->getSize()];
}

Minimizer::~Minimizer()
{
  for (size_t i = 0; i < _estMval.size(); i++) {
    for (size_t j = 0; j < _ids.size(); j++)
      if (_estMval[i][j]) delete[] _estMval[i][j];
    if (_estMval[i]) delete[] _estMval[i];
    }
  _estMval.clear();

  for (size_t i = 0; i < _estSval.size(); i++) {
    for (size_t j = 0; j < _ids.size(); j++) {
      for (int z = 0; z < _Nx; z++)
        if (_estSval[i][j][z]) delete[] _estSval[i][j][z];
      if (_estSval[i][j]) delete[] _estSval[i][j];
      }
    if (_estSval[i]) delete[] _estSval[i];
    }
  _estSval.clear();

  for (size_t i = 0; i < _estCval.size(); i++)
    if (_estCval[i]) delete[] _estCval[i];
  _estCval.clear();

  for (size_t i = 0; i < _ids.size(); i++)
    if (_iteMval[i]) delete[] _iteMval[i];
  delete[] _iteMval;

  for (size_t i = 0; i < _ids.size(); i++) {
    for (int j = 0; j < _Nx; j++)
      if (_iteSval[i][j]) delete[] _iteSval[i][j];
    if (_iteSval[i]) delete[] _iteSval[i];
    }

  delete[] _iteSval;
  delete[] _iteCval;
  _index.clear();
  _ids.clear();
  _mut.clear();
  _w.clear();

}

void Minimizer::setupminimizer(int rep, vector<double> N, RandomGenerator *rg)
{
  _N = N;
  _index.resize(rep); for (int i = 0; i < rep; i++) _index[i] = i+1;
  _rep = rep;
  _nmut = 5;
  _mut.resize(_nmut); for (int i = 0; i < _nmut; i++) _mut[i].resize(_rep,0);
  _rg = rg;
  _w.resize(_rep, 1.0);
  _mutw.resize(_nmut); for (int i = 0; i < _nmut; i++) _mutw[i].resize(_rep,1);
}

double Minimizer::iterate()
{
  //const size_t Msize = _estM.size()-2;
  const size_t Msize = _estM.size();
  double berf = 0;  
  for (size_t es = 0; es < Msize; es++)
    {
      for (size_t fl = 0; fl <_ids.size(); fl++)
        for (int ix = 0; ix < _x->size(); ix++)
          _iteMval[fl][ix] = _estM[es]->Evaluate(_pdf,_ids[fl],_index,_w,ix);
      berf += ERF(_ids.size(), _x->size(), _iteMval, _estMval[es]) / _N[es];
    }

  for (size_t es = 0; es < _estS.size(); es++)
    {
      for (size_t fl = 0; fl <_ids.size(); fl++)
        for (int ix = 0; ix < _x->size(); ix++)
          {
            vector<double> res = _estS[es]->Evaluate(_pdf,_ids[fl],_index,_w,ix);
            for (int l = 0; l < _estS[es]->getRegions(); l++) _iteSval[fl][ix][l] = res[l];
          }
      berf += ERFS(_ids.size(), _x->size(), _estS[es]->getRegions(), _iteSval, _estSval[es]) / _N[es+Msize];
    }
  for (size_t es = 0; es < _estC.size(); es++)
    {
      TMatrixD m = _estC[es]->Evaluate(_pdf,_ids,_index,_w,_x);
      TMatrixD r = m*_invmatrix;

      _iteCval[0] = 0;
      for (int l = 0; l < _estC[es]->getSize(); l++) _iteCval[0] += r(l,l);

      berf += ERFC(1, _iteCval, _estCval[es]) / _N[es+Msize+_estS.size()];      
    }    


  // set mut
  for (int i = 0; i < _nmut; i++)
    for (int j = 0; j < _rep; j++)
      {
        _mut[i][j] = _index[j];
        _mutw[i][j] = _w[j];
      }

  // GA mutation
  for (int i = 0; i < _nmut; i++)
    {
      const double g = _rg->GetRandomUniform();

      int nmut = 4; // 10%
      if (g <= 0.3) nmut = 1; // 30%
      else if (g > 0.3 && g <= 0.6) nmut = 2; // 30%
      else if (g > 0.6 && g <= 0.7) nmut = 3; // 10%

      for (int t = 0; t < nmut; t++)
        {
          const int pos = _rg->GetRandomUniform(_pdf->size()-1)+1;
          _mut[i][_rg->GetRandomUniform(_rep)] = pos;
        }

      // GA for weights -- keep 1 for first mutant
      double wtot = 0;
      for (int t = 0; t < _w.size(); t++)
        {
          _mutw[i][t] = -log(_rg->GetRandomUniform(0,1));
          wtot += _mutw[i][t];
        }
      for (int t = 0; t < _w.size(); t++)
        _mutw[i][t] /= wtot/_rep;
    }

  // Compute ERF
  double *erf = new double[_nmut];
  for (int i = 0; i < _nmut; i++)
    {
      erf[i] = 0.0;
      for (size_t es = 0; es < Msize; es++)
        {
          for (size_t fl = 0; fl <_ids.size(); fl++)
            for (int ix = 0; ix < _x->size(); ix++)
              _iteMval[fl][ix] = _estM[es]->Evaluate(_pdf,_ids[fl],_mut[i],_mutw[i],ix);
          erf[i] += ERF(_ids.size(), _x->size(), _iteMval, _estMval[es]) / _N[es];
        }

      for (size_t es = 0; es < _estS.size(); es++)
        {
          for (size_t fl = 0; fl <_ids.size(); fl++)
            for (int ix = 0; ix < _x->size(); ix++)
              {
                vector<double> res = _estS[es]->Evaluate(_pdf,_ids[fl],_mut[i],_mutw[i],ix);
                for (int l = 0; l < _estS[es]->getRegions(); l++) _iteSval[fl][ix][l] = res[l];
              }
          erf[i] += ERFS(_ids.size(), _x->size(), _estS[es]->getRegions(), _iteSval, _estSval[es]) / _N[es+Msize];
        }

      for (size_t es = 0; es < _estC.size(); es++)
        {
          TMatrixD m = _estC[es]->Evaluate(_pdf,_ids,_mut[i],_mutw[i],_x);
          TMatrixD r = m*_invmatrix;

	  _iteCval[0] = 0;
          for (int l = 0; l < _estC[es]->getSize(); l++) _iteCval[0] += r(l,l);

          erf[i] += ERFC(1, _iteCval, _estCval[es]) / _N[es+Msize+_estS.size()];
        }

    }

  // Selection
  int id = 0;
  double bestchi2 = erf[0];
  for (int i = 0; i < _nmut; i++)
    if (erf[i] < bestchi2)
      {
        bestchi2 = erf[i];
        id = i;
      }

  if (bestchi2 < berf)
    for (int i = 0; i < (int) _rep; i++) { _index[i] = _mut[id][i]; _w[i] = _mutw[id][i]; }
  else bestchi2 = berf;

  delete[] erf;

  return bestchi2;
}

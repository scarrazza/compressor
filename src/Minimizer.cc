#include "Minimizer.hh"
#include "Estimators.hh"
#include "Grid.hh"
#include "utils.hh"
#include <iostream>
using namespace std;

Minimizer::Minimizer(vector<LHAPDF::PDF*> const& pdf,
                     Grid* const& x, double const& Q):
  _nf(3),
  _Q(Q),
  _x(x),
  _pdf(pdf)
{
  // computing estimators
  _estM.push_back(new CentralValue());
  _estM.push_back(new StdDeviation());
  _estM.push_back(new Skewness());
  _estM.push_back(new Kurtosis());
  _estM.push_back(new moment5th());
  _estM.push_back(new moment6th());

  _estS.push_back(new Kolmogorov());

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
  _index.resize(_pdf.size()-1);
  for (size_t r = 0; r < _pdf.size()-1; r++)
    _index[r] = r+1;

  // computing estimators
  for (size_t fl = 0; fl < _ids.size(); fl++)
    for (int ix = 0; ix < _x->size(); ix++)
      for (size_t es = 0; es < _estM.size(); es++)
        _estMval[es][fl][ix] = _estM[es]->Evaluate(_pdf,_ids[fl],_index,_x->at(ix),_Q);

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
          vector<double> res = _estS[es]->Evaluate(_pdf,_ids[fl],_index,_x->at(ix),_Q);
          for (int l = 0; l < _estS[es]->getRegions(); l++) _estSval[es][fl][ix][l] = res[l];
        }

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
}

void Minimizer::setupminimizer(int rep, vector<double> N, RandomGenerator *rg)
{
  _N = N;
  _index.resize(rep); for (int i = 0; i < rep; i++) _index[i] = i+1;
  _rep = rep;
  _nmut = 5;
  _mut.resize(_nmut); for (int i = 0; i < _nmut; i++) _mut[i].resize(_rep,0);
  _rg = rg;
}

double Minimizer::iterate()
{
  const size_t Msize = _estM.size()-2;
  double berf = 0;
  for (size_t es = 0; es < Msize; es++)
    {
      for (size_t fl = 0; fl <_ids.size(); fl++)
        for (int ix = 0; ix < _x->size(); ix++)
          _iteMval[fl][ix] = _estM[es]->Evaluate(_pdf,_ids[fl],_index,_x->at(ix),_Q);
      berf += ERF(_ids.size(), _x->size(), _iteMval, _estMval[es]) / _N[es];
    }

  for (size_t es = 0; es < _estS.size(); es++)
    {
      for (size_t fl = 0; fl <_ids.size(); fl++)
        for (int ix = 0; ix < _x->size(); ix++)
          {
            vector<double> res = _estS[es]->Evaluate(_pdf,_ids[fl],_index,_x->at(ix),_Q);
            for (int l = 0; l < _estS[es]->getRegions(); l++) _iteSval[fl][ix][l] = res[l];
          }
      berf += ERFS(_ids.size(), _x->size(), _estS[es]->getRegions(), _iteSval, _estSval[es]) / _N[es+Msize];
    }

  // set mut
  for (int i = 0; i < _nmut; i++)
    for (int j = 0; j < _rep; j++)
        _mut[i][j] = _index[j];

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
          const int pos = _rg->GetRandomUniform(_pdf.size()-1)+1;
          _mut[i][_rg->GetRandomUniform(_rep)] = pos;
        }
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
              _iteMval[fl][ix] = _estM[es]->Evaluate(_pdf,_ids[fl],_mut[i],_x->at(ix),_Q);
          erf[i] += ERF(_ids.size(), _x->size(), _iteMval, _estMval[es]) / _N[es];
        }

      for (size_t es = 0; es < _estS.size(); es++)
        {
          for (size_t fl = 0; fl <_ids.size(); fl++)
            for (int ix = 0; ix < _x->size(); ix++)
              {
                vector<double> res = _estS[es]->Evaluate(_pdf,_ids[fl],_mut[i],_x->at(ix),_Q);
                for (int l = 0; l < _estS[es]->getRegions(); l++) _iteSval[fl][ix][l] = res[l];
              }
          erf[i] += ERFS(_ids.size(), _x->size(), _estS[es]->getRegions(), _iteSval, _estSval[es]) / _N[es+Msize];
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
    for (int i = 0; i < (int) _rep; i++) _index[i] = _mut[id][i];
  else bestchi2 = berf;

  delete[] erf;

  return bestchi2;
}

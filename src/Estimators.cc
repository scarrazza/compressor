
#include <iostream>
#include "Estimators.hh"
using namespace std;

double CentralValue::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  double res = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    res += pdf[index[i]]->xfxQ(fl,x,Q);
  return res / n;
}

double StdDeviation::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  double sum = 0, sq_sum = 0;

  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q);
      sum += v;
      sq_sum += v*v;
    }

  return sqrt(sq_sum / n - sum/n * sum/n);
}

double Skewness::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v;
    }

  return (sum / n) / (st*st*st);
}

double Kurtosis::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v*v;
    }

  return (sum / n) / (st*st*st*st);
}

double moment5th::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v*v*v;
    }

  return (sum / n) / (st*st*st*st*st);
}

double moment6th::Evaluate(const vector<LHAPDF::PDF *> &pdf, const int &fl,
                              const vector<int> &index, const double &x, const double &Q) const
{
  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  double sum = 0;
  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q) - cv;
      sum += v*v*v*v*v*v;
    }

  return (sum / n) / (st*st*st*st*st*st);
}

vector<double> Kolmogorov::Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                        vector<int> const& index,double const& x, double const& Q) const
{
  vector<double> res(_regions,0);

  CentralValue ecv;
  StdDeviation est;

  double cv = ecv.Evaluate(pdf,fl,index,x,Q);
  double st = est.Evaluate(pdf,fl,index,x,Q);

  const int n = index.size();
  for (int i = 0; i < n; i++)
    {
      const double v = pdf[index[i]]->xfxQ(fl,x,Q);
      if (v <= cv -2*st)
        res[0] += 1;
      else if (v <= cv - st)
        res[1] += 1;
      else if (v <= cv)
        res[2] += 1;
      else if (v <= cv + st)
        res[3] += 1;
      else if (v <= cv + 2*st)
        res[4] += 1;
      else if (v > cv + 2*st)
        res[5] += 1;
      else
        cout << "ERROR" << endl;
    }

  for (int l = 0; l < _regions; l++) res[l] /= (double) n;

  return res;
}

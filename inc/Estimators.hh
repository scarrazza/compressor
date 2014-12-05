#pragma once

#include <vector>
#include <string>
using std::vector;
using std::string;

#include "LHAPDF/LHAPDF.h"

class EstimatorsM
{
protected:
  string _name;
public:
  EstimatorsM(string name): _name(name) {}  
  string  getName() const { return _name; }
  virtual double Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                          vector<int> const& index,double const& x, double const& Q) const = 0;
};

class EstimatorsS
{
protected:
  int _regions;
  string _name;
public:
  EstimatorsS(string name, int r): _regions(r), _name(name) {}
  string  getName() const { return _name; }
  int     getRegions() const { return _regions; }
  virtual vector<double> Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                          vector<int> const& index,double const& x, double const& Q) const = 0;
};

class Kolmogorov: public EstimatorsS
{
public:
  Kolmogorov(): EstimatorsS("Kolmogorov",6) {}
  vector<double> Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                          vector<int> const& index,double const& x, double const& Q) const;
};

class CentralValue: public EstimatorsM
{
public:
  CentralValue(): EstimatorsM("Central Value") {}
  double Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                  vector<int> const& index,double const& x, double const& Q) const;
};

class StdDeviation: public EstimatorsM
{
public:
  StdDeviation(): EstimatorsM("StdDeviation") {}
  double Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                  vector<int> const& index,double const& x, double const& Q) const;
};

class Skewness: public EstimatorsM
{
public:
  Skewness(): EstimatorsM("Skewness") {}
  double Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                  vector<int> const& index,double const& x, double const& Q) const;
};

class Kurtosis: public EstimatorsM
{
public:
  Kurtosis(): EstimatorsM("Kurtosis") {}
  double Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                  vector<int> const& index,double const& x, double const& Q) const;
};

class moment5th: public EstimatorsM
{
public:
  moment5th(): EstimatorsM("5th moment") {}
  double Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                  vector<int> const& index,double const& x, double const& Q) const;
};

class moment6th: public EstimatorsM
{
public:
  moment6th(): EstimatorsM("6th moment") {}
  double Evaluate(vector<LHAPDF::PDF*> const& pdf, int const& fl,
                  vector<int> const& index,double const& x, double const& Q) const;
};

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "randomgenerator.h"
#include "mini.h"
using namespace std;

void ComputeCV(int rep, double *x, double *cv, double *std)
{
  double sum = 0, sq_sum = 0;
  for (size_t i = 0; i < rep; i++)
    {
      const double v = x[i];
      sum += v;
      sq_sum += v*v; 
    }

  *cv = sum / (double) rep;
  *std= sqrt(sq_sum / (double) rep - *cv * *cv);
}

int main(int argc, char **argv)
{
  string filename;
  int rep = 100;
  int trials = 0;
  if (argc > 3)
    {
      rep = atoi(argv[1]);
      filename.assign(argv[3]);
      trials = atoi(argv[2]);
    }
  else
    {
      cout << "\nusage: analyser [desired replicas] [trials] [PDF set name]\n" << endl;
      exit(-1);
    }

  cout << "Working with " << rep << " replicas, with " << trials << " trials" << endl; 

  // load PDF set
  const LHAPDF::PDFSet set(filename.c_str());
  vector<LHAPDF::PDF*> pdf = set.mkPDFs();
  vector<int> fPids = pdf[0]->flavors();

  RandomGenerator *rg = new RandomGenerator(0,rep);
  Mini *min = new Mini(rep,pdf,0);
  double *erfcv = new double[trials];
  double *erfsd = new double[trials];
  double *erfsk = new double[trials];
  double *erfku = new double[trials];
  double *erfko = new double[trials];
  int* index = new int[rep];

  for (int t = 0; t < trials; t++)
    {      
      for (int i = 0; i < rep; i++) index[i] = 0;
      for (int i = 0; i < rep; i++) 
	{
	  bool pass = false;
	  int id;
	  while (pass == false)
	    {
	      pass = true;
	      id = rg->GetRandomUniform(pdf.size()-1)+1;	      
	      for (int j = 0; j < rep; j++)
		if (index[j] == id ) pass = false;	      
	    }
	  index[i] = id;
	}

      double ecv = 0, estd = 0, esk = 0, eku = 0, eko = 0;
      double *res = new double[6];
      for (int f = 0; f < (int) fPids.size(); f++)
	for (int i = 0; i < (int) min->GetX().size(); i++)
	  {
	    double cv = 0, std = 0, sk = 0, kur = 0;
	    min->ComputeEstimators(rep, min->GetX()[i], 1.0, fPids[f], index, 
				   &cv, &std, &sk, &kur, res);	
	    ecv  += pow(min->GetCV(f,i) - cv, 2.0);
	    estd += pow(min->GetSD(f,i) - std, 2.0);
	    esk  += pow(min->GetSK(f,i) - sk, 2.0);
	    eku  += pow(min->GetKU(f,i) - kur, 2.0);
	    
	    for (int l = 0; l < 6; l++)
	      eko += pow(min->GetKO(f,i,l) - res[l], 2.0);	      
	  }

      erfcv[t] = ecv;
      erfsd[t] = estd;
      erfsk[t] = esk;
      erfku[t] = eku;
      erfko[t] = eko;

      delete[] res;
    }

  cout << scientific;

  fstream f;
  
  double cv = 0, std = 0;
  ComputeCV(trials, erfcv, &cv, &std);
  cout << "CV:  " << cv << "\t" << std << endl;

  f.open("cv_set_r.dat", ios::out|ios::app);
  f << fixed << rep << scientific << "\t" << cv << "\t" << std << endl;
  f.close();

  ComputeCV(trials, erfsd, &cv, &std);
  cout << "STD: " << cv << "\t" << std << endl;

  f.open("std_set_r.dat", ios::out|ios::app);
  f << fixed << rep << scientific << "\t" << cv << "\t" << std << endl;
  f.close();

  ComputeCV(trials, erfsk, &cv, &std);
  cout << "SKE: " << cv << "\t" << std << endl;

  f.open("ske_set_r.dat", ios::out|ios::app);
  f << fixed << rep << scientific << "\t" << cv << "\t" << std << endl;
  f.close();

  ComputeCV(trials, erfku, &cv, &std);
  cout << "KUR: " << cv << "\t" << std << endl;

  f.open("kur_set_r.dat", ios::out|ios::app);
  f << fixed << rep << scientific << "\t" << cv << "\t" << std << endl;
  f.close();

  ComputeCV(trials, erfko, &cv, &std);
  cout << "KOL: " << cv << "\t" << std << endl;

  f.open("kol_set_r.dat", ios::out|ios::app);
  f << fixed << rep << scientific << "\t" << cv << "\t" << std << endl;
  f.close();

  delete[] erfcv;
  delete[] erfsd;
  delete[] erfsk;
  delete[] erfku;
  delete[] erfko;
  delete rg;
  delete min;
  delete[] index;

  return 0;
}


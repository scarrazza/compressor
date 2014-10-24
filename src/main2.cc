#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "randomgenerator.h"
#include "mini.h"
using namespace std;

void ComputeCV(int rep, double *x, double& cv, double &md, 
	       double &dn50,double &up50,
	       double &dn68,double &up68,
	       double &dn90,double &up90)
{
  vector<double> xval;  
  double sum = 0, sq_sum = 0;
  for (size_t i = 0; i < rep; i++) {
    const double v = x[i];
    sum += v;
    xval.push_back(v);
  }

  cv = sum / (double) rep;

  std::sort(xval.begin(),xval.end());
  
  int esc = (int) (rep*(1-0.50)/2);  
  up50 = xval[xval.size()-1-esc];
  dn50 = xval[esc];   
  
  esc = (int) (rep*(1-0.68)/2);
  up68 = xval[xval.size()-1-esc];
  dn68 = xval[esc];

  esc = (int) (rep*(1-0.90)/2);
  up90 = xval[xval.size()-1-esc];
  dn90 = xval[esc];

  size_t size = xval.size();
  if (size  % 2 == 0)
    md = (xval[size / 2 - 1] + xval[size / 2]) / 2.0;
  else 
    md = xval[size / 2];  

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
  double *erfcv2 = new double[trials];
  double *erfsd2 = new double[trials];
  double *erfsk2 = new double[trials];
  double *erfku2 = new double[trials];
  double *erfko2 = new double[trials];
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
      double ecv2 = 0, estd2 = 0, esk2 = 0, eku2 = 0, eko2 = 0;
      double *res = new double[6];

      for (int f = 0; f < (int) fPids.size(); f++)
	for (int i = 0; i < (int) min->GetX().size(); i++)
	  {
	    double cv = 0, std = 0, sk = 0, kur = 0;
	    min->ComputeEstimators(rep, min->GetX()[i], 1.0, fPids[f], index, cv, std, sk, kur, res);	
	    ecv  += pow(min->GetCV(f,i) - cv, 2.0);
	    estd += pow(min->GetSD(f,i) - std, 2.0);
	    esk  += pow(min->GetSK(f,i) - sk, 2.0);
	    eku  += pow(min->GetKU(f,i) - kur, 2.0);
	    
	    for (int l = 0; l < 6; l++)
	      eko += pow(min->GetKO(f,i,l) - res[l], 2.0);	      
	    
	    // for relative difference
	    if (fPids[f] == 21 && (min->GetX()[i] < 0.6 && min->GetX()[i] > 1e-3))
	      {
		ecv2  += fabs( (min->GetCV(f,i) - cv)  / min->GetCV(f,i) );
		estd2 += fabs( (min->GetSD(f,i) - std) / min->GetSD(f,i) );
		esk2  += fabs( (min->GetSK(f,i) - sk)  / min->GetSK(f,i) );
		eku2  += fabs( (min->GetKU(f,i) - kur) / min->GetKU(f,i) );
		
		for (int l = 0; l < 6; l++)
		  {
		    if (min->GetKO(f,i,l) != 0)
		      eko2 += fabs( (min->GetKO(f,i,l) - res[l]) / min->GetKO(f,i,l) );	      
		  }
	      }
	  }

      erfcv[t] = ecv;
      erfsd[t] = estd;
      erfsk[t] = esk;
      erfku[t] = eku;
      erfko[t] = eko;

      erfcv2[t] = ecv2;
      erfsd2[t] = estd2;
      erfsk2[t] = esk2;
      erfku2[t] = eku2;
      erfko2[t] = eko2;

      delete[] res;
    }


  fstream f, g;
  f << scientific;
  g << scientific;
  cout << scientific;
  
  double cv = 0, md = 0, up50 = 0, dn50 = 0; 
  double up68 = 0, dn68 = 0, up90 = 0, dn90 = 0;
  
  /////////////////////////////////
  ComputeCV(trials, erfcv, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "CV:  " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
       
  f.open("cv_erf.dat", ios::out|ios::app);  
  f << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  f.close();

  ComputeCV(trials, erfcv2, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "CV2: " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;

  g.open("cv.dat", ios::out|ios::app);  
  g << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  g.close();
  /////////////////////////////////

  /////////////////////////////////
  ComputeCV(trials, erfsd, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "SD:  " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
       
  f.open("sd_erf.dat", ios::out|ios::app);  
  f << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  f.close();

  ComputeCV(trials, erfsd2, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "SD2: " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;

  g.open("sd.dat", ios::out|ios::app);  
  g << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  g.close();
  /////////////////////////////////

  /////////////////////////////////
  ComputeCV(trials, erfsk, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "SK:  " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
       
  f.open("sk_erf.dat", ios::out|ios::app);  
  f << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  f.close();

  ComputeCV(trials, erfsk2, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "SK2: " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;

  g.open("sk.dat", ios::out|ios::app);  
  g << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  g.close();
  /////////////////////////////////

  /////////////////////////////////
  ComputeCV(trials, erfku, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "KU:  " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
       
  f.open("ku_erf.dat", ios::out|ios::app);  
  f << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  f.close();

  ComputeCV(trials, erfku2, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "KU2: " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;

  g.open("ku.dat", ios::out|ios::app);  
  g << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  g.close();
  /////////////////////////////////

  /////////////////////////////////
  ComputeCV(trials, erfko, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "KO:  " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
       
  f.open("ko_erf.dat", ios::out|ios::app);  
  f << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  f.close();

  ComputeCV(trials, erfko2, cv, md, dn50, up50, dn68, up68, dn90, up90);
  cout << "KO2: " << cv << "\t" 
       << md << "\t" << dn50 << "\t" << up50 << "\t" 
       << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;

  g.open("ko.dat", ios::out|ios::app);  
  g << rep << "\t" << cv << "\t" << md << "\t" << dn50 << "\t" << up50 << "\t" 
    << dn68 << "\t" << up68 << "\t" << dn90 << "\t" << up90 << endl;
  g.close();
  /////////////////////////////////

  delete[] erfcv;
  delete[] erfsd;
  delete[] erfsk;
  delete[] erfku;
  delete[] erfko;
  delete[] erfcv2;
  delete[] erfsd2;
  delete[] erfsk2;
  delete[] erfku2;
  delete[] erfko2;

  delete rg;
  delete min;
  delete[] index;

  return 0;
}


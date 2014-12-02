#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "mini.h"
using namespace std;


int main(int argc, char **argv)
{
  string filename;
  int rep = 100;
  int NIte = 0;
  int seed = 0;
  if (argc > 4)
    {
      rep = atoi(argv[1]);
      NIte = atoi(argv[2]);
      filename.assign(argv[3]);
      seed = atoi(argv[4]);
    }
  else
    {
      cout << "\nusage: compressor [desired replicas] [ite] [PDF set name] [seed]\n" << endl;
      exit(-1);
    }

  // create folder
  mkdir(filename.c_str(), 0777);

  // load PDF set
  const LHAPDF::PDFSet set(filename.c_str());
  vector<LHAPDF::PDF*> pdf = set.mkPDFs();
  vector<int> fPids = pdf[0]->flavors();

  // input vector
  int *index = new int[rep];
  for (int i = 0; i < rep; i++) index[i] = i+1;

  // load a subarray
  double e = 0;
  stringstream log("");
  Mini *min = new Mini(rep, pdf, seed);
  for (int i = 0; i < NIte; i++)
    {
      e = min->iterate(index);
      if ( i % 10 == 0)
        {
          cout << fixed << "ITE: " << i << "\t"
               << scientific << e << endl;
          log << fixed << i << "\t"
               << scientific << e << endl;
        }
    }

  cout << fixed << "Final: " << NIte << "\t"
       << scientific << e << endl;
  log << fixed << NIte << "\t"
       << scientific << e << endl;

  // Compute other estimators
  double ecv = 0, estd = 0, esk = 0, eku = 0, eko = 0;
  double ecv2 = 0, estd2 = 0, esk2 = 0, eku2 = 0, eko2 = 0;
  double *res = new double[6];
  for (int f = 0; f < (int) fPids.size(); f++)
    for (int i = 0; i < (int) min->GetX().size(); i++)
      {
	double cv = 0, std = 0, sk = 0, kur = 0;
	min->ComputeEstimators(rep, min->GetX()[i], Q, fPids[f], index, 
			       cv, std, sk, kur, res);	
	ecv  += pow(min->GetCV(f,i) - cv, 2.0);
	estd += pow(min->GetSD(f,i) - std, 2.0);
	esk  += pow(min->GetSK(f,i) - sk, 2.0);
	eku  += pow(min->GetKU(f,i) - kur, 2.0);
	
	for (int l = 0; l < 6; l++)
	  eko += pow(min->GetKO(f,i,l) - res[l], 2.0);	     

	// for relative difference
	if (fPids[f] == 21 || abs(fPids[f]) == 1 || abs(fPids[f]) == 2 ||
	    abs(fPids[f]) == 3 && (min->GetX()[i] < 0.6 && min->GetX()[i] > 1e-3))
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
  delete[] res;

  cout << scientific;
  cout << "CV:  " << ecv << endl;
  cout << "CV2: " << ecv2 << endl;
  cout << "STD: " << estd << endl;
  cout << "CV2: " << estd2 << endl;
  cout << "SKE: " << esk << endl;
  cout << "CV2: " << esk2 << endl;
  cout << "KUR: " << eku << endl;
  cout << "KU2: " << eku2 << endl;
  cout << "KOL: " << eko << endl;
  cout << "KO2: " << eko2 << endl;

  fstream ss;
  stringstream ff("");
  ff << filename << "/output_erf.dat";
  ss.open(ff.str().c_str(), ios::out|ios::app);
  ss << scientific << ecv << "\t" << estd << "\t" << esk << "\t" << eku << "\t" << eko << endl;
  ss.close();

  stringstream ff2("");
  ff2 << filename << "/output_rel.dat";
  ss.open(ff2.str().c_str(), ios::out|ios::app);
  ss << scientific << ecv2 << "\t" << estd2 << "\t" << esk2 << "\t" << eku2 << "\t" << eko2 << endl;
  ss.close();

  // save erf log
  fstream f;
  stringstream a("");
  a << filename << "/erf.dat";
  f.open(a.str().c_str(), ios::out);
  f << log.str() << endl;
  f.close();

  // save replica id
  stringstream b("");
  b << filename << "/replica.dat";
  f.open(b.str().c_str(), ios::out);
  f << filename << "\t" << rep << endl;
  for (int i = 0; i < rep; i++)
    f << index[i] << endl;
  f.close();

  // cleaning memory
  delete min;
  for (int i = 0; i < (int) pdf.size(); i++)
    if (pdf[i]) delete pdf[i];
  pdf.clear();

  return 0;
}


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
  if (argc > 2)
    {
      rep = atoi(argv[1]);
      filename.assign(argv[2]);
    }
  else
    {
      cout << "\nusage: compressor [desired replicas] [PDF set name]\n" << endl;
      exit(-1);
    }

  // create folder
  mkdir(filename.c_str(), 0777);

  // load PDF set
  const LHAPDF::PDFSet set(filename.c_str());
  vector<LHAPDF::PDF*> pdf = set.mkPDFs();

  // general Setup
  const int NIte = 30000;

  // input vector
  int *index = new int[rep];
  for (int i = 0; i < rep; i++) index[i] = i+1;

  // load a subarray
  double e = 0;
  stringstream log("");
  Mini *min = new Mini(rep, pdf);
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

  min->Save(index, filename);

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


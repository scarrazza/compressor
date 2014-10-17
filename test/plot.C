#include <vector>

void plot()
{
  
  
  vector<string> filename;
  filename.push_back("cv_set_r.dat");
  filename.push_back("sd_set_r.dat");
  filename.push_back("ske_set_r.dat");
  filename.push_back("kur_set_r.dat");  
  filename.push_back("kol_set_r.dat");

  TMultiGraph** mg = new TMultiGraph*[filename.size()];

  TLegend *l = new TLegend(0.5,0.67,0.88,0.88);
  l->SetBorderSize(0);
  l->SetFillColor(0);
  l->SetFillStyle(0);

  for (int i = 0; i < (int) filename.size(); i++)
    {
      stringstream file("");
      file << "set_r/" << filename[i];
      stringstream file2("");
      file2 << "set_c/" << filename[i];

      TGraphErrors *g1 = new TGraphErrors(file.str().c_str(),"%lg %lg %lg");
      TGraphErrors *g2 = new TGraphErrors(file2.str().c_str(),"%lg %lg %lg");
      
      g1->SetMarkerStyle(20);
      g1->SetMarkerColor(kBlue);
      g1->SetLineColor(kBlue);
      
      g2->SetMarkerStyle(21);
      g2->SetMarkerColor(kRed);
      g2->SetLineColor(kRed);

      mg[i] = new TMultiGraph();
      mg[i]->Add(g1,"p");
      mg[i]->Add(g2,"p");
      if (i == 0)
	{
	  l->AddEntry(g1,"Random","pl");
	  l->AddEntry(g2,"Compressed","pl");
	}
    }
 
  TCanvas *c = new TCanvas("c","",1200,700);

  c->Divide(3,2);

  vector<string> label;
  label.push_back("Error Function: Central Value");
  label.push_back("Error Function: Standard Deviation");
  label.push_back("Error Function: Skewness");
  label.push_back("Error Function: Kurtosis");
  label.push_back("Error Function: Kolmogorov");

  double cv[] = {4.410692e+00, 1.399956e+01, 8.294280e+02, 2.027068e+05, 7.014082e+00};
  for (int i = 0; i < (int) filename.size(); i++)
    {
      c->cd(i+1);
      c->cd(i+1)->SetLogx();     
      c->cd(i+1)->SetTickx();
      c->cd(i+1)->SetTicky();
      
      mg[i]->SetTitle(label[i].c_str());
      mg[i]->Draw("A");
      mg[i]->GetXaxis()->SetTitle("Number of Replicas");
      mg[i]->GetXaxis()->CenterTitle(true);      
      l->Draw("same");

      TGraph *g = new TGraph(2);
      g->SetPoint(0, 10, cv[i]);
      g->SetPoint(1, 1000, cv[i]);
      g->SetLineColor(kBlue);
      g->SetLineStyle(2);
      g->Draw("same,L");
    }

}

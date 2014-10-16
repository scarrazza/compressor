void plot()
{
  TGraphErrors *g1 = new TGraphErrors("set_r/cv_set_r.dat","%lg %lg %lg");
  TGraphErrors *g2 = new TGraphErrors("set_r/sd_set_r.dat","%lg %lg %lg");

  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(kBlue);

  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(kRed);  

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(g1,"p");
  mg->Add(g2,"p");
 
  TCanvas *c = new TCanvas();
  c->SetLogx();
  mg->Draw("A");
}

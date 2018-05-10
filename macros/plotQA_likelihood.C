#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH3D.h"
#include "TCanvas.h"
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <utility>
#include <functional>
#include <stdexcept>
using namespace std;

std::pair<int,std::string> get_particle(int pid)
{
  if(pid ==  211)   return std::pair<int,std::string>(pid,"piplus");
  if(pid ==  321)   return std::pair<int,std::string>(pid,"Kplus");
  if(pid ==  2212)  return std::pair<int,std::string>(pid,"proton");
  if(pid == -211)   return std::pair<int,std::string>(pid,"piminus");
  if(pid == -321)   return std::pair<int,std::string>(pid,"Kminus");
  if(pid == -2212)  return std::pair<int,std::string>(pid,"antiproton");
  else return std::pair<int,std::string>(-1,"undefined");
}

std::pair<std::string,std::string> get_misIdentifiedParticle(int pid)
{
  if(pid ==  211)   return std::pair<std::string,std::string>("Kplus","proton");
  if(pid ==  321)   return std::pair<std::string,std::string>("piplus","proton");
  if(pid ==  2212)  return std::pair<std::string,std::string>("piplus","Kplus");
  if(pid == -211)   return std::pair<std::string,std::string>("Kminus","antiproton");
  if(pid == -321)   return std::pair<std::string,std::string>("piminus","antiproton");
  if(pid == -2212)  return std::pair<std::string,std::string>("piminus","Kminus");
  else return std::pair<std::string,std::string>("undefined","undefined");
}

void plotQA_likelihood(int pid = 211)
{
  auto particle = get_particle(pid);
  cout << "particle.first = " << particle.first << ", particle.second = " << particle.second << endl;

  auto misIdentified = get_misIdentifiedParticle(pid);
  cout << "misIdentified.first = " << misIdentified.first << ", misIdentified.second = " << misIdentified.second << endl;

  if(particle.first == -1) return;

  string date = "test/";
  string input_dir = Form("/work/eic/xusun/output/modular_rich/%s",date.c_str());
  string inputfile = Form("%sLogLikelihood.root",input_dir.c_str());
  cout << "read in file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  assert(File_InPut);

  string histQA_misIdentified_1st = Form("h_likelihood_%s%s",particle.second.c_str(),misIdentified.first.c_str());
  TH2D *h_misIdentified_1st = (TH2D*)File_InPut->Get(histQA_misIdentified_1st.c_str());
  assert(h_misIdentified_1st);

  string histQA_misIdentified_2nd = Form("h_likelihood_%s%s",particle.second.c_str(),misIdentified.second.c_str());
  TH2D *h_misIdentified_2nd = (TH2D*)File_InPut->Get(histQA_misIdentified_2nd.c_str());
  assert(h_misIdentified_2nd);

  
  TCanvas *c_likelihood = new TCanvas("c_likelihood","c_likelihood",10,10,1600,800);
  c_likelihood->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_likelihood->cd(i_pad+1);
    c_likelihood->cd(i_pad+1)->SetLeftMargin(0.15);
    c_likelihood->cd(i_pad+1)->SetBottomMargin(0.15);
    c_likelihood->cd(i_pad+1)->SetGrid(0,0);
    c_likelihood->cd(i_pad+1)->SetTicks(1,1);
  }

  c_likelihood->cd(1);
  string label_1st = Form("Log(%s) - Log(%s)",particle.second.c_str(),misIdentified.first.c_str());
  h_misIdentified_1st->SetTitle(particle.second.c_str());
  h_misIdentified_1st->SetStats(0);

  h_misIdentified_1st->GetXaxis()->SetTitle("generated momentum (GeV/c)");
  h_misIdentified_1st->GetXaxis()->CenterTitle();
  h_misIdentified_1st->GetXaxis()->SetTitleOffset(1.25);
  h_misIdentified_1st->GetXaxis()->SetNdivisions(505);

  h_misIdentified_1st->GetYaxis()->SetTitle(label_1st.c_str());
  h_misIdentified_1st->GetYaxis()->CenterTitle();
  h_misIdentified_1st->GetYaxis()->SetTitleOffset(1.5);
  h_misIdentified_1st->GetYaxis()->SetNdivisions(505);
  h_misIdentified_1st->Draw("colz");

  c_likelihood->cd(2);
  string label_2nd = Form("Log(%s) - Log(%s)",particle.second.c_str(),misIdentified.second.c_str());
  h_misIdentified_2nd->SetTitle(particle.second.c_str());
  h_misIdentified_2nd->SetStats(0);

  h_misIdentified_2nd->GetXaxis()->SetTitle("generated momentum (GeV/c)");
  h_misIdentified_2nd->GetXaxis()->CenterTitle();
  h_misIdentified_2nd->GetXaxis()->SetTitleOffset(1.25);
  h_misIdentified_2nd->GetXaxis()->SetNdivisions(505);

  h_misIdentified_2nd->GetYaxis()->SetTitle(label_2nd.c_str());
  h_misIdentified_2nd->GetYaxis()->CenterTitle();
  h_misIdentified_2nd->GetYaxis()->SetTitleOffset(1.5);
  h_misIdentified_2nd->GetYaxis()->SetNdivisions(505);
  h_misIdentified_2nd->Draw("colz");

  string outputfile = Form("../figures/likelihood_diff/loglikelihood_difference_%s.pdf",particle.second.c_str());
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_likelihood->SaveAs(outputfile.c_str());
}

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
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

void plotQA_prob(int pid = 211)
{
  auto particle = get_particle(pid);
  cout << "particle.first = " << particle.first << ", particle.second = " << particle.second << endl;

  auto misIdentified = get_misIdentifiedParticle(pid);
  cout << "misIdentified.first = " << misIdentified.first << ", misIdentified.second = " << misIdentified.second << endl;

  if(particle.first == -1) return;

  string date = "May10_2018";
  string inputfile = Form("/work/eic/xusun/output/likelihood/LogLikelihood_%s.root",date.c_str());
  cout << "read in file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  assert(File_InPut);

  string gQA_identified = Form("g_prob_%s%s",particle.second.c_str(),particle.second.c_str());
  TGraph *g_prob_identified = (TGraph*)File_InPut->Get(gQA_identified.c_str());
  assert(g_prob_identified);

  string gQA_misIdentified_1st = Form("g_prob_%s%s",particle.second.c_str(),misIdentified.first.c_str());
  TGraph *g_prob_misIdentified_1st = (TGraph*)File_InPut->Get(gQA_misIdentified_1st.c_str());
  assert(g_prob_misIdentified_1st);

  string gQA_misIdentified_2nd = Form("g_prob_%s%s",particle.second.c_str(),misIdentified.second.c_str());
  TGraph *g_prob_misIdentified_2nd = (TGraph*)File_InPut->Get(gQA_misIdentified_2nd.c_str());
  assert(g_prob_misIdentified_2nd);

  TCanvas *c_probability = new TCanvas("c_probability","c_probability",10,10,800,800);
  c_probability->cd()->SetLeftMargin(0.15);
  c_probability->cd()->SetBottomMargin(0.15);
  c_probability->cd()->SetGrid(0,0);
  c_probability->cd()->SetTicks(1,1);

  TH1D *h_play = new TH1D("h_play","h_play",100,-0.5,19.5);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("generated momemtum (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleOffset(1.25);
  h_play->GetXaxis()->SetNdivisions(505);

  h_play->GetYaxis()->SetTitle("particle identified probability");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitleOffset(1.25);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(-0.05,1.1);
  h_play->Draw("pE");

  g_prob_identified->SetMarkerStyle(24);
  g_prob_identified->SetMarkerColor(kGray+3);
  g_prob_identified->SetMarkerSize(1.4);
  g_prob_identified->Draw("P Same");
  string leg_identified = Form("%s to %s",particle.second.c_str(),particle.second.c_str());

  g_prob_misIdentified_1st->SetMarkerStyle(24);
  g_prob_misIdentified_1st->SetMarkerColor(2);
  g_prob_misIdentified_1st->SetMarkerSize(1.5);
  g_prob_misIdentified_1st->Draw("P Same");
  string leg_misIdentified_1st = Form("%s to %s",particle.second.c_str(),misIdentified.first.c_str());

  g_prob_misIdentified_2nd->SetMarkerStyle(20);
  g_prob_misIdentified_2nd->SetMarkerColor(4);
  g_prob_misIdentified_2nd->SetMarkerSize(1.0);
  g_prob_misIdentified_2nd->Draw("P Same");
  string leg_misIdentified_2nd = Form("%s to %s",particle.second.c_str(),misIdentified.second.c_str());

  TLegend *leg = new TLegend(0.4,0.4,0.8,0.7);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry(g_prob_identified,leg_identified.c_str(),"P");
  leg->AddEntry(g_prob_misIdentified_1st,leg_misIdentified_1st.c_str(),"P");
  leg->AddEntry(g_prob_misIdentified_2nd,leg_misIdentified_2nd.c_str(),"P");
  leg->Draw("same");

  TLine *line = new TLine(-0.5,0.0,19.5,0.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");

  string outputfile = Form("../figures/probability_diff/PID_probability_%s.pdf",particle.second.c_str());
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_probability->SaveAs(outputfile.c_str());
}

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2D.h"
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

void plotQA_PDF()
{
  string date = "May23_2018";
  // string inputfile = Form("/work/eic/xusun/output/database/3mm/database_%s.root.eff.wo.smearing",date.c_str());
  string inputfile = Form("/work/eic/xusun/output/database/database_%s.root",date.c_str());
  cout << "read in file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  assert(File_InPut);

  string key_pion = "h_mPhotonDist_piplus_vx_0_vy_0_mom_2_theta_2_phi_3";
  TH2D *h_PDF_pion = (TH2D*)File_InPut->Get(key_pion.c_str());
  assert(h_PDF_pion);

  string key_kaon = "h_mPhotonDist_Kplus_vx_0_vy_0_mom_2_theta_2_phi_3";
  TH2D *h_PDF_kaon = (TH2D*)File_InPut->Get(key_kaon.c_str());
  assert(h_PDF_kaon);

  string key_proton = "h_mPhotonDist_proton_vx_0_vy_0_mom_2_theta_2_phi_3";
  TH2D *h_PDF_proton = (TH2D*)File_InPut->Get(key_proton.c_str());
  assert(h_PDF_proton);

  TCanvas *c_PDF = new TCanvas("c_PDF","c_PDF",10,10,1500,500);
  c_PDF->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_PDF->cd(i_pad+1);
    c_PDF->cd(i_pad+1)->SetLeftMargin(0.15);
    c_PDF->cd(i_pad+1)->SetBottomMargin(0.15);
    c_PDF->cd(i_pad+1)->SetGrid(0,0);
    c_PDF->cd(i_pad+1)->SetTicks(1,1);
  }

  c_PDF->cd(1);
  h_PDF_pion->SetTitle("pi+ mass hypotheses");
  h_PDF_pion->SetStats(0);

  h_PDF_pion->GetXaxis()->SetTitle("out_x (mm)");
  h_PDF_pion->GetXaxis()->CenterTitle();
  h_PDF_pion->GetXaxis()->SetTitleOffset(1.25);
  h_PDF_pion->GetXaxis()->SetNdivisions(505);

  h_PDF_pion->GetYaxis()->SetTitle("out_y (mm)");
  h_PDF_pion->GetYaxis()->CenterTitle();
  h_PDF_pion->GetYaxis()->SetTitleOffset(1.25);
  h_PDF_pion->GetYaxis()->SetNdivisions(505);

  h_PDF_pion->Draw("colz");

  c_PDF->cd(2);
  h_PDF_kaon->SetTitle("K+ mass hypotheses");
  h_PDF_kaon->SetStats(0);

  h_PDF_kaon->GetXaxis()->SetTitle("out_x (mm)");
  h_PDF_kaon->GetXaxis()->CenterTitle();
  h_PDF_kaon->GetXaxis()->SetTitleOffset(1.25);
  h_PDF_kaon->GetXaxis()->SetNdivisions(505);

  h_PDF_kaon->GetYaxis()->SetTitle("out_y (mm)");
  h_PDF_kaon->GetYaxis()->CenterTitle();
  h_PDF_kaon->GetYaxis()->SetTitleOffset(1.25);
  h_PDF_kaon->GetYaxis()->SetNdivisions(505);

  h_PDF_kaon->Draw("colz");

  c_PDF->cd(3);
  h_PDF_proton->SetTitle("p mass hypotheses");
  h_PDF_proton->SetStats(0);

  h_PDF_proton->GetXaxis()->SetTitle("out_x (mm)");
  h_PDF_proton->GetXaxis()->CenterTitle();
  h_PDF_proton->GetXaxis()->SetTitleOffset(1.25);
  h_PDF_proton->GetXaxis()->SetNdivisions(505);

  h_PDF_proton->GetYaxis()->SetTitle("out_y (mm)");
  h_PDF_proton->GetYaxis()->CenterTitle();
  h_PDF_proton->GetYaxis()->SetTitleOffset(1.25);
  h_PDF_proton->GetYaxis()->SetNdivisions(505);

  h_PDF_proton->Draw("colz");

  string outputfile = "../figures/masshypo_dist/mass_hypotheses.pdf";
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_PDF->SaveAs(outputfile.c_str());

  TCanvas *c_Identified = new TCanvas("c_Identified","c_Identified",10,10,800,800);
  c_Identified->cd();
  c_Identified->cd()->SetLeftMargin(0.15);
  c_Identified->cd()->SetBottomMargin(0.15);
  c_Identified->cd()->SetGrid(0,0);
  c_Identified->cd()->SetTicks(1,1);

  h_PDF_kaon->Draw("colz");
  // h_PDF_proton->Draw("same");
  h_PDF_pion->Draw("same ");

  outputfile = "../figures/masshypo_dist/mass_hypotheses_identify.pdf";
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_Identified->SaveAs(outputfile.c_str());
}

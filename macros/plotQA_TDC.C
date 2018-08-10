#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
using namespace std;

void plotQA_TDC(const int runID = 182, const string mode = "rich")
{
  string inputfile = Form("/Users/xusun/Data/BeamTestData/suite1.0/results/tdc/%sTDC_run%d/sspRich.root",mode.c_str(),runID);
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  assert(File_InPut);

  TCanvas *c_TDC = new TCanvas("c_TDC","c_TDC",10,10,800,800);
  c_TDC->cd();
  c_TDC->cd()->SetLeftMargin(0.15);
  c_TDC->cd()->SetBottomMargin(0.15);
  c_TDC->cd()->SetGrid(0,0);
  c_TDC->cd()->SetTicks(1,1);

  TH1F *h_TDC = new TH1F("h_TDC","h_TDC",5000,-0.5,4999.5);
  TTree * tree_mRICH = (TTree*)File_InPut->Get("data");
  tree_mRICH->Draw("time>>h_TDC","pol==1");

  // TLegend *leg = new TLegend(0.2,0.5,0.5,0.65);
  // leg->SetFillColor(0);
  // leg->SetBorderSize(0);
  // leg->AddEntry(h_PMT,"PMT time","L");
  // leg->AddEntry(h_SiPM,"SiPM time + 1490","L");
  // leg->Draw("same");
}

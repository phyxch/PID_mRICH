#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
using namespace std;

void plotQA_TDC(const int runID = 356, const string mode = "sipm")
{
  // string inputfile = Form("/Users/xusun/Data/BeamTestData/suite1.0/results/tdc/%sTDC_run%d/sspRich.root",mode.c_str(),runID);
  string inputfile = Form("/home/xusun/Data/mRICH/BeamTest/tdc/%sTDC_run%d/sspRich.root",mode.c_str(),runID);
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

  string c_tdc = Form("../../figures/BeamTest_QA/c_TDC_MPPC_%d_allChanel.eps",runID);
  c_TDC->SaveAs(c_tdc.c_str());
}

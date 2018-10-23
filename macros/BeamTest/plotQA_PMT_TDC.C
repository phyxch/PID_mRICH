#include "string"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

void plotQA_PMT_TDC(const int runID = 182)
{
  float tdc_Start = 2000.0;
  float tdc_Stop  = 2050.0;

  float ratio_cut = 2.0; // run dependent
  if(runID == 24) ratio_cut = 1.8; // run dependent
  if(runID == 120) ratio_cut = 1.8; // run dependent
  if(runID == 124) ratio_cut = 1.5; // run dependent
  if(runID == 214) ratio_cut = 2.0; // run dependent
  if(runID == 270) ratio_cut = 3.0;
  if(runID == 314) ratio_cut = 1.5;
  if(runID > 343 && runID < 381) // meson run 344-380
  {
    ratio_cut = 3.0;
    tdc_Start = 490.0;
    tdc_Stop  = 590.0;
  }

  int const NumOfPixel = 33;
  string inputfile = Form("/home/xusun/Data/mRICH/BeamTest/QA/PMT/richTDC_run%d.root",runID);
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH2F *h_mRingImage = (TH2F*)File_InPut->Get("h_mRingImage")->Clone();
  TH1F *h_mTDC[NumOfPixel][NumOfPixel]; // 0 for x-pixel | 1 for y-pixel
  for(int i_pixel_x = 0; i_pixel_x < NumOfPixel; ++i_pixel_x)
  {
    for(int i_pixel_y = 0; i_pixel_y < NumOfPixel; ++i_pixel_y)
    {
      string HistName = Form("h_mTDC_pixelX_%d_pixelY_%d",i_pixel_x,i_pixel_y);
      h_mTDC[i_pixel_x][i_pixel_y] = (TH1F*)File_InPut->Get(HistName.c_str())->Clone();
    }
  }

  TCanvas *c_RingImage = new TCanvas("c_RingImage","c_RingImage",10,10,NumOfPixel*30,NumOfPixel*30);
  c_RingImage->SetLeftMargin(0.15);
  c_RingImage->SetBottomMargin(0.15);
  c_RingImage->SetRightMargin(0.15);
  c_RingImage->SetTicks(1,1);
  c_RingImage->SetGrid(0,0);
  string title = Form("120 GeV proton & run%d",runID);
  if(runID > 343 && runID < 381) title = Form("meson & run%d",runID);
  h_mRingImage->SetTitle(title.c_str());
  h_mRingImage->SetStats(0);
  h_mRingImage->GetXaxis()->SetTitle("pixel ID");
  h_mRingImage->GetXaxis()->CenterTitle();
  h_mRingImage->GetYaxis()->SetTitle("pixel ID");
  h_mRingImage->GetYaxis()->CenterTitle();
  h_mRingImage->Draw("colz");
  string c_ringimage = Form("../../figures/BeamTest_QA/c_RingImage_PMT_%d.eps",runID);
  c_RingImage->SaveAs(c_ringimage.c_str());
  c_ringimage = Form("../../figures/BeamTest_QA/c_RingImage_PMT_%d.png",runID);
  c_RingImage->SaveAs(c_ringimage.c_str());

  // TCanvas *c_TDC = new TCanvas("c_TDC","c_TDC",10,10,NumOfPixel*100,NumOfPixel*100);
  TCanvas *c_TDC = new TCanvas("c_TDC","c_TDC",10,10,NumOfPixel*30,NumOfPixel*30);
  c_TDC->Divide(NumOfPixel,NumOfPixel,0,0);
  for(int i_pixel_x = 0; i_pixel_x < NumOfPixel; ++i_pixel_x)
  {
    for(int i_pixel_y = 0; i_pixel_y < NumOfPixel; ++i_pixel_y)
    {
      double max_counts = h_mTDC[i_pixel_x][i_pixel_y]->GetMaximum();
      int i_pad = NumOfPixel*(NumOfPixel-(i_pixel_y+1)) + i_pixel_x+1;
      // cout << "i_pixel_x = " << i_pixel_x << ", i_pixel_y = " << i_pixel_y << ", i_pad = " << i_pad << endl;
      c_TDC->cd(i_pad);
      h_mTDC[i_pixel_x][i_pixel_y]->SetTitle("");
      h_mTDC[i_pixel_x][i_pixel_y]->SetStats(0);
      float tot = h_mTDC[i_pixel_x][i_pixel_y]->Integral();
      float sig = h_mTDC[i_pixel_x][i_pixel_y]->Integral(tdc_Start,tdc_Stop);
      float bkg = tot - sig;
      float ratio = sig/bkg;
      h_mTDC[i_pixel_x][i_pixel_y]->GetXaxis()->SetRangeUser(tdc_Start-100,tdc_Stop+100);
      h_mTDC[i_pixel_x][i_pixel_y]->GetYaxis()->SetRangeUser(0.1,1e2);

      if(ratio > ratio_cut) h_mTDC[i_pixel_x][i_pixel_y]->SetLineColor(2);
      else h_mTDC[i_pixel_x][i_pixel_y]->SetLineColor(1);

      h_mTDC[i_pixel_x][i_pixel_y]->Draw();
    }
  }
  string c_tdc = Form("../../figures/BeamTest_QA/c_TDC_PMT_%d.eps",runID);
  c_TDC->SaveAs(c_tdc.c_str());
  c_tdc = Form("../../figures/BeamTest_QA/c_TDC_PMT_%d.png",runID);
  c_TDC->SaveAs(c_tdc.c_str());
}

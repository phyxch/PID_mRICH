#include "string"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "../draw.h"

void plotQA_PMT_TimeCuts(const int runID = 182)
{
  float tdc_Start = 2000.0;
  float tdc_Stop  = 2050.0;

  float ratio_cut = 2.0; // run dependent
  int x_OnRing = 10;
  int y_OnRing = 5;
  int x_OffRing = 20;
  int y_OffRing = 10;
  int x_beam = 15;
  int y_beam = 15;

  if(runID == 24)
  {
    ratio_cut = 1.8; // run dependent
    x_OnRing = 10;
    y_OnRing = 3;
    x_OffRing = 20;
    y_OffRing = 10;
    x_beam = 8;
    y_beam = 23;
  }

  if(runID == 120)
  {
    ratio_cut = 1.8; // run dependent
    x_OnRing = 7;
    y_OnRing = 6;
    x_OffRing = 30;
    y_OffRing = 30;
    x_beam = 17;
    y_beam = 14;
  }

  if(runID == 124)
  {
    ratio_cut = 1.5; // run dependent
    x_OnRing = 9;
    y_OnRing = 5;
    x_OffRing = 30;
    y_OffRing = 30;
    x_beam = 24;
    y_beam = 8;
  }

  if(runID == 214)
  {
    ratio_cut = 2.0; // run dependent
    x_OnRing = 4;
    y_OnRing = 11;
    x_OffRing = 30;
    y_OffRing = 30;
    x_beam = 17;
    y_beam = 3;
  }

  if(runID > 343 && runID < 381) // meson run 344-380
  {
    ratio_cut = 3.0;
    tdc_Start = 490.0;
    tdc_Stop  = 590.0;
  }

  int const NumOfPixel = 33;
  string inputfile = Form("/home/xusun/Data/mRICH/BeamTest/QA/richTDC_run%d.root",runID);
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH2F *h_ratio_all = new TH2F("h_ratio_all","h_ratio_all",NumOfPixel,-0.5,32.5,NumOfPixel,-0.5,32.5);
  TH2F *h_ratio_cut = new TH2F("h_ratio_cut","h_ratio_cut",NumOfPixel,-0.5,32.5,NumOfPixel,-0.5,32.5);

  TH2F *h_mean_all = new TH2F("h_mean_all","h_mean_all",NumOfPixel,-0.5,32.5,NumOfPixel,-0.5,32.5);
  TH2F *h_mean_cut = new TH2F("h_mean_cut","h_mean_cut",NumOfPixel,-0.5,32.5,NumOfPixel,-0.5,32.5);

  TH2F *h_sigma_all = new TH2F("h_sigma_all","h_sigma_all",NumOfPixel,-0.5,32.5,NumOfPixel,-0.5,32.5);
  TH2F *h_sigma_cut = new TH2F("h_sigma_cut","h_sigma_cut",NumOfPixel,-0.5,32.5,NumOfPixel,-0.5,32.5);

  TH1F *h_mTDC[NumOfPixel][NumOfPixel]; // 0 for x-pixel | 1 for y-pixel
  float mean_tdc = 0.0;
  float sigma_tdc = 0.0;
  int counter_tdc = 0;
  for(int i_pixel_x = 0; i_pixel_x < NumOfPixel; ++i_pixel_x)
  {
    for(int i_pixel_y = 0; i_pixel_y < NumOfPixel; ++i_pixel_y)
    {
      // if( !(i_pixel_x > 12 && i_pixel_x < 20 && i_pixel_y > 12 && i_pixel_y < 20) )  continue;
      string HistName = Form("h_mTDC_pixelX_%d_pixelY_%d",i_pixel_x,i_pixel_y);
      h_mTDC[i_pixel_x][i_pixel_y] = (TH1F*)File_InPut->Get(HistName.c_str())->Clone();
      float counts_tot = h_mTDC[i_pixel_x][i_pixel_y]->Integral();

      if(counts_tot < 0.5) continue;

      TF1 *f_gaus = new TF1("f_gaus","gaus",0,5000);
      f_gaus->SetParameter(0,100);
      f_gaus->SetParameter(1,0.5*(tdc_Start+tdc_Stop));
      f_gaus->SetParameter(2,0.5*(tdc_Stop-tdc_Start));
      f_gaus->SetRange(tdc_Start,tdc_Stop);
      h_mTDC[i_pixel_x][i_pixel_y]->Fit(f_gaus,"MNQR");

      float norm = f_gaus->GetParameter(0);
      float mean = f_gaus->GetParameter(1);
      float sigma = f_gaus->GetParameter(2);

      int sig_Start = h_mTDC[i_pixel_x][i_pixel_y]->FindBin(mean-3.0*sigma);
      int sig_Stop  = h_mTDC[i_pixel_x][i_pixel_y]->FindBin(mean+3.0*sigma);
      float counts_sig = h_mTDC[i_pixel_x][i_pixel_y]->Integral(sig_Start,sig_Stop);
      float counts_bkg = counts_tot - counts_sig;
      float ratio = counts_sig/counts_bkg;

      h_ratio_all->Fill(i_pixel_x,i_pixel_y,ratio);
      h_mean_all->Fill(i_pixel_x,i_pixel_y,mean);
      h_sigma_all->Fill(i_pixel_x,i_pixel_y,sigma);

      if(ratio > ratio_cut)
      {
	h_ratio_cut->Fill(i_pixel_x,i_pixel_y,ratio);
	h_mean_cut->Fill(i_pixel_x,i_pixel_y,mean);
	h_sigma_cut->Fill(i_pixel_x,i_pixel_y,sigma);
	// cout << "sig_Start = " << sig_Start << ", sig_Stop = " << sig_Stop << ", counts_sig = " << counts_sig << ", counts_tot = " << counts_tot << ", ratio = " << ratio << endl;

	// calculate tdc cuts with sig/bkg > 2.0
	mean_tdc += mean;
	sigma_tdc += sigma;
	counter_tdc++;
      }
    }
  }

  TCanvas *c_TimeCut = new TCanvas("c_TimeCut","c_TimeCut",10,10,900,900);
  c_TimeCut->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    // c_TimeCut->cd(i_pad+1)->SetTopMargin(0.15);
    // c_TimeCut->cd(i_pad+1)->SetBottomMargin(0.15);
    // c_TimeCut->cd(i_pad+1)->SetLeftMargin(0.15);
    c_TimeCut->cd(i_pad+1)->SetRightMargin(0.15);
    c_TimeCut->cd(i_pad+1)->SetTicks(1,1);
    c_TimeCut->cd(i_pad+1)->SetGrid(0,0);
  }

  float yields_beam = h_mTDC[x_beam][y_beam]->GetMaximum();
  float yields_ring = h_mTDC[x_OnRing][y_OnRing]->GetMaximum();
  float yields = (yields_beam > yields_ring) ? yields_beam : yields_ring;

  float mean_tdc_Start = mean_tdc/counter_tdc - 3.0*sigma_tdc/counter_tdc;
  float mean_tdc_Stop  = mean_tdc/counter_tdc + 3.0*sigma_tdc/counter_tdc;

  c_TimeCut->cd(1); // on ring
  h_mTDC[x_OnRing][y_OnRing]->SetTitle("On Ring");
  h_mTDC[x_OnRing][y_OnRing]->GetXaxis()->SetRangeUser(tdc_Start-100,tdc_Stop+100);
  h_mTDC[x_OnRing][y_OnRing]->GetYaxis()->SetRangeUser(0,1.1*yields);
  h_mTDC[x_OnRing][y_OnRing]->Draw();

  TF1 *f_gaus_OnRing = new TF1("f_gaus_OnRing","gaus",0,5000);
  f_gaus_OnRing->SetParameter(0,100);
  f_gaus_OnRing->SetParameter(1,0.5*(tdc_Start+tdc_Stop));
  f_gaus_OnRing->SetParameter(2,0.5*(tdc_Stop-tdc_Start));
  f_gaus_OnRing->SetRange(tdc_Start,tdc_Stop);
  h_mTDC[x_OnRing][y_OnRing]->Fit(f_gaus_OnRing,"MQR");
  PlotLine(mean_tdc_Start,mean_tdc_Start,0,0.5*yields,4,2,2);
  PlotLine(mean_tdc_Stop,mean_tdc_Stop,0,0.5*yields,4,2,2);

  c_TimeCut->cd(2); // off ring
  h_mTDC[x_OffRing][y_OffRing]->SetTitle("Off Ring");
  h_mTDC[x_OffRing][y_OffRing]->GetXaxis()->SetRangeUser(tdc_Start-100,tdc_Stop+100);
  h_mTDC[x_OffRing][y_OffRing]->GetYaxis()->SetRangeUser(0,1.1*yields);
  h_mTDC[x_OffRing][y_OffRing]->Draw();
  PlotLine(mean_tdc_Start,mean_tdc_Start,0,0.5*yields,4,2,2);
  PlotLine(mean_tdc_Stop,mean_tdc_Stop,0,0.5*yields,4,2,2);

  c_TimeCut->cd(3); // off ring
  h_mTDC[x_beam][y_beam]->SetTitle("Beam Spot");
  h_mTDC[x_beam][y_beam]->GetXaxis()->SetRangeUser(tdc_Start-100,tdc_Stop+100);
  h_mTDC[x_beam][y_beam]->GetYaxis()->SetRangeUser(0,1.1*yields);
  h_mTDC[x_beam][y_beam]->Draw();
  PlotLine(mean_tdc_Start,mean_tdc_Start,0,0.5*yields,4,2,2);
  PlotLine(mean_tdc_Stop,mean_tdc_Stop,0,0.5*yields,4,2,2);

  c_TimeCut->cd(4);
  h_ratio_all->SetTitle("sig/bkg");
  h_ratio_all->SetStats(0);
  h_ratio_all->GetXaxis()->SetTitle("pixel ID");
  h_ratio_all->GetXaxis()->CenterTitle();
  h_ratio_all->GetYaxis()->SetTitle("pixel ID");
  h_ratio_all->GetYaxis()->CenterTitle();
  h_ratio_all->Draw("colz");

  c_TimeCut->cd(5);
  h_mean_all->SetTitle("mean");
  h_mean_all->SetStats(0);
  h_mean_all->GetXaxis()->SetTitle("pixel ID");
  h_mean_all->GetXaxis()->CenterTitle();
  h_mean_all->GetYaxis()->SetTitle("pixel ID");
  h_mean_all->GetYaxis()->CenterTitle();
  h_mean_all->Draw("colz");

  c_TimeCut->cd(6);
  h_sigma_all->SetTitle("sigma");
  h_sigma_all->SetStats(0);
  h_sigma_all->GetXaxis()->SetTitle("pixel ID");
  h_sigma_all->GetXaxis()->CenterTitle();
  h_sigma_all->GetYaxis()->SetTitle("pixel ID");
  h_sigma_all->GetYaxis()->CenterTitle();
  h_sigma_all->Draw("colz");

  string str_sig = Form("sig/bkg > %1.1f",ratio_cut);
  c_TimeCut->cd(7);
  h_ratio_cut->SetTitle(str_sig.c_str());
  h_ratio_cut->SetStats(0);
  h_ratio_cut->GetXaxis()->SetTitle("pixel ID");
  h_ratio_cut->GetXaxis()->CenterTitle();
  h_ratio_cut->GetYaxis()->SetTitle("pixel ID");
  h_ratio_cut->GetYaxis()->CenterTitle();
  h_ratio_cut->Draw("colz");

  c_TimeCut->cd(8);
  h_mean_cut->SetTitle("mean");
  h_mean_cut->SetStats(0);
  h_mean_cut->GetXaxis()->SetTitle("pixel ID");
  h_mean_cut->GetXaxis()->CenterTitle();
  h_mean_cut->GetYaxis()->SetTitle("pixel ID");
  h_mean_cut->GetYaxis()->CenterTitle();
  h_mean_cut->Draw("colz");
  string str_mean = Form("mean tdc = %.1f",mean_tdc/counter_tdc);
  plotTopLegend((char*)str_sig.c_str(),0.35,0.55,0.06,1,0.0,42,1,1);
  plotTopLegend((char*)str_mean.c_str(),0.3,0.45,0.06,1,0.0,42,1,1);

  c_TimeCut->cd(9);
  h_sigma_cut->SetTitle("sigma");
  h_sigma_cut->SetStats(0);
  h_sigma_cut->GetXaxis()->SetTitle("pixel ID");
  h_sigma_cut->GetXaxis()->CenterTitle();
  h_sigma_cut->GetYaxis()->SetTitle("pixel ID");
  h_sigma_cut->GetYaxis()->CenterTitle();
  h_sigma_cut->Draw("colz");
  string str_sigma = Form("sigma tdc = %.1f",sigma_tdc/counter_tdc);
  plotTopLegend((char*)str_sig.c_str(),0.35,0.55,0.06,1,0.0,42,1,1);
  plotTopLegend((char*)str_sigma.c_str(),0.3,0.45,0.06,1,0.0,42,1,1);

  string c_timecut = Form("../../figures/BeamTest_QA/c_TimeCuts_PMT_%d.eps",runID);
  c_TimeCut->SaveAs(c_timecut.c_str());

  TH2F *h_mTimeCuts = new TH2F("h_mTimeCuts","h_mTimeCuts",3,-0.5,2.5,800,-0.5,799.5);
  h_mTimeCuts->SetBinContent(1,runID,floor(mean_tdc_Start));
  h_mTimeCuts->SetBinContent(2,runID,ceil(mean_tdc_Stop));
  h_mTimeCuts->SetBinContent(3,runID,runID);
  // cout << "mean_tdc_Start = " << mean_tdc_Start << ", floor = " << floor(mean_tdc_Start) << ", mean_tdc_Stop = " << mean_tdc_Stop << ", ceil = " << ceil(mean_tdc_Stop) << endl;

  string outputfile = Form("/home/xusun/Data/mRICH/BeamTest/QA/richTimeCuts_run%d.root",runID);
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_mTimeCuts->Write();
  File_OutPut->Close();
}

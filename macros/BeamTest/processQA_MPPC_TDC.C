#include "string"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>

#define MAXEDGE 100000

// #define H13700MAP "h13700.txt"
#define H13700MAP "../../include/H13700_180degree_v2.txt"

//MAROC configuration polarity (from ssptest_TDCAll.c)
#define MAROCPOLARITY 1

// global variables for display
int x_mRICH[256];
int y_mRICH[256];
int xp_mRICH[4]; //  PMT coordinates X coordinate of pixel 1 of each mapmt
int yp_mRICH[4]; // Y coordinate of pixel 1 of each mapmt

// translation map MAROC/hamamatsu
int maroc2anode[] = {60,58,59,57,52,50,51,49,44,42,43,41,36,34,35,33,28,26,27,25,20,18,19,17,12,10,11,9,4,2,3,1,5,7,6,8,13,15,14,16,21,23,22,24,29,31,30,32,37,39,38,40,45,47,46,48,53,55,54,56,61,63,62,64};
int maroc2h13700[384];

unsigned int tTrigNum;
double tTrigTime;//long unsigned int tTimestamp;
unsigned int tNedge;

unsigned int tPolarity[MAXEDGE];
unsigned int tChannel[MAXEDGE];
unsigned int tTime[MAXEDGE];
unsigned int tSlot[MAXEDGE];
unsigned int tFiber[MAXEDGE];

void InitDisplay_mRICH();
void ResetEventData();
int GetPMT_mRICH(int slot,int fiber,int asic);
void GenCoord_mRICH(int ipmt, int x1, int y1);
int GetPixel_mRICH(int fiber, int asic, int maroc_channel);

void processQA_MPPC_TDC(const int runID = 672)
{
  int debug = 1;
  int const NumOfPixel = 33;
  // string inputfile = Form("/Users/xusun/Data/BeamTestData/suite1.0/results/tdc/%sTDC_run%d/sspRich.root",mode.c_str(),runID);
  string inputfile = Form("/home/xusun/Data/mRICH/BeamTest/tdc/sipmTDC_run%d/sspRich.root",runID);
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  InitDisplay_mRICH();
  TH2F *h_mRingImage = new TH2F("h_mRingImage","h_mRingImage",NumOfPixel,-0.5,32.5,NumOfPixel,-0.5,32.5);
  TH1F *h_mTDC[NumOfPixel][NumOfPixel]; // 0 for x-pixel | 1 for y-pixel
  for(int i_pixel_x = 0; i_pixel_x < NumOfPixel; ++i_pixel_x)
  {
    for(int i_pixel_y = 0; i_pixel_y < NumOfPixel; ++i_pixel_y)
    {
      string HistName = Form("h_mTDC_pixelX_%d_pixelY_%d",i_pixel_x,i_pixel_y);
      h_mTDC[i_pixel_x][i_pixel_y] = new TH1F(HistName.c_str(),HistName.c_str(),1500,-0.5,1499.5);
    }
  }

  unsigned int pol=MAROCPOLARITY; // 1 falling, 0 rising

  TTree * tree_mRICH = (TTree*)File_InPut->Get("data");
  tree_mRICH->SetBranchAddress("evt",&tTrigNum);
  tree_mRICH->SetBranchAddress("trigtime",&tTrigTime);
  tree_mRICH->SetBranchAddress("nedge",&tNedge);
  tree_mRICH->SetBranchAddress("slot",tSlot);
  tree_mRICH->SetBranchAddress("fiber",tFiber);
  tree_mRICH->SetBranchAddress("ch",tChannel);
  tree_mRICH->SetBranchAddress("pol",tPolarity);
  tree_mRICH->SetBranchAddress("time",tTime);

  int NumOfEvents = tree_mRICH->GetEntries();
  if(NumOfEvents > 50000) NumOfEvents = 50000;
  // int NumOfEvents = 10000;
  printf("NEntries %d\n",NumOfEvents);

  tree_mRICH->GetEntry(0);
  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  {
    if(NumOfEvents>20)if(i_event%(NumOfEvents/10)==0)printf("Processing Event %6d\n",i_event);
    ResetEventData();
    tree_mRICH->GetEntry(i_event);

    if(tNedge>MAXEDGE){
      printf("Event to big: %u edges vs %u array size...skip\n",tNedge,MAXEDGE);
      continue;
    }

    for(unsigned int i_photon = 0; i_photon < tNedge; ++i_photon)
    {
      int slot = tSlot[i_photon];
      int fiber = tFiber[i_photon];
      int channel = tChannel[i_photon];
      int asic = channel/64;
      int pin = channel%64;

      if(tSlot[i_photon] < 3 || tSlot[i_photon] > 7){printf("%s EVT %d Data Error: bad slot %d \n",__FUNCTION__,i_event,slot);continue;}
      if(tFiber[i_photon] > 31){printf("%s EVT %d Data Error: bad fiber %d \n",__FUNCTION__,i_event,fiber);continue;}
      if(tChannel[i_photon] > 191){printf("%s EVT %d Data Error: bad channel %d \n",__FUNCTION__,i_photon,channel); continue;}

      int pmt = GetPMT_mRICH(slot,fiber,asic);
      GenCoord_mRICH(pmt, xp_mRICH[pmt-1], yp_mRICH[pmt-1]);
      int pixel = GetPixel_mRICH(fiber, asic, pin);
      int pixel_x = x_mRICH[pixel-1];
      int pixel_y = y_mRICH[pixel-1];
      if(tPolarity[i_photon] == pol) h_mTDC[pixel_x][pixel_y]->Fill(tTime[i_photon]);

      if(tPolarity[i_photon] == pol && tTime[i_photon] > 500 && tTime[i_photon] < 570) // MPPC
      {
	h_mRingImage->Fill(x_mRICH[pixel-1],y_mRICH[pixel-1]);
      }
    }
  }
  printf("Processed events %d\n",NumOfEvents);

  string outputfile = Form("/home/xusun/Data/mRICH/BeamTest/QA/sipmTDC_run%d.root",runID);
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_mRingImage->Write();
  for(int i_pixel_x = 0; i_pixel_x < NumOfPixel; ++i_pixel_x)
  {
    for(int i_pixel_y = 0; i_pixel_y < NumOfPixel; ++i_pixel_y)
    {
      h_mTDC[i_pixel_x][i_pixel_y]->Write();
    }
  }
  File_OutPut->Close();
}

//----------------------------------
void InitDisplay_mRICH()
{
  int debug = 1;
  double var[1];
  const char * hname = H13700MAP;
  int anode, asic, pin, channel;

  //Right MPPC side (front view)
  xp_mRICH[0]=17;
  yp_mRICH[0]=0;
  xp_mRICH[1]=17;
  yp_mRICH[1]=17;

  //Left MPPC side (front view)
  xp_mRICH[2]=15;
  yp_mRICH[2]=32;
  xp_mRICH[3]=15;
  yp_mRICH[3]=15;

  FILE* fin = fopen(hname,"r");
  if(!fin) return ;
  while(fscanf(fin,"%lf",var)!=EOF){
    anode   = (int)var[0];
    fscanf(fin,"%lf",var);
    asic   = (int)var[0];
    fscanf(fin,"%lf",var);
    pin    = (int)var[0];
    int tmp;
    if(asic==2)tmp=1;
    if(asic==1)tmp=2;

    if(anode<=128){
      channel = asic*64 + pin;
    }else{
      channel = 192+asic*64 + pin;
      // channel = 191+asic*64 + pin;
    }

    if(channel<384){ maroc2h13700[channel]=anode;
      if(debug==1)printf("H13700 anode %4d  asic %2d  pin %4d  -->  ch %4d maroc %4d \n",anode, asic, pin, channel, maroc2h13700[channel]);
      // if(channel == 128) cout << "maroc2h13700 = " << maroc2h13700[channel] << endl;
    }
  }
}

//----------------------------------
void ResetEventData()
{
  tTrigNum=0;
  tTrigTime=0;
  tNedge=0;
  for(unsigned int j=0;j<MAXEDGE;j++)
  {
    tPolarity[j]=-1;
    tChannel[j]=-1;
    tTime[j]=-1;
    tSlot[j]=-1;
    tFiber[j]=-1;
  }
}

//------------------------------
int GetPMT_mRICH(int slot,int fiber,int asic)
{
  if(fiber==0 || fiber==1)return 1;
  if(fiber==2 || fiber==3)return 2;
  if(fiber==4 || fiber==5)return 3;
  if(fiber==6 || fiber==7)return 4;
}

//------------------------------
int GetPixel_mRICH(int fiber, int asic, int maroc_channel)
{
 int k=0;
 if(fiber==1 || fiber==3 || fiber==5 || fiber==7)k=1;
 int i = k*192 + asic*64 + maroc_channel;
//  if(maroc2h13700[i]==0)printf("getpixel fiber %d  asic %d ch %d  -->  ii  %d  %d \n",fiber,asic,maroc_channel,i, maroc2h13700[i]);
 return maroc2h13700[i];
}

//------------------------------
void GenCoord_mRICH(int ipmt, int x1, int y1)
{

  int j;
  int debug=1;
  int rw; // row
  int cm; // column
  for(j=0;j<256;j++)x_mRICH[j]=0;
  for(j=0;j<256;j++)y_mRICH[j]=0;

  for(j=0;j<256;j++){
    rw=(int) j/16.;
    cm=j%16;
    if(ipmt<3){
      x_mRICH[j]=x1+cm; // MPPC
      y_mRICH[j]=y1+rw;
    }else{
      x_mRICH[j]=x1-cm; // MPPC
      y_mRICH[j]=y1-rw;
    }
    // if(debug)if(j==0||j==255)printf("PMT %2d  Pixel %2d  -->  rw %3d  cm  %3d  X %3d Y %3d\n",ipmt, j+1,rw, cm,x_mRICH[j],y_mRICH[j]);
  }
}

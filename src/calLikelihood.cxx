#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <math.h> 
#include "string.h"
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TChain.h>
#include "../include/event.h"
#include "../include/hit.h"
#include "../include/material.h"
#include "../include/calLikelihood.h"

using namespace std;
// using namespace TMath;

calLikelihood::calLikelihood(string inputdatebase, string outputfile)
{
 cout<<endl;
 cout<<"calLikelihood::calLikelihood() ----- Constructor ! ------"<<endl;
 cout<<endl;
 mInPutDataBase = inputdatebase;
 mOutPutFile = outputfile;
 init();
}

calLikelihood::~calLikelihood()
{
 cout<<"calLikelihood::~calLikelihood() ----- Release memory ! ------"<<endl;
 delete mat;
 delete File_InPutDataBase;
 delete File_OutPut;
 delete mTree;

 delete hNEvtvsP;
 delete h_photonDist_piplus;
 delete h_photonDist_piminus;
 delete h_photonDist_Kplus;
 delete h_photonDist_Kminus;
 delete h_photonDist_proton;
 delete h_photonDist_antiproton;
}

int calLikelihood::init()
{
 cout<<"calLikelihood::init() ----- Initialization ! ------"<<endl;

 mat = new material(); //// initialize the material

 cout<<"calLikelihood::init(), create output file: "<< mOutPutFile.c_str() <<endl;
 File_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

 cout<<"calLikelihood::init(), read database file: "<< mInPutDataBase.c_str() <<endl;
 File_InPutDataBase = TFile::Open(mInPutDataBase.c_str());
 hNEvtvsP = (TH2D*) File_InPutDataBase->Get("hNEvtvsP");
 h_photonDist_piplus = (TH3D*)File_InPutDataBase->Get("h_photonDist_piplus");
 h_photonDist_piminus = (TH3D*)File_InPutDataBase->Get("h_photonDist_piminus");
 h_photonDist_Kplus = (TH3D*)File_InPutDataBase->Get("h_photonDist_Kplus");
 h_photonDist_Kminus = (TH3D*)File_InPutDataBase->Get("h_photonDist_Kminus");
 h_photonDist_proton = (TH3D*)File_InPutDataBase->Get("h_photonDist_proton");
 h_photonDist_antiproton = (TH3D*)File_InPutDataBase->Get("h_photonDist_antiproton");

 cout<<"calLikelihood::init(), initialize tree  ; "<<endl;
 mTree = new TTree("LLTreeDst","Tree for Likelihood Analysis");
 mTree->SetDirectory(File_OutPut);
 mTree->Branch("pid",&mDst.pid,"pid/I");
 mTree->Branch("px",&mDst.px,"px/F");
 mTree->Branch("py",&mDst.py,"py/F");
 mTree->Branch("pz",&mDst.pz,"pz/F");
 mTree->Branch("vx",&mDst.vx,"vx/F");
 mTree->Branch("vy",&mDst.vy,"vy/F");
 mTree->Branch("vz",&mDst.vz,"vz/F");
 mTree->Branch("theta",&mDst.theta,"theta/F");
 mTree->Branch("phi",&mDst.phi,"phi/F");
 mTree->Branch("nhit",&mDst.nhit,"nhit/I");
 mTree->Branch("nhitAegl",&mDst.nhitAegl,"nhitAegl/I");
 mTree->Branch("nhitPhoDet",&mDst.nhitPhoDet,"nhitPhoDet/I");
 mTree->Branch("nelSbKCs",&mDst.nelSbKCs,"nelSbKCs/I");
 mTree->Branch("nelGaAsP",&mDst.nelGaAsP,"nelGaAsP/I");
 mTree->Branch("nelGaAs",&mDst.nelGaAs,"nelGaAs/I");
 mTree->Branch("Lpion",&mDst.Lpion,"Lpion/D");
 mTree->Branch("LKaon",&mDst.LKaon,"LKaon/D");
 mTree->Branch("Lproton",&mDst.Lproton,"Lproton/D");
 mTree->SetAutoSave(5000000);

 return 0;
}

int calLikelihood::process_event(event *aevt, hit *ahit)
{
  int pid_gen=0;
  float px_gen=0;
  float py_gen=0;
  float pz_gen=0;
  float vx_gen=0;
  float vy_gen=0;
  float vz_gen=0;
  float p_gen=0;
  float theta_gen=0;
  float phi_gen=0;
  for (unsigned int i=0;i<aevt->get_pid()->size();i++) {
    pid_gen=aevt->get_pid()->at(i);
    px_gen=aevt->get_px()->at(i)/1e3;    //in MeV, convert to GeV
    py_gen=aevt->get_py()->at(i)/1e3;    //in MeV, convert to GeV
    pz_gen=aevt->get_pz()->at(i)/1e3;    //in MeV, convert to GeV
    vx_gen=aevt->get_vx()->at(i)/1e1;    //in mm, convert to cm
    vy_gen=aevt->get_vy()->at(i)/1e1;    //in mm, convert to cm
    vz_gen=aevt->get_vz()->at(i)/1e1;    //in mm, convert to cm
    p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
    theta_gen=acos(pz_gen/p_gen)*DEG;    //in deg
    phi_gen=atan2(py_gen,px_gen)*DEG;    //in deg            
  }  
  
  TH2D *h_database_pion;
  TH2D *h_database_kaon;
  TH2D *h_database_proton;
  
  if(pid_gen > 0)
  {
    int ibin = h_photonDist_piplus->GetZaxis()->FindBin(p_gen);

    h_photonDist_piplus->GetZaxis()->SetRange(ibin,ibin);
    h_database_pion = (TH2D *)h_photonDist_piplus->Project3D("xy");

    h_photonDist_Kplus->GetZaxis()->SetRange(ibin,ibin);
    h_database_kaon = (TH2D *)h_photonDist_Kplus->Project3D("xy");

    h_photonDist_proton->GetZaxis()->SetRange(ibin,ibin);
    h_database_proton = (TH2D *)h_photonDist_proton->Project3D("xy");
  }
  else if(pid_gen < 0)
  {
    int ibin = h_photonDist_piminus->GetZaxis()->FindBin(p_gen);

    h_photonDist_piminus->GetZaxis()->SetRange(ibin,ibin);
    h_database_pion = (TH2D *)h_photonDist_piplus->Project3D("xy");

    h_photonDist_Kminus->GetZaxis()->SetRange(ibin,ibin);
    h_database_kaon = (TH2D *)h_photonDist_Kplus->Project3D("xy");

    h_photonDist_antiproton->GetZaxis()->SetRange(ibin,ibin);
    h_database_proton = (TH2D *)h_photonDist_proton->Project3D("xy");
  }

  int NEvtPiBin = hNEvtvsP->FindBin(2.5,p_gen);
  int NEvtPi = hNEvtvsP->GetBinContent(NEvtPiBin);
  h_database_pion->Sumw2();
  h_database_pion->Scale(1./NEvtPi);

  int NEvtKaonBin = hNEvtvsP->FindBin(1.5,p_gen);
  int NEvtKaon = hNEvtvsP->GetBinContent(NEvtKaonBin);
  h_database_kaon->Sumw2();
  h_database_kaon->Scale(1./NEvtKaon);

  int NEvtProtonBin = hNEvtvsP->FindBin(0.5,p_gen);
  int NEvtProton = hNEvtvsP->GetBinContent(NEvtProtonBin);
  h_database_proton->Sumw2();
  h_database_proton->Scale(1./NEvtProton);
  
  TH2D *h_photonDist_PID = new TH2D("h_photonDist_PID","h_photonDist_PID",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth);
  int nhits = ahit->get_hitn()->size();
  for (int i=0;i<nhits;i++) 
  {
    if(isPhoton(ahit,i) && !isReflection(ahit,i) && isOnPhotonSensor(ahit,i))
    {
      double out_x = ahit->get_out_x()->at(i);
      double out_y = ahit->get_out_y()->at(i);
      h_photonDist_PID->Fill(out_x,out_y);
    }
  }

  mDst.pid = pid_gen;
  mDst.px = px_gen;
  mDst.py = py_gen;
  mDst.pz = pz_gen;
  mDst.vx = vx_gen;
  mDst.vy = vy_gen;
  mDst.vz = vz_gen;
  mDst.theta = theta_gen;
  mDst.phi = phi_gen;
  mDst.nhit = nhits;
  mDst.nhitAegl = 0;
  mDst.nhitPhoDet = 0;
  mDst.nelSbKCs = 0;
  mDst.nelGaAsP = 0;
  mDst.nelGaAs = 0;
  mDst.Lpion = probability(h_database_pion, h_photonDist_PID);
  mDst.LKaon = probability(h_database_kaon, h_photonDist_PID);
  mDst.Lproton = probability(h_database_proton, h_photonDist_PID);
  
  if(mTree) mTree->Fill();
  
  delete h_photonDist_PID;
  delete h_database_pion;
  delete h_database_kaon;
  delete h_database_proton;
  
  return 0;
}

int calLikelihood::end()
{
  cout<<endl;
  cout<<"calLikelihood::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(File_OutPut!= NULL){
    File_OutPut->cd();
    mTree->Write();
    File_OutPut->Close();
  }
  return 0;
}

bool calLikelihood::isPhoton(hit *ahit, int i)
{
  if(ahit->get_pid()->at(i)==0) return true;
  else return false;
}

bool calLikelihood::isReflection(hit *ahit, int i)
{
  if(ahit->get_out_pz()->at(i)<0.) return true;
  else return false;
}

bool calLikelihood::isOnAerogel(hit *ahit, int i)
{
  if(ahit->get_out_z()->at(i)>=55.5 && ahit->get_out_z()->at(i)<=85.5) return true;
  //if(ahit->get_out_z()->at(i)>=50.5 && ahit->get_out_z()->at(i)<=70.5) return true;
  else return false;
}

bool calLikelihood::isOnPhotonSensor(hit *ahit, int i)
{
  double out_z = ahit->get_out_z()->at(i);
  if(out_z > 253.6624 && out_z < 255.1626) return true;
  else return false;
}

double calLikelihood::probability(TH2D *h_database, TH2D *h_photonDist_PID)
{
  double prob=1.0;

  for(unsigned int i_x = 0; i_x < nPads; i_x++)
  {
    for(unsigned int i_y = 0; i_y < nPads; i_y++)
    {
      double k = h_photonDist_PID->GetBinContent(i_x+1,i_y+1); // detected photon number
      double lambda = h_database->GetBinContent(i_x+1,i_y+1); // averaged photon number
      double err_lambda = h_database->GetBinError(i_x+1,i_y+1); // error for specific bin
      if(err_lambda > 0.0) prob *= TMath::PoissonI(k,lambda);
    }
  }
  return log(prob);
}

////// This is the main function 
int main(int argc, char **argv)
{
  TString inputdir = "/work/eic/xusun/output/modular_rich/May02_2018/";

  TString InPutList = "/work/eic/xusun/list/modular_rich/mRICH_PID_test.list";

  string inputdatabase = "/work/eic/xusun/output/modular_rich/test/PDF_database_test.root";
  string outputfile = "/work/eic/xusun/output/modular_rich/test/PID_Likelihood_test.root";
  
  TChain *fevt  = new TChain("generated");
  TChain *fhit  = new TChain("eic_rich");

  if (!InPutList.IsNull())   // if input file is ok
  {
    TString InFo_List ="Open test file list ";
    cout << InFo_List.Data() << endl;
    ifstream in(InPutList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  addfile = inputdir+addfile;
	  fevt->AddFile(addfile.Data(),-1,"generated");
	  fhit->AddFile(addfile.Data(),-1,"eic_rich");
	  Long64_t file_entries = fevt->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      TString InFo_Warning = "WARNING: test file input is problemtic";
      cout << InFo_Warning.Data() << endl;
    }
  }
  
  event *aevt = new event(fevt);  /// declear and save info to branchs for event
  hit *ahit = new hit(fhit);	/// declear and save info to branchs for track
  
  calLikelihood *likelihood = new calLikelihood(inputdatabase,outputfile);
  
  int nevent = (int)fevt->GetEntries();
  cout << "total number of events:  " << nevent << endl;
  for (Int_t i=0;i<nevent;i++) { 
    if(i%100==0) cout << "processing event:  " << i << " ;"<<endl;
    fevt->GetEntry(i);  
    fhit->GetEntry(i);
    likelihood->process_event(aevt, ahit);
  }
  
  likelihood->end();
  delete ahit;
  delete aevt;
  delete fhit;
  delete fevt;
  
  return 0;
}

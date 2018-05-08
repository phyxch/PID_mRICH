/*=========================================================================*
 *                  Likelihood.cxx by Dr. Liang Xue                        *
 *                               04-28-2015                                *
 *-------------------------------------------------------------------------*
 * This is the main program to study the Electron Ion Collider Ring        *
 * Image Cherenkov detector.                                               *
 *                                                                         *
 * 1. read in input File -- FileInpu                                       *
 * 2. generate output histogram/tree file                                  *
 *                                                                         *
 * if "doLikelihoodDB" is defined, i.e. #define doLikelihoodDB, histograms,*
 * hNEvtvsP, hPiXYvsP, hKaonXYvsP, and hProtonXYvsP are generated.         *
 *                                                                         *
 * else hnHitAeglPerEvtvsMom, hPhotonWL, hPhotonE, hnPhotonElvsnHits_SbKCs,*
 * hnPhotonElvsnHits_GaAsP, hnPhotonElvsnHits_GaAs, and LLTreeDst (tree)   *
 * are generated.                                                          *
 *=========================================================================*
 *                            01-01-2016 Ping                              *
 *-------------------------------------------------------------------------*
 * 1. add #define doLikelihoodDB                                           *
 * 2. Assume there's no problem on making file list, comment out           *
 *    if(strstr(FileList,".root")==NULL) {}                                *
 *=========================================================================*
 *                            01-14-2016 Ping                              *
 *-------------------------------------------------------------------------*
 * change position cut to match the updated Modular RICH design (3cm       *
 * thick aerogel)                                                          *
 *=========================================================================*
 *                            01-21-2016 Ping                              *
 *-------------------------------------------------------------------------*
 * change the ranges of x and y axis of hPiXYvsP, hKaonXYvsP, hProtonXYvsP,* 
 * hXY to match the diemsion of Hamamatsu H8500 A photon sensor.          *
 * ========================================================================*/
 
#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <math.h> 
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TChain.h>
#include <TRandom.h>
#include "../include/event.h"
#include "../include/hit.h"
#include "../include/material.h"
#include "../include/likelihood.h"

//#define doLikelihoodDB

using namespace std;
using namespace TMath;


likelihood::likelihood(char *fdb, char *fout): dbf(fdb), outf(fout)
{
 cout<<endl;
 cout<<"likelihood::likelihood() ----- Constructor ! ------"<<endl;
 cout<<endl;
 init();
}

likelihood::~likelihood()
{
 cout<<"likelihood::~likelihood() ----- Release memory ! ------"<<endl;
 delete rd;
 delete mat;
 delete dbf;
 delete outf;
 delete dbfile;
 delete outputfile;
 #ifdef doLikelihoodDB
 delete hNEvtvsP;
 delete hPiXYvsP;
 delete hKaonXYvsP;
 delete hProtonXYvsP;
 #else
 delete mTree;
 delete hNEvtvsP;
 delete hPiXYvsP;
 delete hKaonXYvsP;
 delete hProtonXYvsP;
 delete hnHitAeglPerEvtvsMom;
 delete hPhotonWL;
 delete hPhotonE;
 delete hnPhotonElvsnHits_SbKCs;
 delete hnPhotonElvsnHits_GaAsP;
 delete hnPhotonElvsnHits_GaAs;
 #endif

}

int likelihood::init()
{
 cout<<"likelihood::init() ----- Initialization ! ------"<<endl;

 rd = new TRandom();    /// random number generator
 mat = new material(); //// initialize the material

 cout<<"likelihood::init(), create output file: "<< outf <<endl;
 outputfile = new TFile(outf,"RECREATE");

 #ifdef doLikelihoodDB
 cout<<"likelihood::init(), initialize database histograms ;"<<endl;
 hNEvtvsP = new TH2D("hNEvtvsP","hNEvtvsP",5,0.,5.,65,2.5,15.5);
 /*
 hPiXYvsP = new TH3D("hPiXYvsP","hPiXYvsP",nPads,-44.,44.,nPads,-44.,44.,65,2.5,15.5);
 hKaonXYvsP = new TH3D("hKaonXYvsP","hKaonXYvsP",nPads,-44.,44.,nPads,-44.,44.,65,2.5,15.5);
 hProtonXYvsP = new TH3D("hProtonXYvsP","hProtonXYvsP",nPads,-44.,44.,nPads,-44.,44.,65,2.5,15.5);
 */
 hPiXYvsP = new TH3D("hPiXYvsP","hPiXYvsP",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);
 hKaonXYvsP = new TH3D("hKaonXYvsP","hKaonXYvsP",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);
 hProtonXYvsP = new TH3D("hProtonXYvsP","hProtonXYvsP",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);
 #else 
 cout<<"likelihood::init(), read database file: "<< dbf <<endl;
 dbfile = new TFile(dbf,"READ");
 hNEvtvsP = (TH2D*) dbfile->Get("hNEvtvsP");
 hPiXYvsP = (TH3D *) dbfile->Get("hPiXYvsP");
 hKaonXYvsP = (TH3D *) dbfile->Get("hKaonXYvsP");
 hProtonXYvsP = (TH3D *) dbfile->Get("hProtonXYvsP");

 cout<<"likelihood::init(), initialize tree  ; "<<endl;
 mTree = new TTree("LLTreeDst","Tree for Likelihood Analysis");
 mTree->SetDirectory(outputfile);
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

 cout<<"likelihood::init(), initialize histogram  ; "<<endl;
 hnHitAeglPerEvtvsMom = new TH2D ("hnHitAeglPerEvtvsMom","hnHitAeglPerEvtvsMom",100,0.,10.,200,0.,200.);
 hPhotonWL = new TH1D ("hPhotonWL","hPhotonWL",100,250.,750.);
 hPhotonE = new TH1D ("hPhotonE","hPhotonE",100,0.,5.);
 hnPhotonElvsnHits_SbKCs = new TH2D ("hnPhotonElvsnHits_SbKCs","hnPhotonElvsnHits_SbKCs",50,0.,50.,25,0.,25.);
 hnPhotonElvsnHits_GaAsP = new TH2D ("hnPhotonElvsnHits_GaAsP","hnPhotonElvsnHits_GaAsP",50,0.,50.,25,0.,25.);
 hnPhotonElvsnHits_GaAs = new TH2D ("hnPhotonElvsnHits_GaAs","hnPhotonElvsnHits_GaAs",50,0.,50.,25,0.,25.);
 #endif
 
 return 0;
}

int likelihood::process_event(event *aevt, hit *ahit)
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
  // cout<< "tree generated size: "<< aevt->get_pid()->size() <<";    tree flux size:  "<< ahit->get_hitn()->size() <<endl;
  for (unsigned int i=0;i<aevt->get_pid()->size();i++) {
    // cout << aevt->get_pid()->at(i) << " " << aevt->get_px()->at(i) << " " << aevt->get_py()->at(i) << " " << aevt->get_pz()->at(i) << " " << aevt->get_vx()->at(i) << " " << aevt->get_vy()->at(i) << " " << aevt->get_vz()->at(i) << endl; 
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
  
#ifdef doLikelihoodDB
  if(pid_gen==-211) hNEvtvsP->Fill(0.,p_gen);
  else if(pid_gen==-321) hNEvtvsP->Fill(1.,p_gen);
  else if(pid_gen==2212) hNEvtvsP->Fill(2.,p_gen);
  
  int nhits = ahit->get_hitn()->size();
  for (int i=0;i<nhits;i++) {
    if(isPhoton(ahit,i) && !isReflection(ahit,i) && isOnPhotonSensor(ahit,i)){
      double out_x(0.), out_y(0.);
      Smearing2D(ahit->get_out_x()->at(i), ahit->get_out_y()->at(i), out_x, out_y);
      if(fabs(out_x)>halfWidth||fabs(out_y)>halfWidth) continue;
      double photonE = ahit->get_trackE()->at(i);   /// in MeV (GEANT4 default)
      double wavelength = 1240./(photonE*1.e6);  /// MeV->eV,wavelength in "nm"
      double QE_GaAsP=mat->extrapQE_GaAsP(wavelength);
      if(QE_GaAsP>rd->Uniform(0.,1.)){
	if(pid_gen==-211) hPiXYvsP->Fill(out_x,out_y,p_gen);
	else if(pid_gen==-321) hKaonXYvsP->Fill(out_x,out_y,p_gen);
	else if(pid_gen==2212) hProtonXYvsP->Fill(out_x,out_y,p_gen);
      }
    }
  }
#else
  TH2D *hPiXYDB, *hKaonXYDB, *hProtonXYDB;     //// get database at particular momentum for pi, K, p
  
  int ibin = hPiXYvsP->GetZaxis()->FindBin(p_gen);
  hPiXYvsP->GetZaxis()->SetRange(ibin,ibin);
  hKaonXYvsP->GetZaxis()->SetRange(ibin,ibin);
  hProtonXYvsP->GetZaxis()->SetRange(ibin,ibin);
  hPiXYDB = (TH2D *) hPiXYvsP->Project3D("xy");
  hKaonXYDB = (TH2D *) hKaonXYvsP->Project3D("xy");
  hProtonXYDB = (TH2D *) hProtonXYvsP->Project3D("xy");
  
  int NEvtPiBin = hNEvtvsP->FindBin(0,p_gen);
  int NEvtKaonBin = hNEvtvsP->FindBin(1,p_gen);
  int NEvtProtonBin = hNEvtvsP->FindBin(2,p_gen);
  int NEvtPi = hNEvtvsP->GetBinContent(NEvtPiBin);
  int NEvtKaon = hNEvtvsP->GetBinContent(NEvtKaonBin);
  int NEvtProton = hNEvtvsP->GetBinContent(NEvtProtonBin);
  hPiXYDB->Sumw2();
  hKaonXYDB->Sumw2();
  hProtonXYDB->Sumw2();
  hPiXYDB->Scale(1./NEvtPi);
  hKaonXYDB->Scale(1./NEvtKaon);
  hProtonXYDB->Scale(1./NEvtProton);
  
  //TH2D *hXY = new TH2D("hXY","hXY",nPads,-44.,44.,nPads,-44.,44.);  /// get current XY distribution
  TH2D *hXY = new TH2D("hXY","hXY",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth);  /// get current XY distribution
  double nhitPhoDet(0.0), nhitAerogel(0.0), nel_SbKCs(0.0), nel_GaAsP(0.0), nel_GaAs(0.0);
  int nhits = ahit->get_hitn()->size();
  
  for (int i=0;i<nhits;i++) {
    if(isPhoton(ahit,i) && !isReflection(ahit,i) && isOnAerogel(ahit,i)){
      nhitAerogel++;
      double photonE = ahit->get_trackE()->at(i);   /// in MeV (GEANT4 default)
      double wavelength = 1240./(photonE*1.e6);  /// MeV->eV,wavelength in "nm"
      hPhotonWL->Fill(wavelength);
      hPhotonE->Fill(photonE*1.e6);
    }
    
    if(isPhoton(ahit,i) && !isReflection(ahit,i) && isOnPhotonSensor(ahit,i)){
      nhitPhoDet++;
      double out_x(0.), out_y(0.);
      Smearing2D(ahit->get_out_x()->at(i), ahit->get_out_y()->at(i), out_x, out_y);
      if(fabs(out_x)>halfWidth||fabs(out_y)>halfWidth) continue; //hit outside the photonsensor --> ignore
      
      double photonE = ahit->get_trackE()->at(i);   /// in MeV (GEANT4 default)
      double wavelength = 1240./(photonE*1.e6);  /// MeV->eV,wavelength in "nm"
      
      double QE_SbKCs=mat->extrapQE_SbKCs(wavelength);
      if(QE_SbKCs>rd->Uniform(0.,1.)){
	nel_SbKCs++;
      }
      
      double QE_GaAsP=mat->extrapQE_GaAsP(wavelength);
      if(QE_GaAsP>rd->Uniform(0.,1.)){
	nel_GaAsP++; 
	hXY->Fill(out_x,out_y);
      }
      
      double QE_GaAs=mat->extrapQE_GaAs(wavelength);
      if(QE_GaAs>rd->Uniform(0.,1.)){
	nel_GaAs++;
      }
    }
  }
  
  ///// randomly throw 2 electrons in each event as noise
  for(int i=0; i<2; i++) {
    double out_x = rd->Uniform(-1.,1.)*halfWidth;
    double out_y = rd->Uniform(-1.,1.)*halfWidth;
    hXY->Fill(out_x,out_y);
  }
  
  hnHitAeglPerEvtvsMom->Fill(p_gen,nhitAerogel);
  hnPhotonElvsnHits_SbKCs->Fill(nhitPhoDet,nel_SbKCs);
  hnPhotonElvsnHits_GaAsP->Fill(nhitPhoDet,nel_GaAsP);
  hnPhotonElvsnHits_GaAs->Fill(nhitPhoDet,nel_GaAs);
  
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
  mDst.nhitAegl = nhitAerogel;
  mDst.nhitPhoDet = nhitPhoDet;
  mDst.nelSbKCs = nel_SbKCs;
  mDst.nelGaAsP = nel_GaAsP;
  mDst.nelGaAs = nel_GaAs;
  mDst.Lpion = probability(hPiXYDB, hXY);
  mDst.LKaon = probability(hKaonXYDB, hXY);
  mDst.Lproton = probability(hProtonXYDB, hXY);
  
  if(mTree) mTree->Fill();
  
  delete hXY;
  delete hPiXYDB;
  delete hKaonXYDB;
  delete hProtonXYDB;
#endif
  
  return 0;
}

int likelihood::end()
{
  cout<<endl;
  cout<<"likelihood::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(outputfile != NULL){
    outputfile->cd();
#ifdef doLikelihoodDB
    hNEvtvsP->Write();
    hPiXYvsP->Write();
    hKaonXYvsP->Write();
    hProtonXYvsP->Write();
#else
    mTree->Write();
    hnHitAeglPerEvtvsMom->Write();
    hPhotonWL->Write();
    hPhotonE->Write();
    hnPhotonElvsnHits_SbKCs->Write();
    hnPhotonElvsnHits_GaAsP->Write();
    hnPhotonElvsnHits_GaAs->Write();
#endif
    outputfile->Close();
  }
  return 0;
}

bool likelihood::isPhoton(hit *ahit, int i)
{
  if(ahit->get_pid()->at(i)==0) return true;
  else return false;
}

bool likelihood::isReflection(hit *ahit, int i)
{
  if(ahit->get_out_pz()->at(i)<0.) return true;
  else return false;
}

bool likelihood::isOnAerogel(hit *ahit, int i)
{
  if(ahit->get_out_z()->at(i)>=55.5 && ahit->get_out_z()->at(i)<=85.5) return true;
  //if(ahit->get_out_z()->at(i)>=50.5 && ahit->get_out_z()->at(i)<=70.5) return true;
  else return false;
}

bool likelihood::isOnPhotonSensor(hit *ahit, int i)
{
  if(ahit->get_out_z()->at(i)>=159.5) return true;
  //if(ahit->get_out_z()->at(i)>143 && ahit->get_out_z()->at(i)<147) return true;
  else return false;
}

double likelihood::probability(TH2D *db, TH2D *hXY)
{
  double prob=1.0;
  const double noise = 2./nPads/nPads; ///// noise level is set to be 2 electrons on photonsensor
  
  for(unsigned int i=1; i<=nPads; i++){   ///// photonsensor pad loop in X
    for(unsigned int j=1; j<=nPads; j++)   ///// photonsensot pad loop in Y
      {
	double k = hXY->GetBinContent(i,j); ///// should be integer (mostly 0 or 1)
	double lambda = db->GetBinContent(i,j); ///// get the expected value of Poisson distribution from DB
	if(lambda==0.) prob *= PoissonI(k,noise);
	else if(lambda>0.) prob *= PoissonI(k,lambda);
	else cout<<"DB histogram bin content should not be < 0. !"<<endl;
      }
  }
  return log(prob);
}

void likelihood::Smearing2D(double inx, double iny, double& outx, double& outy)
{
  TF1 *fx = new TF1("fx","1/([1]*sqrt(2*3.1415926))*exp(-1*pow(x-[0],2)/(2.*[1]*[1]))",inx-20.,inx+20.);
  fx->SetParameter(0,inx);
  fx->SetParameter(1,1.); //// 1 mm smearing in photon position x
  outx = fx->GetRandom();
  
  TF1 *fy = new TF1("fy","1/([1]*sqrt(2*3.1415926))*exp(-1*pow(x-[0],2)/(2.*[1]*[1]))",iny-20.,iny+20.);
  fy->SetParameter(0,iny);
  fy->SetParameter(1,1.); //// 1 mm smearing in photon position y
  outy = fy->GetRandom();
  
  delete fx;
  delete fy;
  return ;
}


////// This is the main function 
int main(int argc, char **argv)
{
  cout<<endl;
  cout<<"//////////////////////////////////////////////"<<endl;
  cout<<"This program start at: "<<endl;
  time_t now = time(0);
  char *dt = ctime(&now);
  cout<<dt<<endl;
  char currdir[256];
  cout<<"Current directory is: "<< getcwd(currdir,256)  <<endl;
  cout<<"//////////////////////////////////////////////"<<endl;
  cout<<endl;
  
  if(argc!=3 && argc!=4) return 0;
  
  char *FileInput=0;
  char *FileDB=0;
  char *FileOutput=0;
  
  if(argc==3){
    FileInput = argv[1];
    FileOutput = argv[2];
  }
  else if(argc==4){
    FileInput = argv[1];
    FileDB = argv[2];
    FileOutput = argv[3];
  }
  
  char dbf[256], outputf[256];
  if(FileDB) sprintf(dbf,"%s",FileDB);
  sprintf(outputf,"%s",FileOutput);
  
  TChain *fevt  = new TChain("generated");
  TChain *fhit  = new TChain("eic_rich");
  
  int fileNumber = 0;
  char FileList[512];
  ifstream* inputStream = new ifstream;
  inputStream->open(FileInput);
  if (!(inputStream)) {
    printf("can not open list file\n");
    return 0;
  }
  for (;inputStream->good();) {
    inputStream->getline(FileList,512);
    
    /*--------------------------------------------------------------
      01-01-2016 Ping
      this block gives me error somehow.
      Assume there's no problem on making file list, I comment out
      this block. Now it is working fine in "doLikelihoodDB" macros
      
      if(strstr(FileList,".root")==NULL) {
      printf("%s is not a root-file address!!!\n",FileList);
      continue;
      }
      ----------------------------------------------------------------*/
    if  ( inputStream->good()) {
      TFile *ftmp = new TFile(FileList);
      if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
	printf(" file %s error in opening!!!\n",FileList);
      }
      else {
	printf(" read in file %s\n",FileList);
	fevt->Add(FileList);
	fhit->Add(FileList);
	fileNumber++;
      }
      delete ftmp;
    }
  }
  delete inputStream;
  printf(" files read in %d\n",fileNumber);
  
  event *aevt = new event(fevt);  /// declear and save info to branchs for event
  hit *ahit = new hit(fhit);	/// declear and save info to branchs for track
  
  likelihood *lk = new likelihood(dbf,outputf);
  
  int nevent = (int)fevt->GetEntries();
  cout << "total number of events:  " << nevent << endl;
  for (Int_t i=0;i<nevent;i++) { 
    if(i%100==0) cout << "processing event:  " << i << " ;"<<endl;
    fevt->GetEntry(i);  
    fhit->GetEntry(i);
    lk->process_event(aevt, ahit);
  }
  
  lk->end();
  delete ahit;
  delete aevt;
  delete fhit;
  delete fevt;
  
  return 0;
}

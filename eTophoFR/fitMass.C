#include <iostream>
#include "TMath.h"
#include "TTree.h"
#include "TGraphMultiErrors.h"
//#include "TLeafF16.h"
//#include "TLeafF16.h"
/*
#ifdef __CINT__
#pragma link C++ class TLeafF16+
#endif
*/

void fitNominalMass(RooDataSet *data, RooRealVar x, RooCategory cat, string name, double &eff, double &err, string plotSuf, double xmin, double xmax){


  
  RooRealVar efficiency("efficiency","efficiency",0.9,0,1);

  //voigtian parameters
  //RooRealVar mup("mup","mup",90,70,130);
  RooRealVar mup("mup","mup",86,70,130);
  RooRealVar widthp("widthp","widthp",2.495,2.495,2.495); 
  //RooRealVar sigmap("sigmap","sigmap",5,0,20); //sigma f the resolution function which is Gaus here
  RooRealVar sigmap("sigmap","sigmap",2,0.01,5); //sigma f the resolution function which is Gaus here
  

  RooRealVar muf("muf","muf",90,70,130);
  RooRealVar widthf("widthf","widthf",2.495,2.495,2.495); 
  //RooRealVar sigmaf("sigmaf","sigmaf",5,0,20); //sigma f the resolution function which is Gaus here
  RooRealVar sigmaf("sigmaf","sigmaf",2,0.01,5); //sigma f the resolution function which is Gaus here
 
  //convolute BW with gaussian - Voigt
  //RooBreitWigner signal("signal","signal",x,mu,sigma);
  RooVoigtian signalp("signalp","signalp",x,mup,widthp,sigmap);  

  RooVoigtian signalf("signalf","signalf",x,muf,widthf,sigmaf);  

  RooRealVar alphap("alphap","alphap",-0.01,-10,10);
  RooRealVar alphaf("alphaf","alphaf",-0.08,-10,10);
  RooExponential backgroundp("backgroundp","backgroundp",x,alphap);
  RooExponential backgroundf("backgroundf","backgroundf",x,alphaf);

 
  RooRealVar nsigtot("nsigtot","nsigtot",10000,0,1e+8);
  RooFormulaVar nsigp("nsigp","nsigp","efficiency*nsigtot",RooArgList(efficiency,nsigtot));
  RooFormulaVar nsigf("nsigf","nsigf","(1-efficiency)*nsigtot",RooArgList(efficiency,nsigtot));
  
  RooRealVar nbkgp("nbkgp","nbkgp",100,0,1e+5);
  RooRealVar nbkgf("nbkgf","nbkgf",1000,0,1e+5);
  
  
  RooAddPdf totalp("totalp","totalp",RooArgList(signalp,backgroundp),RooArgList(nsigp,nbkgp));
  RooAddPdf totalf("totalf","totalf",RooArgList(signalf,backgroundf),RooArgList(nsigf,nbkgf));



  RooSimultaneous model("model","",cat);
  model.addPdf(totalp,"passed");
  model.addPdf(totalf,"failed");

  cout<<"Reading tree"<<endl;

  
  
  RooFitResult *results = model.fitTo(*data,RooFit::Extended(true),RooFit::Save(true));


  RooRealVar* effPar = (RooRealVar*) results->floatParsFinal().find("efficiency");
  
  //overwritten below by createIntegral 
  eff = effPar->getVal();
  err = effPar->getError();
  
  RooPlot *framep = x.frame(RooFit::Title("Passed"));
  data->plotOn(framep,RooFit::Cut("cat==cat::passed"));
  model.plotOn(framep,RooFit::Slice(cat,"passed"),RooFit::ProjWData(RooArgSet(cat),*data));
  model.plotOn(framep,RooFit::Slice(cat,"passed"),RooFit::ProjWData(RooArgSet(cat),*data),RooFit::Components("backgroundp"),RooFit::LineColor(kRed));
  

  RooPlot *framef = x.frame(RooFit::Title("Failed"));
  data->plotOn(framef,RooFit::Cut("cat==cat::failed"));
  model.plotOn(framef,RooFit::Slice(cat,"failed"),RooFit::ProjWData(RooArgSet(cat),*data));
  model.plotOn(framef,RooFit::Slice(cat,"failed"),RooFit::ProjWData(RooArgSet(cat),*data),RooFit::Components("backgroundf"),RooFit::LineColor(kRed));
  //model.paramOn(framef, RooFit::Layout(0.55));
  //model.paramOn(framef, RooFit::Layout(0.35));

  RooArgSet *args = new RooArgSet(efficiency,mup,sigmap,muf,sigmaf,nsigp,nsigf);
  // X size of box is from 65% to 99% of Xaxis range, top of box is at 90% of Yaxis range)
  model.paramOn(framef, RooFit::Layout(0.65, 0.99, 0.9), RooFit::Parameters(*args));
  //framef->getAttText()->SetTextSize(0.01);

  //model.paramOn(framef);
  
  

  TCanvas c("c","c",900,500);
  c.Divide(2);
  c.cd(1);
  framep->Draw();
  c.cd(2);
  framef->Draw();
  
  c.Draw();
  c.Print(Form("plots/plot_eta_%s_%s.png",plotSuf.c_str(),name.c_str()));

  ////efficiency within range
  x.setRange("Zee",xmin, xmax);
  //x.setRange("Zee",60,120);
  //RooAbsReal* fsigregion_model = model.createIntegral(x,NormSet(x),Range("Zee"),RooFit::Slice(cat,"passed"),RooFit::Components("signalp")); //The "NormSet(x)" normalizes it to the total number of events to give you the fraction n_signal_region_events/n_total_events


  RooAbsReal* fsigregionPass_model = totalp.createIntegral(x,RooFit::NormSet(x),RooFit::Range("Zee")); 
  RooAbsReal* fsigregionFail_model = totalf.createIntegral(x,RooFit::NormSet(x),RooFit::Range("Zee")); 

  RooAbsReal* fbkgPass_model = backgroundp.createIntegral(x,RooFit::NormSet(x),RooFit::Range("Zee")); 
  RooAbsReal* fbkgFail_model = backgroundf.createIntegral(x,RooFit::NormSet(x),RooFit::Range("Zee")); 
  
  double nsigpass = fsigregionPass_model->getVal()*(nsigp.getVal() + nbkgp.getVal()) - fbkgPass_model->getVal()*nbkgp.getVal();
  double nsigfail = fsigregionFail_model->getVal()*(nsigf.getVal() + nbkgf.getVal()) - fbkgFail_model->getVal()*nbkgf.getVal();
  double eff_integ = nsigpass/(nsigpass+nsigfail);
  cout<<"Actual fitted sigp and sigf "<<nsigp.getVal()<<" "<<nsigf.getVal()<<endl;
  cout<<"Pass : Fail using integ  "<<nsigpass<<" "<<nsigfail<<endl;
  cout<<"Eff comparison : efff : eff_integ "<<eff<<" "<<eff_integ<<endl;


  //eff = effPar->getVal();
  //err = effPar->getError();

  eff = eff_integ;
  err = effPar->getError();

  //return eff;
}



////alternate signal model
void fitAltBkgMass(RooDataSet *data, RooRealVar x, RooCategory cat, string name, double &eff, double &err){


  
  RooRealVar efficiency("efficiency","efficiency",0.9,0,1);

  //voigtian parameters
   RooRealVar mup("mup","mup",90,70,130);
  RooRealVar widthp("widthp","widthp",2.495,2.495,2.495); 
  RooRealVar sigmap("sigmap","sigmap",5,0,20); //sigma f the resolution function which is Gaus here


  RooRealVar muf("muf","muf",90,70,130);
  RooRealVar widthf("widthf","widthf",2.495,2.495,2.495); 
  RooRealVar sigmaf("sigmaf","sigmaf",5,0,20); //sigma f the resolution function which is Gaus here
 
  //convolute BW with gaussian - Voigt
  //RooBreitWigner signal("signal","signal",x,mu,sigma);
  RooVoigtian signalp("signalp","signalp",x,mup,widthp,sigmap);  

  RooVoigtian signalf("signalf","signalf",x,muf,widthf,sigmaf);  

  RooRealVar coef0p("c0p","coefficient #0",1.0,-1.,1);
  RooRealVar coef1p("c1p","coefficient #1",0.1,-1.,1);
  RooRealVar coef2p("c2p","coefficient #2",-0.1,-1.,1);

  RooRealVar coef0f("c0f","coefficient #0",1.0,-1.,1);
  RooRealVar coef1f("c1f","coefficient #1",0.1,-1.,1);
  RooRealVar coef2f("c2f","coefficient #2",-0.1,-1.,1);
  
  RooChebychev backgroundp("backgroundp","backgroundp",x,RooArgList(coef0p, coef1p, coef2p));
  RooChebychev backgroundf("backgroundf","backgroundf",x,RooArgList(coef0f, coef1f, coef2f));

 
  RooRealVar nsigtot("nsigtot","nsigtot",10000,0,1e+8);
  RooFormulaVar nsigp("nsigp","nsigp","efficiency*nsigtot",RooArgList(efficiency,nsigtot));
  RooFormulaVar nsigf("nsigf","nsigf","(1-efficiency)*nsigtot",RooArgList(efficiency,nsigtot));
  
  RooRealVar nbkgp("nbkgp","nbkgp",100,0,1e+5);
  RooRealVar nbkgf("nbkgf","nbkgf",1000,0,1e+5);
  
  
  RooAddPdf totalp("totalp","totalp",RooArgList(signalp,backgroundp),RooArgList(nsigp,nbkgp));
  RooAddPdf totalf("totalf","totalf",RooArgList(signalf,backgroundf),RooArgList(nsigf,nbkgf));



  RooSimultaneous model("model","",cat);
  model.addPdf(totalp,"passed");
  model.addPdf(totalf,"failed");

  cout<<"Reading tree"<<endl;

  
  
  RooFitResult *results = model.fitTo(*data,RooFit::Extended(true),RooFit::Save(true));


  RooRealVar* effPar = (RooRealVar*) results->floatParsFinal().find("efficiency");
  
  eff = effPar->getVal();
  err = effPar->getError();
  
  RooPlot *framep = x.frame(RooFit::Title("Passed"));
  data->plotOn(framep,RooFit::Cut("cat==cat::passed"));
  model.plotOn(framep,RooFit::Slice(cat,"passed"),RooFit::ProjWData(RooArgSet(cat),*data),RooFit::LineColor(kGreen));
  model.plotOn(framep,RooFit::Slice(cat,"passed"),RooFit::ProjWData(RooArgSet(cat),*data),RooFit::Components("backgroundp"),RooFit::LineColor(kRed));
  

  RooPlot *framef = x.frame(RooFit::Title("Failed"));
  data->plotOn(framef,RooFit::Cut("cat==cat::failed"));
  model.plotOn(framef,RooFit::Slice(cat,"failed"),RooFit::ProjWData(RooArgSet(cat),*data),RooFit::LineColor(kGreen));
  model.plotOn(framef,RooFit::Slice(cat,"failed"),RooFit::ProjWData(RooArgSet(cat),*data),RooFit::Components("backgroundf"),RooFit::LineColor(kRed));
  //model.paramOn(framef, RooFit::Layout(0.55));
  //model.paramOn(framef, RooFit::Layout(0.35));
  RooArgSet *args = new RooArgSet(efficiency,mup,sigmap,muf,sigmaf,nsigp,nsigf);

  // X size of box is from 65% to 99% of Xaxis range, top of box is at 90% of Yaxis range)
    model.paramOn(framef, RooFit::Layout(0.65, 0.99, 0.9), RooFit::Parameters(*args));
  //framef->getAttText()->SetTextSize(0.01);
  
  //model.paramOn(framef);
  
  

  TCanvas c("c","c",900,500);
  c.Divide(2);
  c.cd(1);
  framep->Draw();
  c.cd(2);
  framef->Draw();
  
  c.Draw();
  c.Print(Form("plots/plot_eta_altbkg_%s.png",name.c_str()));

  //return eff;
}



void callfitMass(){

  //invariant x
  RooRealVar x("selPairMass","selPairMass",60,120);
  RooRealVar xRed("selPairMass","selPairMass",80,100);
  RooRealVar pt("selphoPt","selphoPt",0,5000);
  RooRealVar eta("selphoEta","selphoEta",-10,10);
  RooRealVar phi("selphoPhi","selphoPhi",-10,10);
  RooRealVar haspix("selphoHasPix","selphoHasPix",-10,10);


  RooCategory cat("cat","cat");
  cat.defineType("passed");
  cat.defineType("failed");
  
  ///get data from the rootfile
  //TFile *fin = TFile::Open("minitree_data.root");
  TFile *fin = TFile::Open("minitree_data95_out.root");
  //TFile *fin = TFile::Open("minitree_dataB_out.root");
  TTree *tin = (TTree*)fin->Get("minitree");
  
  
  cout<<"Now importing data"<<endl;
  
  //RooDataSet *dataFull = new RooDataSet("dataFull","All data",tin,RooArgSet(x,pt,eta,phi,haspix));
  //RooDataSet *dataFull = new RooDataSet("dataFull","All data",tin,RooArgSet(x,pt,eta,phi));
  RooDataSet *dataFull = new RooDataSet("dataFull","All data",RooArgSet(x,pt,eta,phi,haspix),RooFit::Import(*tin));
  
  
  cout<<"Imported data"<<endl;

  /*
  const int ncats = 6;
  double fr[ncats], fr_err[ncats];
  string cut[ncats] = {"fabs(selphoEta)>0 && fabs(selphoEta)<0.1", "fabs(selphoEta)>0.1 && fabs(selphoEta)<0.4", "fabs(selphoEta)>0.4 && fabs(selphoEta)<0.8", "fabs(selphoEta)>0.8 && fabs(selphoEta)<1.0","fabs(selphoEta)>1 && fabs(selphoEta)<1.4442", "fabs(selphoEta)>1.566 && fabs(selphoEta)<2.5"};
  string name[ncats] = {"EB1", "EB2", "EB3", "EB4", "EB5", "EE1"};

  double etaArr[ncats] = {0.1, 0.4,0.8,1,1.4442,2.5};
  */


  string commonCut = "selphoPt > 200" ; ////11Nov, 2020 - after studying pT dependence, beyond 200, its quite flat

  ///just in the EB - abs bin in eta
  /*const int ncats = 5;
  double fr[ncats], fr_err[ncats];
  double fr_altbkg[ncats], fr_altbkg_err[ncats], fr_altmass[ncats], fr_altmass_err[ncats], fr_sys_err[ncats];
  double fr_tot_err[ncats]; 
  double etalerr[ncats], etaherr[ncats];  

  string cut[ncats] = {"fabs(selphoEta)>0 && fabs(selphoEta)<0.1", "fabs(selphoEta)>0.1 && fabs(selphoEta)<0.4", "fabs(selphoEta)>0.4 && fabs(selphoEta)<0.8", "fabs(selphoEta)>0.8 && fabs(selphoEta)<1.0","fabs(selphoEta)>1 && fabs(selphoEta)<1.4442"};
  string name[ncats] = {"EB1", "EB2", "EB3", "EB4", "EB5"};
  double etaArr[ncats] = {0.1, 0.4,0.8,1,1.4442};
  
  int plusEtaBin = 0;
  */

  /*
  const int ncats = 6;
  double fr[ncats], fr_err[ncats];
  double fr_altbkg[ncats], fr_altbkg_err[ncats], fr_altmass[ncats], fr_altmass_err[ncats], fr_sys_err[ncats];
  double fr_tot_err[ncats]; 
  double etalerr[ncats], etaherr[ncats];  

  string cut[ncats] = {"fabs(selphoEta)>0 && fabs(selphoEta)<0.1", "fabs(selphoEta)>0.1 && fabs(selphoEta)<0.2", "fabs(selphoEta)>0.2 && fabs(selphoEta)<0.4", "fabs(selphoEta)>0.4 && fabs(selphoEta)<0.6", "fabs(selphoEta)>0.6 && fabs(selphoEta)<1", "fabs(selphoEta)>1 && fabs(selphoEta)<1.4442"};
  string name[ncats] = {"EB1", "EB2", "EB3", "EB4", "EB5","EB6"};
  double etaArr[ncats] = {0.1, 0.2, 0.4,0.6,1,1.4442};
  
  int plusEtaBin = 0;
  */

  /*
  const int ncats = 8;
  double fr[ncats], fr_err[ncats];
  double fr_altbkg[ncats], fr_altbkg_err[ncats], fr_altmass[ncats], fr_altmass_err[ncats], fr_sys_err[ncats];
  double fr_tot_err[ncats]; 
  double etalerr[ncats], etaherr[ncats];  

  string cut[ncats] = {"fabs(selphoEta)>0 && fabs(selphoEta)<0.05", "fabs(selphoEta)>0.05 && fabs(selphoEta)<0.1", "fabs(selphoEta)>0.1 && fabs(selphoEta)<0.2", "fabs(selphoEta)>0.2 && fabs(selphoEta)<0.4", "fabs(selphoEta)>0.4 && fabs(selphoEta)<0.6", "fabs(selphoEta)>0.6 && fabs(selphoEta)<0.8", "fabs(selphoEta)>0.8 && fabs(selphoEta)<1", "fabs(selphoEta)>1 && fabs(selphoEta)<1.4442"};
  string name[ncats] = {"EB1", "EB2", "EB3", "EB4", "EB5","EB6", "EB7", "EB8"};
  double etaArr[ncats] = {0.05, 0.1, 0.2, 0.4,0.6,0.8,1,1.4442};
  
  int plusEtaBin = 0;
  */

  const int ncats = 5;
  double fr[ncats], fr_err[ncats];
  double fr_altbkg[ncats], fr_altbkg_err[ncats], fr_altmass[ncats], fr_altmass_err[ncats], fr_sys_err[ncats];
  double fr_tot_err[ncats]; 
  double etalerr[ncats], etaherr[ncats];  

  string cut[ncats] = {"fabs(selphoEta)>0 && fabs(selphoEta)<0.1", "fabs(selphoEta)>0.1 && fabs(selphoEta)<0.4", "fabs(selphoEta)>0.4 && fabs(selphoEta)<0.8", "fabs(selphoEta)>0.8 && fabs(selphoEta)<1.1", "fabs(selphoEta)>1.1 && fabs(selphoEta)<1.4442"};
  string name[ncats] = {"EB1", "EB2", "EB3", "EB4", "EB5"};
  double etaArr[ncats] = {0.1, 0.4,0.8,1.1,1.4442};
  
  int plusEtaBin = 0;

  ///just in the EB - abs bin in eta
  /*
  const int ncats = 10;
  double fr[ncats], fr_err[ncats];
  double fr_altbkg[ncats], fr_altbkg_err[ncats], fr_altmass[ncats], fr_altmass_err[ncats], fr_sys_err[ncats];
  double fr_tot_err[ncats];
  double etalerr[ncats], etaherr[ncats];


  string cut[ncats] = {"selphoEta>-1.4442 && selphoEta<-1", "selphoEta>-1 && selphoEta<-0.8", "selphoEta>-0.8 && selphoEta<-0.4", "selphoEta>-0.4 && selphoEta<-0.1", "selphoEta>-0.1 && selphoEta<0", "selphoEta>=0 && selphoEta<0.1", "selphoEta>0.1 && selphoEta<0.4", "selphoEta>0.4 && selphoEta<0.8", "selphoEta>0.8 && selphoEta<1.0","selphoEta>1 && selphoEta<1.4442"};
  string name[ncats] = {"EB1", "EB2", "EB3", "EB4", "EB5", "EB6", "EB7", "EB8", "EB9", "EB10"};
  double etaArr[ncats] = {-1.4444, -1, -0.8, -0.4, -0.1, 0.1, 0.4,0.8,1,1.4442};
  
  int plusEtaBin = 5;
  */

  /*
  const int ncats = 5;
  double fr[ncats], fr_err[ncats];
  double fr_altbkg[ncats], fr_altbkg_err[ncats], fr_altmass[ncats], fr_altmass_err[ncats], fr_sys_err[ncats];
  double fr_tot_err[ncats];
  double etalerr[ncats], etaherr[ncats];
  
  string cut[ncats] = {"fabs(selphoEta)<1.4442 && selphoPt>=70&&selphoPt<100", "fabs(selphoEta)<1.4442 && selphoPt>=100&&selphoPt<150","fabs(selphoEta)<1.4442 && selphoPt>=150&&selphoPt<200",
		       "fabs(selphoEta)<1.4442 && selphoPt>=200&&selphoPt<300","fabs(selphoEta)<1.4442 && selphoPt>=300&&selphoPt<500"};
  
  string name[ncats] = {"pt1", "pt2", "pt3", "pt4", "pt5"};
  double etaArr[ncats] = {100, 150, 200, 300, 500};
  
  int plusEtaBin = 5;
  */


  for(int icat=0; icat<ncats; icat++){
    
    RooDataSet *datap = (RooDataSet*)dataFull->reduce(RooArgSet(x,pt,eta,phi,haspix),Form("%s && selphoHasPix && %s",cut[icat].c_str(),commonCut.c_str()));
    
    //RooDataSet *datap = (RooDataSet*)dataFull->reduce(RooArgSet(x,pt,eta,phi,haspix));
    datap->Print();
    RooDataSet *dataf = (RooDataSet*)dataFull->reduce(RooArgSet(x,pt,eta,phi,haspix),Form("%s && !selphoHasPix && %s",cut[icat].c_str(), commonCut.c_str()));
    dataf->Print();
    
    cout<<"============================================"<<endl;
    cout<<"icat "<<icat<<endl;
    cout<<"Entries in passing and failing "<< datap->numEntries()<<" "<<dataf->numEntries()<<endl;

    ///combined data
    RooDataSet *data = new RooDataSet("data","data",RooArgSet(x),RooFit::Index(cat),RooFit::Import("passed",*datap),RooFit::Import("failed",*dataf));
    
    ///Nominal 
    double eff = 0; double err = 0;
    fitNominalMass(data, x, cat, name[icat], eff, err, "nominal", 60, 120);

    fr[icat] = (1-eff)/eff;
    fr_err[icat] = err;

    //////////////////Alt bkg
    double effAltBkg = 0; double errAltBkg = 0;
    fitAltBkgMass(data, x, cat, name[icat], effAltBkg, errAltBkg);
    double fr_altbkg = (1-effAltBkg)/effAltBkg;
    fr_altbkg_err[icat] = fabs(fr[icat]-fr_altbkg);
    
    ///////////Alternate mass window of 80 to 100
    /*RooDataSet *datapAltMass = (RooDataSet*)dataFull->reduce(RooArgSet(x,pt,eta,phi,haspix),Form("%s && selphoHasPix && selPairMass>80 && selPairMass < 100",cut[icat].c_str()));
    datapAltMass->Print();
    RooDataSet *datafAltMass = (RooDataSet*)dataFull->reduce(RooArgSet(x,pt,eta,phi,haspix),Form("%s && !selphoHasPix && selPairMass>80 && selPairMass < 100",cut[icat].c_str()));
    datafAltMass->Print();
    
    cout<<"============================================"<<endl;
    cout<<"icat "<<icat<<endl;
    cout<<"Alt Mass : Entries in passing and failing "<< datapAltMass->numEntries()<<" "<<datafAltMass->numEntries()<<endl;

    ///combined data
    double effAltMass = 0; double errAltMass = 0;
    RooDataSet *dataAltMass = new RooDataSet("dataAltMass","dataAltMass",RooArgSet(xRed),RooFit::Index(cat),RooFit::Import("passed",*datapAltMass),RooFit::Import("failed",*datafAltMass));
    fitNominalMass(dataAltMass, xRed, cat, name[icat], effAltMass, errAltMass, "altMass");
    double fr_altmass = (1-effAltMass)/effAltMass;
    fr_altmass_err[icat] = fr_altmass;
    */
    double effAltMass = 0; double errAltMass = 0;
    fitNominalMass(data, x, cat, name[icat], effAltMass, errAltMass, "altMass", 80, 100);
    double fr_altmass = (1-effAltMass)/effAltMass;
    fr_altmass_err[icat] = fabs(fr_altmass - fr[icat]);

    /////////Now full combination
    fr_tot_err[icat] = sqrt( pow( pow(fr_err[icat],2) + fr_altbkg_err[icat],2) + pow(fr_altmass_err[icat],2) );

    
    fr_sys_err[icat] = sqrt( pow( fr_altbkg_err[icat],2) + pow(fr_altmass_err[icat],2) );



    /////////Now full combination
    fr_tot_err[icat] = sqrt( pow( pow(fr_err[icat],2) + fr_altbkg_err[icat],2) + pow(fr_altmass_err[icat],2) );

    
    fr_sys_err[icat] = sqrt( pow( fr_altbkg_err[icat],2) + pow(fr_altmass_err[icat],2) );

    if(icat!=0) etalerr[icat] = fabs(etaArr[icat-1] - etaArr[icat])/2.;
    else etalerr[icat] = 0;
    
    if(icat!=(ncats-1)) 
      etaherr[icat] = fabs(etaArr[icat+1] - etaArr[icat])/2.;
    else etaherr[icat] = 0;

    cout<<"================fr  : stat : Altbkg : altmass : "<<fr[icat]<<" "<<fr_altbkg_err[icat]<<" "<<fr_altmass_err[icat]<<endl;
  }//for(int icat=0; icat<ncats; icat++)  



  

  /*TGraphErrors *g = new TGraphErrors(ncats,etaArr,fr, etaerr, fr_err);
  g->SetMarkerStyle(20);
  g->Draw("APE");
  */
  
  TFile *fout = new TFile("histFR.root","recreate");

  ///store simple graph with total errors
  TCanvas c("c","c",900,500);
  TGraphAsymmErrors *gs = new TGraphAsymmErrors(ncats,etaArr,fr, etalerr, etaherr, fr_tot_err, fr_tot_err);
  gs->SetMarkerStyle(20);
  gs->Draw("APE");
  c.Draw();
  c.Print("plots/plot_eta_fr_simple.png");
  c.Print("plots/plot_eta_fr_simple.C");


  
  TCanvas c1("c1","c1",900,500);
  
  TGraphMultiErrors* g = new TGraphMultiErrors("g", "", ncats, etaArr, fr, etalerr, etaherr, fr_err, fr_err);
  g->GetXaxis()->SetTitle("#eta");
  g->GetYaxis()->SetTitle("e#rightarrow#gamma fake  rate");
  
  g->AddYError(ncats, fr_sys_err, fr_sys_err);
  //g->AddYError(ncats, fr_tot_err, fr_tot_err);
  g->SetMarkerStyle(20);
  g->SetLineColor(1);
  g->GetAttLine(0)->SetLineColor(8);
  g->GetAttLine(0)->SetLineWidth(3);
  g->GetAttLine(1)->SetLineColor(9);
  g->GetAttFill(1)->SetFillColor(9);

  //g->GetAttLine(2)->SetLineColor(46);
  //g->GetAttFill(2)->SetFillColor(46);

  // Graph and x erros drawn with "APS"
  // Stat Errors drawn with "Z"
  // Sys Errors drawn with "5 s=0.5"
  //gme->Draw("APS ; Z ; 5 s=0.5");
  g->Draw("APS ; ; 5 s=0.5"); //Z:Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
  //g->Draw("APS ; ; 5 s=0.5 ; 5 s=0.5"); //Z:Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
  //g->Draw("P ; ; 5 s=0.5 same" );
  g->Draw("P ;  same");


  c1.Draw();
  c1.Print("plots/plot_eta_fr.png");
  c1.Print("plots/plot_eta_fr.C");


  //////-ve
  TCanvas cM("cM","cM",900,500);
  TGraphMultiErrors* gM = new TGraphMultiErrors("gM", "", ncats/2, etaArr, fr, etalerr, etaherr, fr_err, fr_err);
  gM->GetXaxis()->SetTitle("#eta");
  gM->GetYaxis()->SetTitle("e#rightarrow#gamma fake  rate");
  
  gM->AddYError(ncats/2, fr_altbkg_err, fr_altbkg_err);
  gM->SetMarkerStyle(20);
  gM->SetLineColor(1);
  gM->GetAttLine(0)->SetLineColor(8);
  gM->GetAttLine(0)->SetLineWidth(3);
  gM->GetAttLine(1)->SetLineColor(9);
  gM->GetAttFill(1)->SetFillColor(9);
  // Graph and x erros drawn with "APS"
  // Stat Errors drawn with "Z"
  // Sys Errors drawn with "5 s=0.5"
  //gme->Draw("APS ; Z ; 5 s=0.5");
  gM->Draw("APS ; ; 5 s=0.5"); //Z:Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
  gM->Draw("P ;  same");


  cM.Draw();
  cM.Print("plots/plot_fr_etaM.png");
  cM.Print("plots/plot_fr_etaM.pdf");
  cM.Print("plots/plot_fr_etaM.C");

  ////+ve
  TCanvas cP("cP","cP",900,500);
  TGraphMultiErrors* gP = new TGraphMultiErrors("gP", "", ncats/2, &etaArr[plusEtaBin], &fr[plusEtaBin], &etalerr[plusEtaBin], &etaherr[plusEtaBin], &fr_err[plusEtaBin], &fr_err[plusEtaBin]);
  gP->GetXaxis()->SetTitle("#eta");
  gP->GetYaxis()->SetTitle("e#rightarrow#gamma fake  rate");
  
  gP->AddYError(ncats/2, &fr_altbkg_err[plusEtaBin], &fr_altbkg_err[plusEtaBin]);
  gP->SetMarkerStyle(20);
  gP->SetLineColor(1);
  gP->GetAttLine(0)->SetLineColor(8);
  gP->GetAttLine(0)->SetLineWidth(3);
  gP->GetAttLine(1)->SetLineColor(9);
  gP->GetAttFill(1)->SetFillColor(9);
  // Graph and x erros drawn with "APS"
  // Stat Errors drawn with "Z"
  // Sys Errors drawn with "5 s=0.5"
  //gme->Draw("APS ; Z ; 5 s=0.5");
  gP->Draw("APS ; ; 5 s=0.5"); //Z:Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
  gP->Draw("P ;  same");


  cP.Draw();
  cP.Print("plots/plot_fr_etaP.png");
  cP.Print("plots/plot_fr_etaP.pdf");
  cP.Print("plots/plot_fr_etaP.C");



  
  fout->cd();
  g->Write();
  gM->Write();
  gP->Write();
  fout->Write();
  fout->Close();

  for(int i=0; i<g->GetN();i++){
    cout<<"Stat : sys "<<g->GetErrorY(i,0)<<" "<<g->GetErrorY(i,1)<<endl;
  }

}//void callfitMass()
  

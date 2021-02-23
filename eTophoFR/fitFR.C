#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>
#include "TString.h"
#include <TH2.h>
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TSystem.h"
//#include "TProfile"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include <iostream>
using namespace std;


TH1D *hpr;
TH1D *hpi;
TH1D *hdata;
TGraphMultiErrors* gr;

double chiSqArr[10000];
int icall = 0;
double parValues[4][1000000]; ///for 3 parameters

TMatrixD xval_gr;
TMatrixD data_gr;
TMatrixD data_gr_T;
TMatrixD data_err_gr;

bool fitLin = false;
//bool fitLin = true;


Double_t fitLinFun(Double_t *x, Double_t *par){

  double consta = par[0];
  double slope = par[1];

  double f = consta + slope * x[0];
  return f;
}

Double_t fitQuadFun(Double_t *x, Double_t *par){

  double consta = par[0];
  double slope = par[1];
  double quad = par[2];

  double f = consta + slope * x[0] + quad * x[0] * x[0];
  return f;
}


void getTGraphMatrix(){
  
  cout<<"Inside getTGraphMatrix"<<endl;
  const Int_t nbins = gr->GetN();

  TMatrixD xval_mat(1,nbins);
  TMatrixD data_mat(1,nbins);
  TArrayD data_arr(nbins);
  TArrayD xval_arr(nbins);
  TMatrixD data_mat_tr(nbins,1);
  
  for (int i=0;i<nbins; i++) {
    double xval, yval;
    gr->GetPoint(i,xval,yval);
    data_arr[i] = yval;
    xval_arr[i] = xval;
  }


  xval_mat.SetMatrixArray(xval_arr.GetArray());
  data_mat.SetMatrixArray(data_arr.GetArray());
  

  data_mat_tr.Transpose(data_mat);

  //xval_mat.Print();
  
  //data_mat.Print();
  //data_mat_tr.Print();
  
  ////now error matrix
  TMatrixD data_err_mat(nbins,nbins);
  TArrayD data_err_arr(nbins*nbins);
  int ntot = 0;
  for (int i=0;i<nbins; i++) {
    for (int j=0;j<nbins; j++) {
      
      ///both stat and sys for diagonal
      if(i==j) data_err_arr[ntot] = pow(gr->GetErrorY(i,0),2) + pow(gr->GetErrorY(i,1),2) ;
      else data_err_arr[ntot] = gr->GetErrorY(i,1)*gr->GetErrorY(j,1); //only sys
      //else data_err_arr[ntot] = 0; //only sys
      ntot++;
    }
  }

  data_err_mat.SetMatrixArray(data_err_arr.GetArray());
  //data_err_mat.Print();

  ////now copy all the matrices here
  xval_gr.ResizeTo(xval_mat);
  xval_gr=xval_mat;

  data_gr.ResizeTo(data_mat);
  data_gr=data_mat;
  
  data_gr_T.ResizeTo(data_mat_tr);
  data_gr_T=data_mat_tr;
  
  data_err_gr.ResizeTo(data_err_mat);
  data_err_gr=data_err_mat;
  

  //data_err_gr.Print();


}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  //cout<<"===================CALL=================="<<endl;
  //getTGraphMatrix();
  //data_err_gr.Print();

  const Int_t nbins = xval_gr.GetNoElements();
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  
  TMatrixD data_mat(1,nbins);
  TArrayD data_arr(nbins);
  TMatrixD data_mat_tr(nbins,1);


  double chi2_simple = 0;
  
  for (int i=0;i<nbins; i++) {
    double xval[1] = {xval_gr[0][i]}; ///make sure its not the transpose - use the 1xnbins matrix 
    
    
    double hypVal = -999;
    if(fitLin) 
      hypVal = fitLinFun(xval,par);

    if(!fitLin)
      hypVal = fitQuadFun(xval,par);

    ////now form the FR - hyp array
    data_arr[i] = data_gr[0][i] - hypVal;///use the 1xnbins matrix

    chi2_simple += data_arr[i]*data_arr[i]/data_err_gr[i][i]; //already squard
    //cout<<"i : x : y : yhyp :  err : dif : "<<i <<" "<<xval[0]<<" "<<data_gr[0][i]<<" "<<hypVal<<" "<<data_err_gr[i][i]<<" "<<data_arr[i]<<endl;
  }

  data_mat.SetMatrixArray(data_arr.GetArray());

  //cout<<"Printing diff matrix"<<endl;
  //data_mat.Print();

  data_mat_tr.Transpose(data_mat);

  
  //cout<<"Printing transpose"<<endl;
  //data_mat_tr.Print();
  
  ///makee a copy of this error matrix because .Invert() inverts the actual matrix and heence in the next call, it takes the inverted one
  TMatrixD data_err_gr_inv = data_err_gr;
  data_err_gr_inv.Invert();
  
  //cout<<"Priting the inverted error matrix "<<endl;
  //data_err_gr_inv.Print();
  
  ////now get the chi2
  TMatrix c(nbins,1); ///
  c.Mult(data_err_gr_inv,data_mat_tr);
  
  //cout<<"Print intemediary matrix"<<endl;
  //c.Print();
  
  
  TMatrix chi2(1,1); ///
  chi2.Mult(data_mat,c);
  

  f = chi2[0][0];

  //if(fitLin) cout<<"par0 : par1 : chi2 : chi2_simple : "<<par[0]<<" "<<par[1]<<" "<<f<<" "<<chi2_simple<<endl;
  //if(!fitLin) cout<<"par0 : par1 : par 2 : chi2 : chi2_simple : "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<f<<" "<<chi2_simple<<endl;
  chiSqArr[icall] = chi2[0][0];
  parValues[0][icall] = par[0];
  parValues[1][icall] = par[1];
  if(!fitLin)   parValues[2][icall] = par[2];

  icall = icall+1;
  
  
  //std::cout<<"chisq "<<chisq<<endl;
}


/////////////////////////Gett and Generate multivariate gaussian/////////////////////////
////https://root.cern.ch/doc/v606/MultivariateGaussianTest_8C_source.html
RooDataSet*  getMVG(RooArgList xVec, RooArgList muVec, TMatrixDSym covMat, int ngen){

  RooMultiVarGaussian mvg("mvg", "mvg", xVec, muVec, covMat);
 
  // --------------------
  // make a toy dataset
  RooDataSet* data = mvg.generate(xVec, ngen);
 
  data->Print();
  /*
  ///print the events
  for(Int_t i=0; i < data->numEntries(); i++)
    {
      double x0 = (data->get(i))->getRealValue("x0");
      double x1 = (data->get(i))->getRealValue("x1");
      double x2 = (data->get(i))->getRealValue("x2");

      //cout<<"SJ!!!!x0 : x1 : x2 "<<x0<<" "<<x1 <<" "<<x2<<endl;
    }
  */
  
  return data;

}

map<string, TH1F*> generateMVG(RooDataSet *data, TF1 *f1, TGraphMultiErrors *gr, TFile *fout, int fitLin){
  
  int nGpoints = gr->GetN();
  int nPars = f1->GetNpar();

  fout->cd();
  map<string, TH1F*> hmap;

  for(int ipar=0; ipar<nPars; ipar++){
    
    double xval = f1->GetParameter(ipar);
    double xmin = -xval*3;
    double xmax = xval*3;

    if(xmax<xmin){
      xmax = xmin;
      xmin = -xmax;
    }

    hmap[Form("par%d",ipar)] = new TH1F(Form("par%d",ipar), Form("par%d",ipar), 500, xmin, xmax);
  }

  for(int ix=0; ix<nGpoints; ix++){
    
    hmap[Form("xval%d",ix)] = new TH1F(Form("xval%d",ix), Form("xval%d",ix), 500, -0.1, 0.2);
  }

  //TF1 *fsys = (TF1*)f1->Clone(); ///somehow this doesnt work. it doesnt change fsys
  
  TF1 *fsys = new TF1("f",fitLinFun,0,1.5,2);
  if(!fitLin){
    fsys = new TF1("f",fitQuadFun,0,1.5,3);
  }

  for(Int_t i=0; i < data->numEntries(); i++){
    
    double x0 = (data->get(i))->getRealValue("x0");
    double x1 = (data->get(i))->getRealValue("x1");
    
    double x2;
    
    if(!fitLin) x2 = (data->get(i))->getRealValue("x2");
    
    hmap[Form("par0")]->Fill(x0);
    hmap[Form("par1")]->Fill(x1);
    if(!fitLin) hmap[Form("par2")]->Fill(x2);

    
    fsys->SetParameters(x0, x1);
    
    if(!fitLin){

      
      fsys->SetParameters(x0, x1, x2);

      //cout<<"x0 : x1 : x2 : FR sys "<<x0<<" "<<x1<<" "<<x2<<" "<<fsys->Eval(1.)<<endl;
    }
    
    
    for(int ip=0; ip<nGpoints; ip++){
      double x,y;
      gr->GetPoint(ip,x,y);
      double fr = fsys->Eval(x);

      hmap[Form("xval%d",ip)]->Fill(fr);

      //cout<<"FR for sys: f1 : frsys at "<<ip<<" is "<<f1->Eval(x)<<" "<<fr<<endl;
      
    }//for(int ip=0; ip->GetN(); ip++)
      
  }//for(Int_t i=0; i < data->numEntries(); i++) 

  fout->cd();
  for(map<string,TH1F*>::iterator it = hmap.begin(); it != hmap.end(); ++it) {
    hmap[it->first]->Write();
  }
  
  return hmap;

}

///////////////////////////////////////////////////////////////////////
void fitFR(){

  TFile *fin = TFile::Open("histFR.root");
  //gr = (TGraphMultiErrors*)fin->Get("gM");
  //gr = (TGraphMultiErrors*)fin->Get("gP");
  gr = (TGraphMultiErrors*)fin->Get("g");

  ///fout
  TFile *fout = new TFile("histo_out.root","RECREATE");
  
  const Int_t nbins = gr->GetN();

  getTGraphMatrix();

  /*
  data_gr.Print();

  data_err_gr.Print();
  
  ///Mulltiply
  TMatrix c(nbins,1); ///
  c.Mult(data_err_gr,data_gr_T);
  
  TMatrix d(1,1); ///
  d.Mult(data_gr,c);

  c.Print();
  d.Print();
  */

  int nparsToFit = 2;
  if(!fitLin) nparsToFit = 3;

  TMinuit *gMinuit = new TMinuit(nparsToFit);  //initialize TMinuit with a maximum of 2 params
  
  gMinuit->SetFCN(fcn);
  
  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg); ///http://polywww.in2p3.fr/activites/physique/glast/workbook/pages/advanced_GRUG/GRUGminuit.pdf

  // Set starting values and step sizes for parameters
  //static Double_t vstart[2] = {0.12,0.05};
  //static Double_t step[2] = {0.01,0.01};
  //gMinuit->mnparm(0, "Const", vstart[0], step[0], -0.3,0.3,ierflg);
  //gMinuit->mnparm(1, "Slope", vstart[1], step[1], -0.5,0.5,ierflg);

  //static Double_t vstart[2] = {0.12,-0.054};
  if(fitLin){
    static Double_t vstart[2] = {5.84055e-02,-2.37529e-02}; ///abs eta bin = 5
    static Double_t step[2] = {0.01,0.01};
    
    gMinuit->mnparm(0, "Const", vstart[0], step[0], -0.3,0.3,ierflg);
    gMinuit->mnparm(1, "Slope", vstart[1], step[1], -0.5,0.5,ierflg);
  }


  if(!fitLin){
    static Double_t vstart[3] = {5.27120e-02, -2.49742e-03, -1.45032e-02}; ///abs eta bin = 5
    static Double_t step[3] = {0.01,0.001, 0.01};
    /*
    gMinuit->mnparm(0, "Const", vstart[0], step[0], -0.6,0.6,ierflg);
    gMinuit->mnparm(1, "Slope", vstart[1], step[1], -0.5,0.5,ierflg);
    gMinuit->mnparm(2, "Quad", vstart[2], step[2], -0.5,0.5,ierflg);
    */

    gMinuit->mnparm(0, "Const", vstart[0], step[0], -0.2,0.2,ierflg);
    gMinuit->mnparm(1, "Slope", vstart[1], step[1], -0.01,0.01,ierflg);
    gMinuit->mnparm(2, "Quad", vstart[2], step[2], -0.2,0.2,ierflg);

  }


  ///Minimization strategy
  // 1 standard
  // 2 try to improve minimum (slower)
  arglist[0]=2;
  gMinuit->mnexcm("SET STR",arglist,1,ierflg);
  



  // Now ready for minimization step
  //arglist[0] = 500;
  arglist[0] = 5000;  ///no. of calls to migrad
  arglist[1] = 1;

  // Run the simplex minimizer to get close to the minimum
  gMinuit->mnexcm("SIMPLEX",arglist ,2,ierflg);
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  //gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);

  ////plot the likelihood scan
  //arglist[0] = 1;  // TMinuit starts from 1 
  arglist[0] = 2;  // TMinuit starts from 1 
  arglist[1] = 30;
  gMinuit->mnexcm("SCAN",arglist,2,ierflg);

  // get TGraph object
  TGraph * grLL = dynamic_cast<TGraph *>(gMinuit->GetPlot() );

  /*
  TCanvas cL("cL","cL",900,500);
  gr->Draw("AL");
  cL.Draw();
  cL.Print("plots/plot_fr_fit_LikeP1.png");
  cL.Print("plots/plot_fr_fit_LikeP1.C");
  */


  /// Print results
  cout << "\nPrint results from minuit\n";
  double fParamVal[3];
  double fParamErr[3];
  gMinuit->GetParameter(0,fParamVal[0],fParamErr[0]);
  cout << "a1=" << fParamVal[0] << "\n";
  cout << "a1Err=" << fParamErr[0] << "\n";


  gMinuit->GetParameter(1,fParamVal[1],fParamErr[1]);
  cout << "a2=" << fParamVal[1] << "\n";
  cout << "a2Err=" << fParamErr[1] << "\n";

  if(!fitLin){
    gMinuit->GetParameter(2,fParamVal[2],fParamErr[2]);
    cout << "a3=" << fParamVal[2] << "\n";
    cout << "a3Err=" << fParamErr[2] << "\n";
  }

  int dim = 2;
  if(!fitLin) dim = 3;

  TMatrixDSym covMat(dim);
  
  TArrayD covMatArray(dim*dim);

  ///get th cov matrix
  if(fitLin){
    int ntot = 0;
    Double_t matrix[2][2];
    gMinuit->mnemat(&matrix[0][0],2);

    for(int i=0; i<2; i++){
      cout<<"|";
      for(int j=0; j<2; j++){
	covMatArray[ntot] = matrix[i][j];
	cout<<" "<<matrix[i][j];
	ntot++;
      }//for(int i=0; i<2; i++)
      cout<<"|"<<endl;
    }//for(int j=0; j<2; j++)

  }

  if(!fitLin){

    int ntot = 0;
    Double_t matrix[3][3];
    gMinuit->mnemat(&matrix[0][0],3);

    for(int i=0; i<3; i++){
      cout<<"|";
      for(int j=0; j<3; j++){
	covMatArray[ntot] = matrix[i][j];
	cout<<" "<<matrix[i][j];
	ntot++;
      }//for(int i=0; i<3; i++)
      cout<<"|"<<endl;
    }//for(int j=0; j<3; j++)

  }

  covMat.SetMatrixArray(covMatArray.GetArray());  

  ////get MV gaussian

  RooArgList xVec;
  RooArgList muVec;

  
  // make the observable and means
  
  RooRealVar* x;
  RooRealVar* mu_x;
  for (int i = 0; i < dim; i++) {
    char* name = Form("x%d", i);
    x = new RooRealVar(name, name, 0, -10,10);
    xVec.add(*x);
    
    char* mu_name = Form("mu_x%d",i);
    mu_x = new RooRealVar(mu_name, mu_name, fParamVal[i], -10,10);
    muVec.add(*mu_x);

  }

  RooDataSet* data = getMVG(xVec,muVec,covMat,10000);
  /////end of getting MVG

  //TF1 *f1 = new TF1("f",fitLinFun,-1.5,0,2);
  TF1 *f1 = new TF1("f",fitLinFun,0,1.5,2);
  f1->SetParameters(fParamVal[0], fParamVal[1]);
  
  if(!fitLin){

    //f1 = new TF1("f",fitLinFun,0,1.5,3);
    f1 = new TF1("f",fitQuadFun,0,1.5,3);
    f1->SetParameters(fParamVal[0], fParamVal[1], fParamVal[2]);
    //0.04813347  0.0074153  -0.02120253 --->taken by doing a fit: scipy.optimize.curve_fit.
    //f1->SetParameters(0.04813347, 0.0074153,  -0.02120253);

    //0.05247852 -0.00197266 -0.01522909 ---> uncorrrelated scipy
    //    f1->SetParameters(0.05247852, -0.00197266, -0.01522909);
  }///if(!fitLin)

  ////get the systematics the band to this.
  int nGpoints = gr->GetN();
  TGraph *grHighErr = new TGraph(nGpoints);
  TGraph *grLowErr = new TGraph(nGpoints);

  map<string, TH1F*> hmap = generateMVG(data, f1, gr, fout, fitLin);
  
  
  for(int ip=0; ip<nGpoints; ip++){
    double x,y;
    gr->GetPoint(ip,x,y);  
    double highErr = f1->Eval(x) + hmap[Form("xval%d",ip)]->GetRMS();
    double lowErr = f1->Eval(x) - hmap[Form("xval%d",ip)]->GetRMS();
    
    grHighErr->SetPoint(ip, x, highErr);
    grLowErr->SetPoint(ip, x, lowErr);
  }

  ////Clone does not work: https://root-forum.cern.ch/t/tobject-clone-and-tf1/8085 so use copy constructor
  TF1 *fhighErr = new TF1("fHE",fitLinFun,0,1.5,2);
  TF1 *flowErr = new TF1("fLE",fitLinFun,0,1.5,2);
  if(!fitLin){
    fhighErr = new TF1("fHE",fitQuadFun,0,1.5,3);
    flowErr = new TF1("fHE",fitQuadFun,0,1.5,3);
  }

  ///use a simplle TF1->Fit() this time to fit
  grHighErr->Fit(fhighErr,"RIVE");
  grLowErr->Fit(flowErr,"RIVE");

  grHighErr->SetMarkerColor(2);
  grHighErr->SetMarkerStyle(20);
  grLowErr->SetMarkerColor(4);
  grLowErr->SetMarkerStyle(20);

  fhighErr->SetLineStyle(9);
  fhighErr->SetLineColor(7);

  flowErr->SetLineStyle(9);
  flowErr->SetLineColor(7);
  
  ////get the fit values as low and high band
  double etaVals[100], fr_highErr[100], fr_lowErr[100], xerr[100], fr_med[100];
  for(int ip=0; ip<nGpoints; ip++){
    double x,y;
    gr->GetPoint(ip,x,y);
    double highErr = fhighErr->Eval(x) - f1->Eval(x) ;
    double lowErr = -flowErr->Eval(x) + f1->Eval(x) ;
    
    etaVals[ip] = x;
    fr_highErr[ip] = highErr;
    fr_lowErr[ip] = lowErr;
    xerr[ip] = 0;
    fr_med[ip] = f1->Eval(x);
  }

  TGraphAsymmErrors *gSys = new TGraphAsymmErrors(nGpoints,etaVals,fr_med, xerr,xerr, fr_lowErr,fr_highErr);
  gSys->SetFillColor(32);
  gSys->SetFillStyle(3005);
  gSys->SetLineColor(32);
  //g1->SetMarkerStyle(22);

  
  ///end of getting systematic band

  
  TCanvas c("c","c",900,500);
  gr->SetMarkerStyle(20);
  gr->SetLineColor(1);
  gr->GetAttLine(0)->SetLineColor(8);
  gr->GetAttLine(0)->SetLineWidth(3);
  gr->GetAttLine(1)->SetLineColor(9);
  gr->GetAttFill(1)->SetFillColor(9);


  ///draw the one with shade
  gr->Draw("APS ; ; 5 s=0.5");
  gr->Draw("P ;  same");
  f1->Draw("same");
  
  ///Sys band
  ///for checks of the exact fits on the low and high band
  /*grHighErr->Draw("Psame");
  grLowErr->Draw("Psame");
  fhighErr->Draw("same");
  flowErr->Draw("same");
  */
  ///not drawing it before because then the axis is set by TGraph 
  gSys->Draw("3C");  
  gr->Draw("PS ; ; 5 s=0.5");
  gr->Draw("P ;  same");
  f1->Draw("same");
  

  ///

  TLegend *leg = new TLegend(0.6080402,0.7125436,0.8994975,0.8954704,NULL,"brNDC");
  leg->AddEntry(gr,"Measured","lfp");
  leg->AddEntry(f1,"Fit","l");
  leg->AddEntry(gSys,"68% Unc. band","f");
  leg->Draw();


  string suf = "lin";
  if(!fitLin) suf = "quad";

  c.Draw();
  c.Print(Form("plots/plot_fr_fit_%s.png",suf.c_str()));
  c.Print(Form("plots/plot_fr_fit_%s.C",suf.c_str()));


  ///////////Finally printing the median values and the low and high Sys bands
  double fParamVal_HE[3], fParamVal_LE[3];
  fParamVal_HE[0] = fhighErr->GetParameter(0);
  fParamVal_HE[1] = fhighErr->GetParameter(1);

  fParamVal_LE[0] = flowErr->GetParameter(0);
  fParamVal_LE[1] = flowErr->GetParameter(1);


  if(!fitLin){

    fParamVal_HE[2] = fhighErr->GetParameter(2);
    fParamVal_LE[2] = flowErr->GetParameter(2);
  }


  if(fitLin){

    cout<<"Median fit eqn : "<<fParamVal[0]<<" + "<<fParamVal[1]<<"|x|"<<endl;
    cout<<"HighSys fit eqn : "<<fParamVal_HE[0]<<" + "<<fParamVal_HE[1]<<"|x|"<<endl;
    cout<<"LowSys fit eqn : "<<fParamVal_LE[0]<<" + "<<fParamVal_LE[1]<<"|x|"<<endl;
  }

  if(!fitLin){

    cout<<"Median fit eqn : "<<fParamVal[0]<<" + "<<fParamVal[1]<<"|x| + "<<fParamVal[2]<<"x*x"<<endl;
    cout<<"HighSys fit eqn : "<<fParamVal_HE[0]<<" + "<<fParamVal_HE[1]<<"|x| + "<<fParamVal_HE[2]<<"x*x"<<endl;
    cout<<"LowSys fit eqn : "<<fParamVal_LE[0]<<" + "<<fParamVal_LE[1]<<"|x| + "<<fParamVal_LE[2]<<"x*x"<<endl;
  }

  
  //////////////////////////////////////////////////////////////
  ///////////////////get the correlation///////////////////////
  gMinuit->SetErrorDef(1); ///SetErrorDef(N2) for N-sigma error
  //TGraph *gr0 = (TGraph *)gMinuit->Contour(80,0,1); ///Show 80 points on contour, show correlation of param 0 with param 1
  TGraph *gr0 = (TGraph *)gMinuit->Contour(200,0,1); ///Show 80 points on contour, show correlation of param 0 with param 1
  gr0->SetLineColor(kRed);
  
  c.Clear();
  gr0->SetMarkerStyle(20);
  gr0->GetXaxis()->SetTitle("parameter 0");
  gr0->GetYaxis()->SetTitle("parameter 1");

  gr0->Draw("alp");
  c.Draw();
  c.Print(Form("plots/plot_corr_0vs1_%s.png",suf.c_str()));
  c.Print(Form("plots/plot_corr_0vs1_%s.C",suf.c_str()));

  ///between 1 &2
  gr0 = (TGraph *)gMinuit->Contour(200,1,2); ///Show 80 points on contour, show correlation of param 0 with param 1
  gr0->SetLineColor(kRed);
  
  c.Clear();
  gr0->SetMarkerStyle(20);
  gr0->GetXaxis()->SetTitle("parameter 1");
  gr0->GetYaxis()->SetTitle("parameter 2");

  gr0->Draw("alp");
  c.Draw();
  c.Print(Form("plots/plot_corr_0vs1_%s.png",suf.c_str()));
  c.Print(Form("plots/plot_corr_0vs1_%s.C",suf.c_str()));


  ///between 0 &2
  gr0 = (TGraph *)gMinuit->Contour(200,0,2); ///Show 80 points on contour, show correlation of param 0 with param 1
  gr0->SetLineColor(kRed);
  
  c.Clear();
  gr0->SetMarkerStyle(20);
  gr0->GetXaxis()->SetTitle("parameter 0");
  gr0->GetYaxis()->SetTitle("parameter 2");

  gr0->Draw("alp");
  c.Draw();
  c.Print(Form("plots/plot_corr_0vs2_%s.png",suf.c_str()));
  c.Print(Form("plots/plot_corr_0vs2_%s.C",suf.c_str()));

  ///Likelihood scan 

  ////par 0
  c.Clear();
  c.cd();
  gr0 = new TGraph(icall,parValues[0],chiSqArr);
  gr0->SetMarkerStyle(20);
  gr0->Draw("ap");

  c.Draw();
  c.Print(Form("plots/plot_LLscan_0_%s.png",suf.c_str()));
  c.Print(Form("plots/plot_LLscan_0_%s.C",suf.c_str()));

  ////par 1
  c.Clear();
  c.cd();
  gr0 = new TGraph(icall,parValues[1],chiSqArr);
  gr0->SetMarkerStyle(20);
  gr0->Draw("ap");

  c.Draw();
  c.Print(Form("plots/plot_LLscan_1_%s.png",suf.c_str()));
  c.Print(Form("plots/plot_LLscan_1_%s.C",suf.c_str()));


  ////par 2
  c.Clear();
  c.cd();
  gr0 = new TGraph(icall,parValues[2],chiSqArr);
  gr0->SetMarkerStyle(20);
  gr0->Draw("ap");

  c.Draw();
  c.Print(Form("plots/plot_LLscan_2_%s.png",suf.c_str()));
  c.Print(Form("plots/plot_LLscan_2_%s.C",suf.c_str()));

  /*
  ///Likelihood scan
  arglist[0] = 0;
  gMinuit->mnexcm("SCAN", arglist
		  ,1,ierflg);
  
  */

  // Print results
  //Double_t amin,edm,errdef;
  //Int_t nvpar,nparx,icstat;
  //gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  

}

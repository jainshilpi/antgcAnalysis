#define Analyse_cxx
#include "Analyse.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"

#include <iostream>

#include "Utilities.h"


//bool debug = true;
bool debug = false;

void Analyse::Loop(string fnamein, string fnameout, bool isData, double xsec, int itype, string idType)
{
  
  using namespace std;
  

  
  if(debug) cout<<"inside Loop"<<endl;
  
  cout<<"File name, xsec "<<fnamein<<" "<<xsec<<endl;

  
  bool debugL1 = false;
  //bool debugL1 = true;



  const double pi = 4.*atan(1.);

  if(debug){
    
    cout<<"opening file now"<<endl;
    }


  TFile *fin = TFile::Open(fnamein.c_str());
  if(debug){
    
    cout<<"opened file"<<endl;
    }

  
  TTree *tin = (TTree*)fin->Get("Z2eeTree");
  if(debug){
    
    cout<<"openedtree"<<endl;
  }

  Init(tin);
  if(debug){
    
    cout<<"initiated tree"<<endl;
  }
  

  Int_t selphoHasPix, selphoIsMatchEle, selTotPairs;
  Float_t selphoElePt, selphoEleEta, selphoElePhi, selphoEleSCEta, selphoEleSCPhi;
  Float_t selphoPt, selphoEta, selphoPhi, selPairMass, selphoSCEta, selphoSCPhi, selelePt, seleleEta, selelePhi, seleleSCEta, seleleSCPhi;
  
  TFile *fmintree = new TFile( ("minitree_"+fnameout).c_str(),"RECREATE");
  TTree *tmintree = new TTree("minitree","Tree for limit calculation");
  tmintree->Branch("run", &run);
  tmintree->Branch("lumi", &lumis);
  tmintree->Branch("event", &event, "event/L");
  tmintree->Branch("selphoPt", &selphoPt, "selphoPt/F");
  tmintree->Branch("selphoEta", &selphoEta, "selphoEta/F");
  tmintree->Branch("selphoPhi", &selphoPhi, "selphoPhi/F");
  tmintree->Branch("selphoSCEta", &selphoSCEta, "selphoSCEta/F");
  tmintree->Branch("selphoSCPhi", &selphoSCPhi, "selphoSCPhi/F");
  tmintree->Branch("selphoHasPix", &selphoHasPix, "selphoHasPix/I");

  tmintree->Branch("selphoElePt", &selphoElePt, "selphoElePt/F");
  tmintree->Branch("selphoEleEta", &selphoEleEta, "selphoEleEta/F");
  tmintree->Branch("selphoElePhi", &selphoElePhi, "selphoElePhi/F");
  tmintree->Branch("selphoEleSCEta", &selphoEleSCEta, "selphoEleSCEta/F");
  tmintree->Branch("selphoEleSCPhi", &selphoEleSCPhi, "selphoEleSCPhi/F");

  tmintree->Branch("selTotPairs", &selTotPairs, "selTotPairs/I");

  tmintree->Branch("selphoIsMatchEle", &selphoIsMatchEle, "selphoIsMatchEle/I");


  tmintree->Branch("selelePt", &selelePt, "selelePt/F");
  tmintree->Branch("seleleEta", &seleleEta, "seleleEta/F");
  tmintree->Branch("selelePhi", &selelePhi, "selelePhi/F");
  tmintree->Branch("seleleSCEta", &seleleSCEta, "seleleSCEta/F");
  tmintree->Branch("seleleSCPhi", &seleleSCPhi, "seleleSCPhi/F");




  tmintree->Branch("selPairMass", &selPairMass, "selPairMass/F");

  map<string,TH1F*>hmap;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<100000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    bool isdoublephoHLT70 = HLTPho>>23&1;
    if(!isdoublephoHLT70)
      continue;
    
    ///to be away from the signal region
    if(pfMET > 150)
      continue;
    
    
    //select electron first
    int nGoodEle = 0;
    int eleArr[100];
    int* eleIndArr = selectEle(nGoodEle, eleArr);


    int nGoodPho = 0;
    int phoArr[100];
    int* phoIndArr = selectPho(nGoodPho, phoArr, idType);

    if(debug || debugL1) cout<<""<<endl;    
    if(debug || debugL1) cout<<"======"<<endl;

    if(debug || debugL1) cout<<"good ele : good pho "<<nGoodEle<<" "<<nGoodPho<<endl;
    
    //int elePhoInd[100];
    //int usedPho[100];
    
    //std::pair <int,int> elePhoPair;
    
    vector<std::pair <int,int>> elePhoPair;
    
    vector<int> usedPho;
    for(int iele=0; iele<nGoodEle; iele++){

      int eleInd = *(eleIndArr + iele);
      
      double tmpmassDiff = 999;
      bool foundPair = false;
      int tmpelePhoInd = -99;
      
      if(debug || debugL1) cout<<"iele : pT "<<iele<<" "<<elePt->at(eleInd)<<endl;
      
      TLorentzVector *ele = new TLorentzVector();
      ele->SetPtEtaPhiM(elePt->at(eleInd), eleEta->at(eleInd), elePhi->at(eleInd), 0.511/1000.);
      //ele->SetPtEtaPhiE(elePt->at(eleInd), eleEta->at(eleInd), elePhi->at(eleInd), eleEn->at(eleInd));
      
      for(int ipho=0; ipho<nGoodPho; ipho++){

	int phoInd = *(phoIndArr+ipho);
	
	if( std::find(usedPho.begin(), usedPho.end(), phoInd) != usedPho.end() ) //i.e. if found this pho in the list of paired 
	  {
	    if(debugL1){
	      cout<<"Found this photon in the list, photons index is  "<<phoInd<<endl;
	    }

	  continue;
	  }
	
	double phoSCEta = ecalSCeta->at(phoDirectEcalSCindex->at(phoInd));
	double phoSCPhi = ecalSCphi->at(phoDirectEcalSCindex->at(phoInd));

	TLorentzVector *pho = new TLorentzVector();
	pho->SetPtEtaPhiM(phoEt->at(phoInd), phoSCEta, phoSCPhi, 0);
	//pho->SetPtEtaPhiM(phoEt->at(phoInd), phoSCEta, phoSCPhi, 0.511/1000);
	//pho->SetPtEtaPhiE(phoEt->at(phoInd), phoSCEta, phoSCPhi, phoE->at(phoInd));
	

	TLorentzVector zP = *ele + * pho;
	double mass = zP.M();
	if(debug) cout<<"Mass "<<mass<<endl;

	double massDiff = fabs(91.18-mass);
	
	if( mass>60 && mass<120 && tmpmassDiff>massDiff )
	  {
	    foundPair = true;
	    tmpmassDiff = massDiff;
	    //elePhoInd[iele] = phoInd;

	    tmpelePhoInd = phoInd;


	    //usedPho[ipho] = 1;
	  }
	//cout<<"Mass "<<mass<<endl;

      }//for(int ipho=0; ipho<nGoodPho; ipho++)
      
      if(foundPair){
	if(debugL1) 
	  cout<<"Found pair .... paired ele "<<eleInd<<" with pho "<<tmpelePhoInd<<endl;


	usedPho.push_back(tmpelePhoInd);
	elePhoPair.push_back(make_pair (eleInd,tmpelePhoInd));
      }
    }// for(int iele=0; iele<nGoodEle; iele++)
    
    ///now all the pairs are formed - loop over pairs to put into num and deno
    int totPair = elePhoPair.size();
    selTotPairs = totPair;
    
    for(int ipair=0; ipair<totPair; ipair++){
      

      int eleInd = elePhoPair[ipair].first;
      int phoInd = elePhoPair[ipair].second;
      
      if(debugL1){
	std::cout<<"Pairs are ele pho "<<eleInd<<" "<<phoInd<<endl;
	cout<<" ele pT : pho pT : "<<elePt->at(eleInd)<<" "<<phoEt->at(phoInd)<<endl;

      }


      TLorentzVector *ele = new TLorentzVector();
      ele->SetPtEtaPhiM(elePt->at(eleInd), eleEta->at(eleInd), elePhi->at(eleInd), 0.511/1000.);
      //ele->SetPtEtaPhiE(elePt->at(eleInd), eleEta->at(eleInd), elePhi->at(eleInd), eleEn->at(eleInd));
      
      double eleSCEta = ecalSCeta->at(eleDirectEcalSCindex->at(eleInd));
      double eleSCPhi = ecalSCphi->at(eleDirectEcalSCindex->at(eleInd));

      double phoSCEta = ecalSCeta->at(phoDirectEcalSCindex->at(phoInd));
      double phoSCPhi = ecalSCphi->at(phoDirectEcalSCindex->at(phoInd));
      
      
      TLorentzVector *pho = new TLorentzVector();
      pho->SetPtEtaPhiM(phoEt->at(phoInd), phoSCEta, phoSCPhi, 0);
      //pho->SetPtEtaPhiM(phoEt->at(phoInd), phoSCEta, phoSCPhi, 0.511/1000);
      //pho->SetPtEtaPhiE(phoEt->at(phoInd), phoSCEta, phoSCPhi, phoE->at(phoInd)); 
      
      
      TLorentzVector zP = *ele + *pho;
      double mass = zP.M();

      if(debugL1){
	cout<<"Mass is "<<mass<<endl;
      }

      //cout<<"Mass is "<<mass<<endl;
      //cout<<"Ele: Px Py Pz E "<<ele->Px()<<" "<<ele->Py()<<" "<<ele->Pz()<<" "<<ele->Energy()<<endl;
      //cout<<"Pho: Px Py Pz E "<<pho->Px()<<" "<<pho->Py()<<" "<<pho->Pz()<<" "<<pho->Energy()<<endl;

      int phopixVeto = phoQualityBits->at(phoInd)>>0&1; 	 
      int phoeleVeto = phoQualityBits->at(phoInd)>>1&1; 

      selphoHasPix = phopixVeto; //has pixel seed
      
      selphoPt = pho->Pt();
      selphoEta = pho->Eta();
      selphoPhi = pho->Phi();
      selphoSCEta = phoSCEta;
      selphoSCPhi = phoSCPhi;
      
      selPairMass = mass;      

      selelePt = ele->Pt();
      seleleEta = ele->Eta();
      selelePhi = ele->Phi();


      seleleSCEta = eleSCEta;
      seleleSCPhi = eleSCPhi;

      ///check if this photon is actually reconstructed as an electron
      ////10th June, 2020 - check which electron is tracker driven and ECAL driven
      ///consider only those photons which are matched to ECAL tracker electrons
      int matchedIele = -99;
      bool foundEleMatch = isPhoAnEle(phoInd, matchedIele);
      
      selphoIsMatchEle = (foundEleMatch==true) ? 1 : 0;
      if(selphoIsMatchEle){
	
	selphoElePt = elePt->at(matchedIele);
	selphoEleEta = eleEta->at(matchedIele);
	selphoElePhi = elePhi->at(matchedIele);

	double eleSCEta = ecalSCeta->at(eleDirectEcalSCindex->at(matchedIele));
	double eleSCPhi = ecalSCphi->at(eleDirectEcalSCindex->at(matchedIele));
	
	//cout<<"MatchedElle Ind is "<<matchedIele<<" pt "<<elePt->at(matchedIele)<<" pt ffrom tree "<<selphoElePt<<endl;
	selphoEleSCEta = eleSCEta;
	selphoEleSCPhi = eleSCPhi;

      }
      

      tmintree->Fill();
      
    }

    
  }//for (Long64_t jentry=0; jentry<nentries;jentry++)
  
  fmintree->cd();
  tmintree->Write();
  fmintree->Write();
   
}


int* Analyse::selectEle(int &nGoodEle, int *eleIndArr){
  
  //int eleIndArr[100];
  
  int igoodele = 0;
  
  if(debug) cout<<"inside sel ele "<<endl;
  
  for(int iele=0; iele<nEle; iele++){
    
    int eleID = eleIDbit->at(iele)>>3&1;
    
    if(!eleID)
      continue;
    
    if(elePt->at(iele) < minelePt)
      continue;
    
    eleIndArr[igoodele] = iele;
    if(debug) cout<<"iele : pT "<<iele<<" "<<elePt->at(iele)<<endl;
    igoodele++;
    
  }

  nGoodEle = igoodele;
  return eleIndArr;
}


int* Analyse::selectPho(int &nGoodPho, int *phoIndArr, string idType){

  //int phoIndArr[100];
  int igoodpho = 0;
  
  double maxPhoPt = -999;
  for(int ipho=0; ipho<nPho; ipho++){
    
    if(phoEt->at(ipho) < minphoPtcut) continue;

    ///Fall v2 ID - Medium
    
    double cutHoE[2] = {0.02197, 0.0326};
    double cutSieie[2] = {0.01015, 0.0272};
    double cutNeIso[2][3] = { {1.189, 0.01512, 2.259e-5},  //const, linear, quadratic
			      {2.718, 0.0117, 2.3e-5}
    };


    double cutPhIso[2][2] = { {2.08, 0.004017},  //const, linear
			      {3.867, 0.0037}
    };
    
    double cutChIso[2] = {1.146, 1.056}; //cut in EE is made up following the EGM twiki and diff between EB cuts wrt POG 
    
    double cutSieieMin = 0.001;
    double cutSipipMin = 0.001;
    double cutEMip = 4.9;

    double cutSeedTime = 3;
    
    ///in addition reject phi regions in the EE

    int reg = 0;
    double phoSCEta = ecalSCeta->at(phoDirectEcalSCindex->at(ipho));
    double phoSCPhi = ecalSCphi->at(phoDirectEcalSCindex->at(ipho));

    if( fabs(phoSCEta)<1.4442 )
      reg = 0;
    
    else
      if( fabs(phoSCEta)>1.566 && fabs(phoSCEta)<2.5 ) 
	reg = 1;
    
      else
	continue;

    /*
    int isRealPho = !phopixVeto[ipho];
    
    if(selReal){
      if(!isRealPho)
	continue;
    }
    
    else{
      if(isRealPho)
	continue;

    }
    */
    /*
    if( phoHoverE[ipho] > cutHoE[reg] ) continue;
    if( phoSieie[ipho] > cutSieie[reg] ) continue;
    if( phoNhiso[ipho] > (cutNeIso[reg][0]+cutNeIso[reg][1]*phoPt[ipho]+cutNeIso[reg][2]*pow(phoPt[ipho],2)) ) continue;
    if( phoPiso[ipho] > (cutPhIso[reg][0]+cutPhIso[reg][1]*phoPt[ipho]) ) continue;
    if( phoChwiso[ipho] > cutChIso[reg] ) continue;
    
    if( phoSieie[ipho] < cutSieieMin ) continue;
    if( phoSipip[ipho] < cutSipipMin ) continue;
    
    //if( phoMIPTotEnergy[ipho] > cutEMip ) continue;
    if( fabs(phoSeedTime[ipho]) > cutSeedTime ) continue;
    */
    /*
    0 - 70% eff
1 - 80% eff
2 - 90% eff
3 - 95% eff
    */
    
    int phoID = pho_BDT_ID_bits->at(ipho)>>2&1; ///90

    if(idType=="70") phoID = pho_BDT_ID_bits->at(ipho)>>0&1; ///90
    if(idType=="80") phoID = pho_BDT_ID_bits->at(ipho)>>1&1; ///90
    if(idType=="90") phoID = pho_BDT_ID_bits->at(ipho)>>2&1; ///90
    if(idType=="95") phoID = pho_BDT_ID_bits->at(ipho)>>3&1; ///90

    if(!phoID) continue;
    
    phoIndArr[igoodpho] = ipho;
    igoodpho++;
    
  }//for(int ipho=0; ipho<nPho; ipho++)

  nGoodPho = igoodpho;

  return phoIndArr;
  
}

/////////////////////////////////////////////////////////////
bool Analyse::isPhoAnEle(int ipho, int &matchedIele){
  
  double pi = 4*atan(1.);
  
  bool foundMatch = false;
  matchedIele = -99;

  for(int iele=0; iele<nEle; iele++){
    
    double eleSCEta = ecalSCeta->at(eleDirectEcalSCindex->at(iele));
    double eleSCPhi = ecalSCphi->at(eleDirectEcalSCindex->at(iele));
    
    double phoSCEta = ecalSCeta->at(phoDirectEcalSCindex->at(ipho));
    double phoSCPhi = ecalSCphi->at(phoDirectEcalSCindex->at(ipho));
    
    double dEta = fabs(eleSCEta-phoSCEta);
    double dPhi = fabs(eleSCPhi-phoSCPhi);
    if(dPhi>pi) dPhi = 2*pi - dPhi;
    
    double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
    
    //if(dEta < 0.02 && dPhi < 0.02){
    if(dR<0.05){

      foundMatch = true;
      matchedIele = iele;
      break;
    }//if(dEta < 0.05 && dPhi < 0.05)
  }//for(int iele=0; iele<nEle; iele++)
  
  return foundMatch;
}

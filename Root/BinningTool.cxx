/****************************************************
*****************************************************
Binning Tool
Author: Brock Moir

This class handles the creation and manipulation of the analysis bins.  

*****************************************************
*****************************************************/

#include <Smearing/BinningTool.h>

#include <TFile.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
// this is needed to distribute the algorithm to the workers
ClassImp(BinningTool);

BinningTool :: BinningTool ()
{

}


BinningTool ::~BinningTool()
{
}


void BinningTool :: SetBins (std::vector <double> eta, 
			     std::vector <double> phi, 
			     std::vector <double> pt)
{
  
  etabin = eta;
  phibin = phi;
  ptbin = pt;
  nbins = 0;
}

void BinningTool :: SetBins ()
{

  //double eta[]={-5., -2.6,-1.,0.0,1.,2.6, 5.};
  double eta[]={-4.9,-3.2,-1.7,-1,.0,1,1.7,3.2,4.9};

  double phi[]={-TMath::Pi(), -TMath::Pi()/2.0, 0,TMath::Pi()/2.0, TMath::Pi()};
  //double phi[]={-TMath::Pi(), TMath::Pi()};

  double pt[]={0.0,45000.0,70000.0, 100000000000.0};
  //double pt[]={30000.0,50000.0,150000.0, 1000000000000,0};

  for(unsigned int i=0; i < sizeof (eta) / sizeof (double); i++){
    etabin.push_back(eta[i]);
  }

  for(unsigned int i=0; i < sizeof (phi) / sizeof (double); i++){
    phibin.push_back(phi[i]);
  }

  for(unsigned int i=0; i < sizeof (pt) / sizeof (double); i++){
    ptbin.push_back(pt[i]);
  }

}

void BinningTool :: Initialize()
{

  for(unsigned int i=1; i < etabin.size(); i++){
    std::string Eta = "Eta";
    std::ostringstream iEta;
    iEta << i-1;
    Eta+=iEta.str();
    if(i==1 || i == etabin.size()-1) {
      binnames.push_back(Eta+"Phi0Pt0");
      nbins++;
    }else{
      for(unsigned int j=1; j < phibin.size(); j++){
	std::string Phi = "Phi";
	std::ostringstream iPhi;
        iPhi << j-1;
        Phi+=iPhi.str();
        for(unsigned int k=1; k < ptbin.size(); k++){
	  std::string Pt = "Pt";
	  std::ostringstream iPt;
          iPt << k-1;
          Pt+=iPt.str(); 
          binnames.push_back(Eta+Phi+Pt);
	  nbins++;
        }
      }
    }
  } 

  for(unsigned int i=1; i < etabin.size(); i++){
    std::string Eta = "Eta";
    std::ostringstream iEta;
    iEta << i-1;
    Eta+=iEta.str(); 
    etanames.push_back(Eta);
  }
  for(unsigned int i=1; i < phibin.size(); i++){
    std::string Phi = "Phi";
    std::ostringstream iPhi;
    iPhi << i-1;
    Phi+=iPhi.str(); 
    phinames.push_back(Phi);
  }
  for(unsigned int i=1; i < ptbin.size(); i++){
    std::string Pt = "Pt";
    std::ostringstream iPt;
    iPt << i-1;
    Pt+=iPt.str(); 
    ptnames.push_back(Pt);

  }

}

std::string BinningTool :: GetBinName (double eta, double phi, 
				       double pt)//, int run)
{

  std::vector <int> index = this->GetBinIndices(eta, phi, pt);//, run);
  return binnames[index[0]];
}


std::vector <int> BinningTool :: GetBinIndices (double eta, double phi, 
				double pt)
{
  int i=0, j=0, k=0;
  for(unsigned int ieta=1; ieta < etabin.size(); ieta++){
    if(eta <= etabin[ieta]){
      i = ieta-1;
      break;
    }
  }  

  for(unsigned int iphi=1; iphi < phibin.size(); iphi++){
    if(phi <= phibin[iphi]){
      j = iphi-1;
      break;
    }
  }

  for(unsigned int ipt=1; ipt < ptbin.size(); ipt++){ 
    if(pt <= ptbin[ipt]){
      k = ipt-1;
      break;
    }
  }

  if(i==0 || i==etabin.size()-2) j=0, k=0;

  std::vector <int> indices;
  std::string BinName;
  std::string Eta = "Eta";
  std::string Phi = "Phi";
  std::string Pt = "Pt";
  std::ostringstream iEta, iPhi, iPt;
  iEta << i; iPhi << j; iPt << k;
  Eta+=iEta.str(); Phi+=iPhi.str(); Pt+=iPt.str(); 
  BinName=Eta+Phi+Pt;

  for(unsigned int l=0; l <binnames.size(); l++){
    if(binnames[l] == BinName){
      indices.push_back(l);
      indices.push_back(i);
      indices.push_back(j);
      indices.push_back(k);
      return indices;
    }
  }

  return indices;
}

int BinningTool::GetBinIndex(int ieta, int iphi, int ipt){
  std::string histname = etanames[ieta]+phinames[iphi]+ptnames[ipt];
  for(int i=0; i<binnames.size(); i++){
    if (histname == binnames[i]) return i;
  }

}

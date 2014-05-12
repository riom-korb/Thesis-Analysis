/****************************************************
*****************************************************
Unfolding_Tool
Author: Brock Moir

This class produces the smearing histogram for each bin.  

*****************************************************
*****************************************************/

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include <ROOUnfold/RooUnfold.h>
#include <ROOUnfold/RooUnfoldResponse.h>
#include <ROOUnfold/RooUnfoldBayes.h>

#include <Smearing/UnfoldingTool.h>

#include <TH2.h>
#include <TH1.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

// this is needed to distribute the algorithm to the workers
ClassImp(UnfoldingTool);


UnfoldingTool :: UnfoldingTool ()
{

};


float UnfoldingTool::Chi2(TH1F* h1, TH1F* h2){
  Int_t n_bins  = h1->GetNbinsX();
  float chi2 = 0;
  float Pvalue = 0;
  float ndf = -1.;
  float x1, x2;

  for(int i=0; i<n_bins; i++){
    x1 = h1->GetBinContent(i+1);
    x2 = h2->GetBinContent(i+1);
    if(x1+x2!=0.0) {
      chi2 += (((x1-x2)*(x1-x2))/(x1+x2));   
      ndf += 1.;
    }
  }
  
  Pvalue =  this->igf(ndf*0.5, chi2 * 0.5);
  Pvalue /= tgamma(ndf*0.5);
  return 1.-Pvalue;

}  

double UnfoldingTool::igf(double S, double Z)
{
    if(Z < 0.0)
    {
	return 0.0;
    }
    double Sc = (1.0 / S);
    Sc *= pow(Z, S);
    Sc *= exp(-Z);
 
    double Sum = 1.0;
    double Nom = 1.0;
    double Denom = 1.0;
 
    for(int I = 0; I < 200; I++)
    {
	Nom *= Z;
	S++;
	Denom *= S;
	Sum += (Nom / Denom);
    }
 
    return Sum * Sc;
}

double UnfoldingTool::KSTest(TH1F* h1, TH1F* h2)
{
  double ks = 0., diff = 0.;
  double sum1 = 0, sum2 = 0;

  for (int i=0; i<h1->GetNbinsX(); i++){
    sum1 += h1->GetBinContent(i+1);
    sum2 += h2->GetBinContent(i+1);
    diff = fabs(sum1-sum2);
    if (diff > ks) ks = diff;
  }

  double sum = 0;
  for (int j=1; j<400; j++){
    sum += 2.*pow(-1., j-1.)*exp(-2. * (double) j * ks * (double) j * ks);
  }

  return sum;
}

/*


TH1F *UnfoldingTool :: BayesIter(TH1F* data, TH1F* mc, TH1F* guess){

  //prepare the histograms
  Int_t n_bins  = data->GetNbinsX();
  Int_t nzero   = mc->FindBin(0); 
  Double_t xmin = data->GetXaxis()->GetXmin();
  Double_t xmax = data->GetXaxis()->GetXmax();
  TH1F *res = new TH1F("res", "res", n_bins, xmin, xmax);

  //variables for unfolding
  float n_not[n_bins], P_not[n_bins], 
        n_hat[n_bins], P_hat[n_bins],
        eff[n_bins], 
        theta[n_bins][n_bins], lambda[n_bins][n_bins],
        denominator[n_bins], 
        chi2=100, sum1=0, sum2=0;

  //build up the smearing matrix
  for(int i=0; i<n_bins; i++){
    for(int j=0; j<n_bins; j++){
      if(i-j+nzero<n_bins && i-j+n_bins>0) lambda[i][j] = 
	                                   mc->GetBinContent(i-j+nzero);
      else lambda[i][j]=0;
    }
  }  

  //ititialize n_not & P_not 
  for(int i=0; i<n_bins; i++){
    n_not[i] = guess->GetBinContent(i+1);
    sum1 += n_not[i];
  }
  for(int i=0; i<n_bins; i++){
    if(sum1 == 0) P_not[i]=1.0/100.0;
    else P_not[i] = n_not[i]/sum1;
  }

  while(chi2>0.5){
    //calculate theta's denominator, efficiency
    for(int i=0; i<n_bins; i++){
      sum1=0, sum2=0;
      for(int j=0; j<n_bins; j++){
        sum1 += lambda[j][i] * P_not[j];
        sum2 += lambda[j][i];
      }
      denominator[i] = sum1;
      eff[i]         = sum2;
    }

    //calculate theta
    for(int i=0; i<n_bins; i++){
      for(int j=0; j<n_bins; j++){
        if(denominator[j] == 0.0) theta[i][j] = 0;
        else theta[i][j] = (lambda[j][i] * P_not[i]) / denominator[j];
      }
    }

    //calculate n_hat, P_hat
    sum2=0;
    for(int i=0; i<n_bins; i++){
      if(eff[i] == 0) n_hat[i] = 0;
      else{
        sum1=0;
        for(int j=0; j<n_bins; j++){
          sum1 += theta[i][j] * data->GetBinContent(j+1);
        }
        n_hat[i] = sum1 / eff[i];
        sum2 += n_hat[i];
      }    }
    for(int i=0; i<n_bins; i++) P_hat[i] = n_hat[i] / sum2;

    //calculate chi2, and carry over the iteration
    chi2 = this->Chi2(P_hat, P_not, n_bins);
    for(int i=0; i<n_bins; i++){
      n_not[i] = n_hat[i];
      P_not[i] = P_hat[i];
    }
    std::cout<<chi2<<std::endl;
  }

  for(int i=0; i<n_bins; i++){
    res->SetBinContent(i+1, P_hat[i]);
  }

  return res;
}
*/
EL::StatusCode UnfoldingTool::MakeSmearingHistBayes(TH1F * hData, TH1F * hMC, TH1F * hSmear, int niter, int rebinNum)
 {
  hSmear->Reset();
  hData->Rebin(rebinNum);
  hMC->Rebin(rebinNum);
  hSmear->Rebin(rebinNum);

  Int_t n = hMC->GetNbinsX();
  Double_t xmin = hMC->GetXaxis()->GetXmin();
  Double_t xmax = hMC->GetXaxis()->GetXmax();
  std::ostringstream t_name;
  t_name << hMC->GetName() << "mat";
  TString name = t_name.str();
  TH2F *mat = new TH2F(name, name, n, xmin, xmax, n, xmin, xmax);
  int nzero=hMC->FindBin(0);
  for(int i=0; i<n; i++){
     for(int j=0; j<n; j++){
       if(i-j+nzero<n){
         mat->SetBinContent(i+1, j+1, hMC->GetBinContent(i-j+nzero));
         mat->SetBinError(i+1, j+1, hMC->GetBinError(i-j+nzero));
       }
     }
  }

  Root::RooUnfoldResponse *R = new Root::RooUnfoldResponse(0, 0, mat);
  Root::RooUnfoldBayes *unfold = new Root::RooUnfoldBayes(R, hData, niter);
  // hSmear = (TH1F*) unfold->Hreco();
   
  for(int i=0; i<n; i++){
    hSmear->SetBinContent(i, unfold->Hreco()->GetBinContent(i));
    hSmear->SetBinError(i, unfold->Hreco()->GetBinError(i));
  }

   return EL::StatusCode::SUCCESS;    
 }

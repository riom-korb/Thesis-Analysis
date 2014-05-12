/****************************************************
*****************************************************
Make_Bins
Author: Brock Moir

This class iterates over the direct photon events calculating 
Z and matching it to the correct bin.  

*****************************************************
*****************************************************/

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include <Smearing/Make_Bins.h>
#include <Smearing/BinningTool.h>

#include <TFile.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>


#define XXX std::cout << " I am here: " << __FILE__ << ":" << __LINE__ << std::endl;

// this is needed to distribute the algorithm to the workers
ClassImp(Make_Bins);

Make_Bins :: Make_Bins ()
{
 
  m_dosmearing = false;
  m_smearingfile = 0;

  m_analysistree = 0;
  m_outputStreamName = "output";
  m_bintool = 0;
 
};

EL::StatusCode Make_Bins :: setupJob (EL::Job& job)
{
 
  EL::OutputStream stream (m_outputStreamName.Data());
  job.outputAdd	(stream);
  job.useD3PDReader ();
 
  return EL::StatusCode::SUCCESS;
};

EL::StatusCode Make_Bins :: changeInput (bool firstFile)
{
 
  m_analysistree= new AnalysisTree(wk()->tree());
 
  return EL::StatusCode::SUCCESS;
};

EL::StatusCode Make_Bins :: initialize ()
{
 
  TH1::SetDefaultSumw2();  
  TH2::SetDefaultSumw2();

  outputFile= wk()->getOutputFile (m_outputStreamName.Data());

  if(m_dosmearing) m_smearingfile = new TFile(m_smearpath.c_str());

  m_bintool->Initialize();

  this->BookHistos();
  this->BookExtras();
 
  return EL::StatusCode::SUCCESS;
};



EL::StatusCode Make_Bins :: execute ()
{
  // This method is called for each iteration

  m_analysistree->GetEntry(wk()->treeEntry());

  double weight    = m_analysistree->weight,
         jet_pt    = (*m_analysistree->jet_pt)[0],
         jet_eta   = (*m_analysistree->jet_eta)[0], 
         jet_phi   = (*m_analysistree->jet_phi)[0],
         gamma_pt  = (*m_analysistree->ph_pt)[0],
         gamma_eta = (*m_analysistree->ph_eta)[0], 
         gamma_phi = (*m_analysistree->ph_phi)[0],
         MET_pt    = m_analysistree->MET_RefFinal_et,
         MET_phi   = m_analysistree->MET_RefFinal_phi,
         Res       = m_analysistree->Res,
         thresh    = m_analysistree->TrigThreshold,
         prescale  = m_analysistree->TrigPrescale,
         scale20   = 12.7722,
         scale40   = 176.526,
         scale60   = 999.755,
         scale80   = 4587.2,
         scale150  = 4587.2;

  std::vector <int> index = m_bintool->GetBinIndices(jet_eta, jet_phi, jet_pt);
  std::string binname = m_bintool->GetBinName(jet_eta, jet_phi, jet_pt);

  if( jet_pt/1000. < 30) return EL::StatusCode::SUCCESS;  
  if(thresh==150) return EL::StatusCode::SUCCESS;

  if(! m_analysistree->isSimulation){
    if(thresh==20){
      weight = 1./scale20;
      h_ptph_trig20->Fill(gamma_pt/1000., weight);
      h_ptph_trig20_up->Fill(gamma_pt/1000., 1./scale20);
    }if(thresh==40){
      weight = 1./scale40;
      h_ptph_trig40->Fill(gamma_pt/1000., weight);
      h_ptph_trig40_up->Fill(gamma_pt/1000., 1./scale40);
    }if(thresh==60){
      weight = 1./scale60;
      h_ptph_trig60->Fill(gamma_pt/1000., weight);
      h_ptph_trig60_up->Fill(gamma_pt/1000., 1./scale60);
    }if(thresh==80){
      weight = 1./scale80;
      h_ptph_trig80->Fill(gamma_pt/1000., weight);
      h_ptph_trig80_up->Fill(gamma_pt/1000., 1./scale80);
    }
  }

  hist_allbins[index[0]]->Fill(Res, weight);
  hist_allres->Fill(Res, weight);

  h_ptph_bins[index[0]]->Fill(gamma_pt/1000., weight);

  hist_etabins[index[1]]->Fill(Res, weight);
  if(index[1] != 0 and index[1] != m_bintool->etanames.size()-1){
    hist_phibins[index[2]]->Fill(Res, weight);
    hist_ptbins[index[3]]->Fill(Res, weight);
  }  

  hphijetmet->Fill(jet_phi, MET_phi, weight);
  hptjetph  ->Fill(jet_pt/1000., gamma_pt/1000., weight); 

  h_jetpt   ->Fill(jet_pt/1000., weight); 
  h_phpt    ->Fill(gamma_pt/1000., weight); 
  h_phpt_up ->Fill(gamma_pt/1000., weight/prescale); 
  h_diffpt  ->Fill((jet_pt- gamma_pt)/1000., weight); 

  h_jeteta  ->Fill(jet_eta, weight); 
  h_pheta   ->Fill(gamma_eta, weight);

  h_jetphphi->Fill(jet_phi-gamma_phi, weight);
  h_jetphphi->Fill(gamma_phi-jet_phi, weight);

  if(index[0] == 0){
    h_lumi_freq->Fill(thresh);
    h_unweighted->Fill(Res);
  }

  hphijetph->Fill(jet_phi, gamma_phi, weight);
  
  return EL::StatusCode::SUCCESS;
};



EL::StatusCode Make_Bins :: finalize ()
{
  // do everything you need to do before the job on the
  // worker node ends.  so far I have no example of what
  // that would be, but the hook is here just in case.
  return EL::StatusCode::SUCCESS;
};


void Make_Bins :: BookHistos()
{

  int nbins=100;
  double xmin=-2, xmax=2;

  for(unsigned int i=0; i < m_bintool->binnames.size(); i++){
    TH1F *htemp = new TH1F((m_bintool->binnames[i]).c_str(), (m_bintool->binnames[i]).c_str(), nbins, xmin, xmax);
    hist_allbins.push_back(htemp);
    wk()->addOutput(htemp);

    TString phptname = "ptph_" + (m_bintool->binnames[i]);
    TH1F *htemp1 = new TH1F(phptname, phptname, nbins, 0, 150);
    h_ptph_bins.push_back(htemp1);
    wk()->addOutput(htemp1);
  }

  for(unsigned int i=0; i < m_bintool->etanames.size(); i++){
    TH1F *htemp = new TH1F((m_bintool->etanames[i]).c_str(), (m_bintool->etanames[i]).c_str(), nbins, xmin, xmax);
    hist_etabins.push_back(htemp);
    wk()->addOutput(htemp);
  }

  for(unsigned int i=0; i < m_bintool->phinames.size(); i++){
    TH1F *htemp = new TH1F((m_bintool->phinames[i]).c_str(), (m_bintool->phinames[i]).c_str(), nbins, xmin, xmax);
    hist_phibins.push_back(htemp);
    wk()->addOutput(htemp);
  }

  for(unsigned int i=0; i < m_bintool->ptnames.size(); i++){
    TH1F *htemp = new TH1F((m_bintool->ptnames[i]).c_str(), (m_bintool->ptnames[i]).c_str(), nbins, xmin, xmax);
    hist_ptbins.push_back(htemp);
    wk()->addOutput(htemp);
  }

  hist_allres    = new TH1F("allres", "allres", 100, -2, 2);
  wk()->addOutput(hist_allres);

  h_MET = new TH1F("MET", "MET", 100, 0, 100000);
  wk()->addOutput(h_MET);

  h_sm_MET = new TH1F("sm_MET", "sm_MET", 100, 0, 100000);
  wk()->addOutput(h_sm_MET);

}

void Make_Bins :: BookExtras()
{
  float ptmax = 150;
  
  hphijetmet = new TH2F("phijetmet", "phijetmet", 100, -3.14, 3.14, 100, -3.14, 3.14); 
  wk()->addOutput(hphijetmet);

  hptjetph   = new TH2F("ptjetph", "ptjetph", 100, 30, ptmax, 100, 30, ptmax);
  wk()->addOutput(hptjetph);

  hphijetph   = new TH2F("phijetph", "phijetph", 100, -3.14, 3.14, 100, -3.14, 3.14);
  wk()->addOutput(hphijetph);

  h_jetpt = new TH1F("jetpt", "jetpt", 100, 0, ptmax);
  wk()->addOutput(h_jetpt);

  h_phpt = new TH1F("phpt", "phpt", 100, 0, ptmax);
  wk()->addOutput(h_phpt);

  h_phpt_up = new TH1F("phpt_up", "phpt_up", 100, 0, ptmax);
  wk()->addOutput(h_phpt_up);

  h_diffpt = new TH1F("diffpt", "diffpt", 100, -ptmax, ptmax);
  wk()->addOutput(h_diffpt);

  h_jeteta = new TH1F("jeteta", "jeteta", 100, -4.9, 4.9);
  wk()->addOutput(h_jeteta);

  h_pheta = new TH1F("pheta", "pheta", 100, -4.9, 4.9);
  wk()->addOutput(h_pheta);

  h_jetphphi = new TH1F("jetphphi", "jetphphi", 100, 3., 3.14);
  wk()->addOutput(h_jetphphi);

  h_ptph_trig20 = new TH1F("ptph_trig20", "ptph_trig20", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig20);
  h_ptph_trig40 = new TH1F("ptph_trig40", "ptph_trig40", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig40);
  h_ptph_trig60 = new TH1F("ptph_trig60", "ptph_trig60", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig60);
  h_ptph_trig80 = new TH1F("ptph_trig80", "ptph_trig80", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig80);
  h_ptph_trig150 = new TH1F("ptph_trig150", "ptph_trig150", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig150);

  h_ptph_trig20_up = new TH1F("ptph_trig20_up", "ptph_trig20_up", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig20_up);
  h_ptph_trig40_up = new TH1F("ptph_trig40_up", "ptph_trig40_up", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig40_up);
  h_ptph_trig60_up = new TH1F("ptph_trig60_up", "ptph_trig60_up", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig60_up);
  h_ptph_trig80_up = new TH1F("ptph_trig80_up", "ptph_trig80_up", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig80_up);
  h_ptph_trig150_up = new TH1F("ptph_trig150_up", "ptph_trig150_up", 100, 0., ptmax);
  wk()->addOutput(h_ptph_trig150_up);

  h_lumi_freq = new TH1F("lumi_freq", "lumi_freq", 100, 0, 160);
  wk()->addOutput(h_lumi_freq);
  h_unweighted = new TH1F("unweighted", "unweighted", 100, -2, 2);
  wk()->addOutput(h_unweighted);
}

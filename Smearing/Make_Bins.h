#ifndef Smearing_Make_Bins_H
#define Smearing_Make_Bins_H

#include <EventLoop/Algorithm.h>
#include <Smearing/AnalysisTree.h>
#include <Smearing/BinningTool.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>

class Make_Bins : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  std::vector <TH1F*>  hist_allbins;   //!
  std::vector <TH1F*>  hist_etabins; //!
  std::vector <TH1F*>  hist_phibins; //!
  std::vector <TH1F*>  hist_ptbins;  //!

  std::vector <TH1F*>  h_ptph_bins; //!

  TH1F*  h_ptph_trig20; //!
  TH1F*  h_ptph_trig40; //!
  TH1F*  h_ptph_trig60; //!
  TH1F*  h_ptph_trig80; //!
  TH1F*  h_ptph_trig150; //!

  TH1F*  h_ptph_trig20_up; //!
  TH1F*  h_ptph_trig40_up; //!
  TH1F*  h_ptph_trig60_up; //!
  TH1F*  h_ptph_trig80_up; //!
  TH1F*  h_ptph_trig150_up; //!

  TH1F* h_lumi_freq; //!
  TH1F* h_unweighted; //!

  TH1F* hist_allres; //!
  TH1F* hsmear; //!
  TH1F * h_MET; //!
  TH1F * h_sm_MET; //!
  TH1F * h_jetpt; //!
  TH1F * h_phpt; //!
  TH1F * h_phpt_up; //! 
  TH1F * h_diffpt; //!

  TH1F * h_jetphphi; //!

  TH1F * h_jeteta; //!
  TH1F * h_pheta; //!

  TH2F* hphijetmet; //!
  TH2F* hptjetph; //!
  TH2F* hphijetph; //!

  TFile *outputFile; //!
  TFile *m_smearingfile; 
  AnalysisTree *m_analysistree; //!
  BinningTool *m_bintool;
  TString m_outputStreamName; 

  std::string m_smearpath;

  bool m_dosmearing; 
  // this is a standard constructor
  Make_Bins ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode finalize ();

  void BookHistos();
  void BookExtras();

  // this is needed to distribute the algorithm to the workers
  ClassDef(Make_Bins, 1);
};

#endif

#ifndef Smearing_SmearJets_H
#define Smearing_SmearJets_H

#include <EventLoop/Algorithm.h>

#include <Smearing/AnalysisTree.h>
#include <Smearing/BinningTool.h>

#include <TH1.h>
#include <TH2.h>

class SmearJets : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  std::vector <TH2F*>  histos2d;
  std::vector <std::vector <TH1F*> >  histos;

  std::vector <TH2F*>  hist_etabins;
  std::vector <std::vector <TH1F*> >  h_etabins;

  std::vector <TH2F*>  hist_phibins;
  std::vector <std::vector <TH1F*> >  h_phibins;

  std::vector <TH2F*>  hist_ptbins;
  std::vector <std::vector <TH1F*> >  h_ptbins;

  std::vector <TH1F*> allres;
  std::vector <TH1F*> allresbin;
  std::vector <TH1F*> histav;

  std::vector <TH1F*> hchi2;
  std::vector <TH1F*> hNDF;
  std::vector <TH1F*> hchi2NDF;

  std::vector <TH1F*> hchi2_pt;
  std::vector <TH1F*> hNDF_pt;
  std::vector <TH1F*> hchi2NDF_pt;

  std::vector <TH1F*> hchi2_phi;
  std::vector <TH1F*> hNDF_phi;
  std::vector <TH1F*> hchi2NDF_phi;

  std::vector <TH1F*> hchi2_eta;
  std::vector <TH1F*> hNDF_eta;
  std::vector <TH1F*> hchi2NDF_eta;

  std::vector <TH1F*> htjetpt;
  std::vector <TH1F*> htdiffpt;

  TH2F* hjetpt;
  TH2F* hdiffpt;

  TFile *outputFile;
  TFile *m_smearingfile;
  AnalysisTree *m_analysistree;
  BinningTool *m_bintool;
  TString m_outputStreamName;
  int m_niter;
  std::string m_smearpath;

  TH1F *hsmear; //!
  TH1F *hallsmear; //!
  TH2F *hallres; //!
  TH2F *hallresbin; //!

  TH1F* hallchi2; //!
  TH1F* hallchi2NDF; //!
  TH1F* hallNDF; //!

  TH1F* hallchi2bin; //!
  TH1F* hallchi2NDFbin; //!
  TH1F* hallNDFbin; //!

  // this is a standard constructor
  SmearJets ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode finalize ();

  void BookHistos();

  // this is needed to distribute the algorithm to the workers
  ClassDef(SmearJets, 1);
};

#endif

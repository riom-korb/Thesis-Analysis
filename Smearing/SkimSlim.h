#ifndef JetResponse_SkimSlim_H
#define JetResponse_SkimSlim_H

#include <D3PDReader/Event.h>
#include <EventLoop/Algorithm.h>
#include <GoodRunsLists/TGoodRunsListReader.h>
#include <TrigRootAnalysis/TrigDecisionToolD3PD.h>
#include <JetSelectorTools/TJetCleaningTool.h>
#include <egammaAnalysisUtils/EnergyRescaler.h>
#include <ApplyJetCalibration/ApplyJetCalibration.h>
#include <PileupReweighting/TPileupReweighting.h>

#include <Smearing/BinningTool.h>

#include <vector>

class SkimSlim : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
private:

  D3PD::TrigDecisionToolD3PD *m_tdt; //!

public:

  bool m_directph;

 // event related variables
 D3PDReader::Event *event; //!
 int trigPre;
 float trigThresh;

 // object related variables
 TString m_jetName;
 TString m_METName;
 D3PDReader::ElectronD3PDObject   *final_electrons; //!
 D3PDReader::PhotonD3PDObject     *final_photons; //!
 D3PDReader::JetD3PDObject        *final_jets; //!
 D3PDReader::MuonD3PDObject       *final_muons; //!
 D3PDReader::TauD3PDObject        *final_taus; //!
 float elptcut;
 float phptcut;
 float jetptcut;
 float muptcut;
 float tauptcut;

 // output tree variables
 TFile *outputFile; //!
 TTree *outputTree; //!
 D3PDReader::PhotonD3PDObject     *out_photons; //!
 D3PDReader::JetD3PDObject        *out_jets; //!
 float res;
 float MClumi;
 float MCpileup;
 float MCeventweight;
 float event_weight;
 std::vector<TString> cutflow_event_names; 
 std::vector<TString> cutflow_directph_names;
 TH1F *cutflow_event; //! 
 TH1F *cutflow_directph; //! 


 // tools
 Root::TGoodRunsList m_grl;
 Root::TJetCleaningTool *m_JetCleaningTool;
 Root::TPileupReweighting *m_pileupTool;
 eg2011::EnergyRescaler* energyRescaler; //!  
 JetCalibrationTool *m_JES;

 // configure the algorithm 
 TString m_outputStreamName;
 bool m_isData;
 bool m_doTrigger;
 float m_deltaPhiB2B;
 std::vector<std::string> m_investigateTrig;
 std::vector<double> m_trigThreshold;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  SkimSlim ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode finalize ();

  void getElectrons();
  void getPhotons();
  void getJets();
  void getMuons();
  void getTaus();
  float deltaR(float, float, float, float);
  void BookCutflow();
  void FillCutflow_event(TString);
  void FillCutflow_directph(TString);

  // this is needed to distribute the algorithm to the workers
  ClassDef(SkimSlim, 1);
};

#endif

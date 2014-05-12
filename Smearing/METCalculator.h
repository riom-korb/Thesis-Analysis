#ifndef Smearing_METCalculator_H
#define Smearing_METCalculator_H

#include <EventLoop/Algorithm.h>

#include <Smearing/AnalysisTree.h>
#include <Smearing/BinningTool.h>

#include <TH1.h>
#include <TH2.h>

class METCalculator : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  TFile *outputFile;
  TFile *m_smearingfile;
  AnalysisTree *m_analysistree;
  BinningTool *m_bintool;
  TString m_outputStreamName;
  std::string m_smearpath;

  TH1F *hsmear; //!
  TH1F * h_MET; //!
  TH1F * h_sm_MET; //!

  // this is a standard constructor
  METCalculator ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode finalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(METCalculator, 1);
};

#endif

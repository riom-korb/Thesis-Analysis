#ifndef Smearing_UnfoldingTool_H
#define Smearing_UnfoldingTool_H

#include <EventLoop/Algorithm.h>
#include "EventLoop/StatusCode.h"

#include <TH1.h>
#include <TF1.h>

class UnfoldingTool : public EL::Algorithm
{

public:

  UnfoldingTool ();

  float Chi2(TH1F*, TH1F*);
  double igf(double, double);
  double KSTest(TH1F*, TH1F*);

  TH1F* BayesIter(TH1F*, TH1F*, TH1F*);

  EL::StatusCode MakeSmearingHistBayes(TH1F*, TH1F*, TH1F*, int, int);

  // this is needed to distribute the algorithm to the workers
  ClassDef(UnfoldingTool, 1);
};

#endif

#ifndef Smearing_BinningTool_H
#define Smearing_BinningTool_H

#include <TNamed.h>
#include <TFile.h>
#include <TH1.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

class BinningTool : public TNamed{

  public:

    // constructor, destructor
    BinningTool ();

    ~BinningTool();

    // Bins should be set for every instance of BinningTool
    void SetBins(std::vector <double>, std::vector <double>, 
		 std::vector <double>);

    void SetBins();

    void Initialize();

    // input eta phi pt and run, get a bin name
    std::string GetBinName(double, double, double);

    // input eta phi pt and run, get the bin indices
    std::vector <int> GetBinIndices(double, double, double);

    int GetBinIndex(int, int, int);

    //private:

    // variables needed by this tool
    std::vector <double>      etabin;
    std::vector <double>      phibin;
    std::vector <double>      ptbin;

    std::vector <std::string> binnames;
    std::vector <std::string> etanames;
    std::vector <std::string> phinames;
    std::vector <std::string> ptnames;

    int nbins;

  // this is needed to distribute the algorithm to the workers
  ClassDef(BinningTool, 1);
};

#endif

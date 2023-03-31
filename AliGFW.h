/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
Class steers the initialization and calculation of n-particle correlations. Uses recursive function, all terms are calculated only once.
Latest version includes the calculation of any number of gaps and any combination of harmonics (including eg symmetric cumulants, etc.)
If used, modified, or distributed, please aknowledge the author of this code.
*/
#ifndef AliGFW__H
#define AliGFW__H
#include "AliGFWCumulant.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <complex>
using std::vector;
using std::complex;
using std::string;
class AliGFW {
 public:
  struct Region {
    int Nhar, Npar, NpT;
    vector<int> NparVec;
    double EtaMin=-999;
    double EtaMax=-999;
    int BitMask=1;
    string rName="";
    bool operator<(const Region& a) const {
      return EtaMin < a.EtaMin;
    };
    void PrintStructure() {printf("%s: eta [%f.. %f].",rName.c_str(),EtaMin,EtaMax); };
  };
  struct CorrConfig {
    vector<vector<int>> Regs {};
    vector<vector<int>> Hars {};
    vector<int> Overlap;
    vector<int> ptInd;
    bool pTDif=false;
    string Head="";
  };
  AliGFW();
  ~AliGFW();
  vector<Region> fRegions;
  vector<AliGFWCumulant> fCumulants;
  void AddRegion(string refName, int lNhar, int lNpar, double lEtaMin, double lEtaMax, int lNpT=1, int BitMask=1);
  void AddRegion(string refName, int lNhar, int *lNparVec, double lEtaMin, double lEtaMax, int lNpT=1, int BitMask=1);
  int CreateRegions();
  void Fill(double eta, int ptin, double phi, double weight, int mask, double secondWeight=-1);
  void Clear();
  AliGFWCumulant GetCumulant(int index) { return fCumulants.at(index); };
  CorrConfig GetCorrelatorConfig(string config, string head = "", bool ptdif=false);
  complex<double> Calculate(CorrConfig corconf, int ptbin, bool SetHarmsToZero, bool DisableOverlap=false);
public:
  bool fInitialized;
  complex<double> TwoRec(int n1, int n2, int p1, int p2, int ptbin, AliGFWCumulant*, AliGFWCumulant*, AliGFWCumulant*);
  complex<double> RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, int ptbin, vector<int> &hars, vector<int> &pows); //POI, Ref. flow, overlapping region
  complex<double> RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, int ptbin, vector<int> &hars); //POI, Ref. flow, overlapping region
  void AddRegion(Region inreg) { fRegions.push_back(inreg); };
  Region GetRegion(int index) { return fRegions.at(index); };
  int FindRegionByName(string refName);
  //Calculating functions:
  complex<double> Calculate(int poi, int ref, vector<int> hars, int ptbin=0); //For differential, need POI and reference
  complex<double> Calculate(int poi, vector<int> hars); //For integrated case
  //Operations on strings. Equivalent to TString operations, but one to rid of root dependence
  int s_index(string &instr, const string &pattern, const int &spos=0);
  bool s_contains(string &instr, const string &pattern);
  void s_replace(string &instr, const string &pattern1, const string &pattern2, const int &spos=0);
  void s_replace_all(string &instr, const string &pattern1, const string &pattern2);
  bool s_tokenize(string &instr, string &substr, int &spos, const string &delim);
};
#endif

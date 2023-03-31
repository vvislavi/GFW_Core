/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
Class steers the initialization and calculation of n-particle correlations. Uses recursive function, all terms are calculated only once.
Latest version includes the calculation of any number of gaps and any combination of harmonics (including eg symmetric cumulants, etc.)
If used, modified, or distributed, please aknowledge the author of this code.
*/

#include "GFW.h"
GFW::GFW():
  fInitialized(false)
{
};

GFW::~GFW() {
  for(auto pItr = fCumulants.begin(); pItr != fCumulants.end(); ++pItr)
    pItr->DestroyComplexVectorArray();
};

void GFW::AddRegion(string refName, int lNhar, int lNpar, double lEtaMin, double lEtaMax, int lNpT, int BitMask) {
  if(lNpT < 1) {
    printf("Number of pT bins cannot be less than 1! Not adding anything.\n");
    return;
  };
  if(lEtaMin >= lEtaMax) {
    printf("Eta min. cannot be more than eta max! Not adding...\n");
    return;
  };
  if(refName=="") {
    printf("Region must have a name!\n");
    return;
  };
  Region lOneRegion;
  lOneRegion.Nhar = lNhar; //Number of harmonics
  lOneRegion.Npar = lNpar; //Number of powers
  lOneRegion.NparVec = vector<int>{}; //if powers defined, then set this to empty vector
  lOneRegion.EtaMin = lEtaMin; //Min. eta
  lOneRegion.EtaMax = lEtaMax; //Max. eta
  lOneRegion.NpT = lNpT; //Number of pT bins
  lOneRegion.rName = refName; //Name of the region
  lOneRegion.BitMask = BitMask; //Bit mask
  AddRegion(lOneRegion);
};
void GFW::AddRegion(string refName, int lNhar, int *lNparVec, double lEtaMin, double lEtaMax, int lNpT, int BitMask) {
  if(lNpT < 1) {
    printf("Number of pT bins cannot be less than 1! Not adding anything.\n");
    return;
  };
  if(lEtaMin >= lEtaMax) {
    printf("Eta min. cannot be more than eta max! Not adding...\n");
    return;
  };
  if(refName=="") {
    printf("Region must have a name!\n");
    return;
  };
  Region lOneRegion;
  lOneRegion.Nhar = lNhar; //Number of harmonics
  lOneRegion.Npar = 0; //If vector with powers defined, set this to zero
  lOneRegion.NparVec = vector<int>{};//lNparVec; //vector with powers for each harmonic
  for(int i=0;i<lNhar;i++) lOneRegion.NparVec.push_back(lNparVec[i]);
  lOneRegion.EtaMin = lEtaMin; //Min. eta
  lOneRegion.EtaMax = lEtaMax; //Max. eta
  lOneRegion.NpT = lNpT; //Number of pT bins
  lOneRegion.rName = refName; //Name of the region
  lOneRegion.BitMask = BitMask; //Bit mask
  AddRegion(lOneRegion);
};
int GFW::CreateRegions() {
  if(fRegions.size()<1) {
    printf("No regions set. Skipping...\n");
    return 0;
  };
  int nRegions=0;
  for(auto pItr=fRegions.begin(); pItr!=fRegions.end(); pItr++) {
    GFWCumulant *lCumulant = new GFWCumulant();
    if(pItr->NparVec.size()) {
      lCumulant->CreateComplexVectorArrayVarPower(pItr->Nhar, pItr->NparVec, pItr->NpT);
    } else {
      lCumulant->CreateComplexVectorArray(pItr->Nhar, pItr->Npar, pItr->NpT);
    };
    fCumulants.push_back(*lCumulant);
    ++nRegions;
  };
  if(nRegions) fInitialized=true;
  return nRegions;
};
void GFW::Fill(double eta, int ptin, double phi, double weight, int mask, double SecondWeight) {
  if(!fInitialized) CreateRegions();
  if(!fInitialized) return;
  for(int i=0;i<(int)fRegions.size();++i) {
    if(fRegions.at(i).EtaMin<eta && fRegions.at(i).EtaMax>eta && (fRegions.at(i).BitMask&mask))
      fCumulants.at(i).FillArray(ptin,phi,weight,SecondWeight);
  };
};
complex<double> GFW::TwoRec(int n1, int n2, int p1, int p2, int ptbin, GFWCumulant *r1, GFWCumulant *r2, GFWCumulant *r3) {
  complex<double> part1 = r1->Vec(n1,p1,ptbin);
  complex<double> part2 = r2->Vec(n2,p2,ptbin);
  complex<double> part3 = r3?r3->Vec(n1+n2,p1+p2,ptbin):complex<double>(0.,0.);
  complex<double> formula = part1*part2-part3;
  return formula;
};
complex<double> GFW::RecursiveCorr(GFWCumulant *qpoi, GFWCumulant *qref, GFWCumulant *qol, int ptbin, vector<int> &hars) {
  vector<int> pows;
  for(int i=0; i<(int)hars.size(); i++)
    pows.push_back(1);
  return RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows);
};

complex<double> GFW::RecursiveCorr(GFWCumulant *qpoi, GFWCumulant *qref, GFWCumulant *qol, int ptbin, vector<int> &hars, vector<int> &pows) {
  if((pows.at(0)!=1) && qol) qpoi=qol; //if the power of POI is not unity, then always use overlap (if defined).
  //Only valid for 1 particle of interest though!
  if(hars.size()<2) return qpoi->Vec(hars.at(0),pows.at(0),ptbin);
  if(hars.size()<3) return TwoRec(hars.at(0), hars.at(1),pows.at(0),pows.at(1), ptbin, qpoi, qref, qol);
  int harlast=hars.at(hars.size()-1);
  int powlast=pows.at(pows.size()-1);
  hars.erase(hars.end()-1);
  pows.erase(pows.end()-1);
  complex<double> formula = RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows)*qref->Vec(harlast,powlast);
  int lDegeneracy=1;
  int harSize = (int)hars.size();
  for(int i=harSize-1;i>=0;i--) {
  //checking if current configuration is a permutation of the next one.
  //Need to have more than 2 harmonics though, otherwise it doesn't make sense.
    if(i>2) { //only makes sense when we have more than two harmonics remaining
      if(hars.at(i) == hars.at(i-1) && pows.at(i) == pows.at(i-1)) {//if it is a permutation, then increase degeneracy and continue;
        lDegeneracy++;
        continue;
      };
    }
    hars.at(i)+=harlast;
    pows.at(i)+=powlast;
    //The issue is here. In principle, if i=0 (dif), then the overlap is only qpoi (0, if no overlap);
    //Otherwise, if we are not working with the 1st entry (dif.), then overlap will always be from qref
    //One should thus (probably) make a check if i=0, then qovl=qpoi, otherwise qovl=qref. But need to think more
    //-- This is not aplicable anymore, since the overlap is explicitly specified
    complex<double> subtractVal = RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows);
    if(lDegeneracy>1) { subtractVal *= lDegeneracy; lDegeneracy=1; };
    formula-=subtractVal;
    hars.at(i)-=harlast;
    pows.at(i)-=powlast;

  };
  hars.push_back(harlast);
  pows.push_back(powlast);
  return formula;
};
void GFW::Clear() {
  for(auto ptr = fCumulants.begin(); ptr!=fCumulants.end(); ++ptr) ptr->ResetQs();
};
GFW::CorrConfig GFW::GetCorrelatorConfig(string config, string head, bool ptdif) {
  //First remove all ; and ,:
  s_replace_all(config,","," ");
  s_replace_all(config,";"," ");
  s_replace_all(config,"| ","|");
  //If pT-bin is provided, then look for & remove space before "(" (so that it's clean afterwards)
  while(s_index(config," (")>-1) s_replace_all(config," (","(");
  //Then make sure we don't have any double-spaces:
  while(s_index(config,"  ")>-1) s_replace_all(config,"  "," ");
  vector<int> regs;
  vector<int> hars;
  int sz1=0;
  int szend=0;
  string ts, ts2;
  CorrConfig ReturnConfig;
  //Fetch region descriptor
  if(!s_tokenize(config,ts,szend,"{")) {
    printf("Could not find any harmonics!\n");
    return ReturnConfig;
  };
  szend=0;
  int counter=0;
  while(s_tokenize(config,ts,szend,"{")) {
    counter++;
    ReturnConfig.Regs.push_back(vector<int>{});
    ReturnConfig.Hars.push_back(vector<int>{});
    ReturnConfig.Overlap.push_back(-1); //initially, assume no overlap
    //Check if there's a particular pT bin I should be using here. If so, store it (otherwise, it's bin 0)
    int ptbin=-1;
    int sz2=0;
    if(s_contains(ts,"(")) {
      if(!s_contains(ts,")")) {printf("Missing \")\" in the configurator. Returning...\n"); return ReturnConfig; };
      sz2 = s_index(ts,"(");
      sz1=sz2+1;
      s_tokenize(ts,ts2,sz1,")");
      ptbin=stoi(ts2);
      ts.erase(sz2,(sz1-sz2+1));
      szend-=(sz1-sz2); //szend also becomes shorter
      //also need to remove this from config now:
      sz2 = s_index(config,"(");
      sz1 = s_index(config,")");
      config.erase(sz2,sz1-sz2+1);
    };
    ReturnConfig.ptInd.push_back(ptbin);
    sz1=0;
    //Fetch regions
    while(s_tokenize(ts,ts2,sz1," ")) {
      if(sz1>=szend) break;
      bool isOverlap=s_contains(ts2,"|");
      if(isOverlap) ts2.erase(0,1); //If overlap, remove the delimiter |
      int ind=FindRegionByName(ts2);
      if(ts2 == " " || ts2 == "") continue;
      if(ind<0) {
        printf("Could not find region named %s!\n",ts2.c_str());
        break;
      };
      if(!isOverlap)
        ReturnConfig.Regs.at(counter-1).push_back(ind);
      else ReturnConfig.Overlap.at((int)ReturnConfig.Overlap.size()-1) = ind;
    };
    string harstr;
    s_tokenize(config,harstr,szend,"}");
    int dummys=0;
    //Fetch harmonics
    while(s_tokenize(harstr,ts,dummys," ")) ReturnConfig.Hars.at(counter-1).push_back(stoi(ts));
  };
  ReturnConfig.Head = head;
  ReturnConfig.pTDif = ptdif;
  // ReturnConfig.pTbin = ptbin;
  return ReturnConfig;
};

complex<double> GFW::Calculate(int poi, int ref, vector<int> hars, int ptbin) {
  GFWCumulant *qref = &fCumulants.at(ref);
  GFWCumulant *qpoi = &fCumulants.at(poi);
  GFWCumulant *qovl = qpoi;
  return RecursiveCorr(qpoi, qref, qovl, ptbin, hars);
};
// complex<double> GFW::Calculate(CorrConfig corconf, int ptbin, bool SetHarmsToZero, bool DisableOverlap) {
//    vector<int> ptbins;
//    for(int i=0;i<(int)corconf.size();i++) ptbins.push_back(ptbin);
//    return Calculate(corconf,ptbins,SetHarmsToZero,DisableOverlap);
// }
complex<double> GFW::Calculate(CorrConfig corconf, int ptbin, bool SetHarmsToZero, bool DisableOverlap) {
  if(corconf.Regs.size()==0) return complex<double>(0,0); //Check if we have any regions at all
  // if(ptbins.size()!=corconf.Regs.size()) {printf("Number of pT-bins is not the same as number of subevents!\n"); return complex<double>(0,0); };
  complex<double> retval(1,0);
  int ptInd;
  for(int i=0;i<(int)corconf.Regs.size();i++) { //looping over all regions
    if(corconf.Regs.at(i).size()==0)  return complex<double>(0,0); //again, if no regions in the current subevent, then quit immediatelly
    ptInd = corconf.ptInd.at(i); //for i=0 (potentially, POI)
    if(ptInd<0) ptInd = ptbin;
    // int ptbin = ptbins.at(i);
    //picking up the indecies of regions...
    int poi = corconf.Regs.at(i).at(0);
    int ref = (corconf.Regs.at(i).size()>1)?corconf.Regs.at(i).at(1):corconf.Regs.at(i).at(0);
    int ovl = corconf.Overlap.at(i);
    //and regions themselves
    GFWCumulant *qref = &fCumulants.at(ref);
    GFWCumulant *qpoi = &fCumulants.at(poi);
    if(!qref->IsPtBinFilled(ptInd)) return complex<double>(0,0); //if REF is not filled, don't even continue. Could be redundant, but should save little CPU time
    if(!qpoi->IsPtBinFilled(ptInd)) return complex<double>(0,0);//if POI is not filled, don't even continue. Could be redundant, but should save little CPU time
    GFWCumulant *qovl=0;
    //Check if in the ref. region we have enough particles (no. of particles in the region >= no of harmonics for subevent)
    int sz1 = corconf.Hars.at(i).size();
    if(poi!=ref) sz1--;
    if(qref->GetN() < sz1) return complex<double>(0,0);
    //Then, figure the overlap
    if(ovl > -1) //if overlap is defined, then (unless it's explicitly disabled)
      qovl = DisableOverlap?0:&fCumulants.at(ovl);
    else if(ref==poi) qovl = qref; //If ref and poi are the same, then the same is for overlap. Only, when OL not explicitly defined
    if(SetHarmsToZero) for(int j=0;j<(int)corconf.Hars.at(i).size();j++) corconf.Hars.at(i).at(j) = 0;
    retval *= RecursiveCorr(qpoi, qref, qovl, ptInd, corconf.Hars.at(i));
  }
  return retval;
};

complex<double> GFW::Calculate(int poi, vector<int> hars) {
  GFWCumulant *qpoi = &fCumulants.at(poi);
  return RecursiveCorr(qpoi, qpoi, qpoi, 0, hars);
};
int GFW::FindRegionByName(string refName) {
  for(int i=0;i<(int)fRegions.size();i++) if(fRegions.at(i).rName == refName) return i;
  return -1;
};
//String processing:
int GFW::s_index(string &instr, const string &pattern, const int &spos) {
  return instr.find(pattern,spos);
};
bool GFW::s_contains(string &instr, const string &pattern) {
  return (s_index(instr,pattern)>-1);
};
void GFW::s_replace(string &instr, const string &pattern1, const string &pattern2, const int &spos) {
  int lpos = s_index(instr,pattern1,spos);
  if(lpos<0) return;
  instr.replace(lpos,pattern1.size(),pattern2);
};
void GFW::s_replace_all(string &instr, const string &pattern1, const string &pattern2) {
  int lpos=s_index(instr,pattern1);
  while(lpos>-1) { s_replace(instr,pattern1,pattern2,lpos); lpos=s_index(instr,pattern1,lpos); };
};
bool GFW::s_tokenize(string &instr, string &subs, int &spos, const string &delim) {
  if(spos<0 || spos>=(int)instr.size()) {spos=-1; subs=""; return false;};
  int lpos = s_index(instr,delim,spos);
  if(lpos<0) lpos=instr.size();
  subs = instr.substr(spos,lpos-spos);
  spos=lpos+1;
  return true;
}

#include <vector>
#include <stdlib.h>
#include "GFW.h"
using std::vector;
using std::make_pair;
using std::pair;
typedef vector<GFW::CorrConfig> CorrConfigs;
//Macro to read in angles. Returns a vector<vector<double> >, where each vector<double> corresponds to tracks in a given event.
//This is only for the test/example, as in "real" life, one would instead loop over all the (relevant) tracks in the event.
vector<vector<double> >ReadAngles(string inFile) {
  FILE *fin = fopen(inFile.c_str(),"r");
  double rVal;
  vector<vector<double> > retVec;
  bool isEOF=false;
  while(!isEOF) {
    int nTracks = rand() % 100 + 101; //Number of trakcs in each event anywhere between 100 and 200
    retVec.push_back(vector<double>{});
    int rvInd=retVec.size()-1;
    for(int nRead=0; nRead<nTracks; nRead++) {
      isEOF = fscanf(fin,"%lf, ",&rVal)<0;
      if(isEOF) break;
      retVec[rvInd].push_back(rVal);
    };
  };
  fclose(fin);
  return retVec;
};
pair<double, double> CalculateSingleConfig(GFW *inGFW, const GFW::CorrConfig &inConfig, const int &ptBin) {
  double dnx = inGFW->Calculate(inConfig,ptBin,true).real(); //First, set all harmonics to zero -- this gives a (weighted) number of combinations
  if(dnx==0) return make_pair(0,0); //If no pairs, then nothing to return
  double val = inGFW->Calculate(inConfig,ptBin,false).real(); //Then calculate the (unnormalised) correlation function
  pair<double, double> retVal = make_pair(val/dnx, dnx); //Here we already normalize the function by number of combinations. This is a little bit redundant, because when we average over events, we want to weight it again by dnx. But I leave it here as is for readability
  return retVal;
};
//Main loop utilizing GFW framework
int main() {
  //First, read in the input angles
  auto lEvents = ReadAngles("PhiAngles.dat");
  //Create GFW object
  GFW *fGFW = new GFW();
  //Initialize regions. Each region has identifier ("pos", "neg", "whatever"), eta range, number of pT bins (>0), and bit mask (>0). The preferred way to do this is:
  fGFW->AddRegion("FullReg",-1.0,1.0,1,1);
  //Here the dimensions of Q-vector are calculated automatically based on the configurations that are fetched when calling GetCorrelatorConfig. Other methods are left as legacy:
  //Legacy no. 1: powers are provided as a vector<int>. For their calculation, see GFWPowersArray.
  vector<int> PowersAsVector = {5, 0, 4, 4, 3, 0, 3};
  fGFW->AddRegion("FullReg",PowersAsVector,-1,1,1,1);
  //Legacy no. 2: powers as array. Similar to previous case, but now also need to specify number of entries
  int PowersAsArray[] = {3, 0, 2, 2, 3, 0, 3};
  int nPowersAsArray=7;
  fGFW->AddRegion("PosSide",nPowersAsArray,PowersAsArray,0,1,1,1);
  //Lgeacy no. 3: least CPU-efficient, b/c many unused Qs will be calculated. Maximum power for each harmonic considered
  int PowersAsInt=3; //Maximum power is 3
  int nPowersAsInt=7; //Maximum harmonic is 6 (-> 3+3)
  fGFW->AddRegion("NegSide",nPowersAsInt,PowersAsInt,-1,0,1,1);
  //Let us also define a separate region of POI's, let's say, pT-differential and for specific spiecies
  int nPtBins=3; //Let's assume 3 different pT bins
  int BitMask=2; //Different bit mask. This is so that we can first identify the particle and then can decide whether to fill the relevant vector or not
  fGFW->AddRegion("PIDNeg",-1,0,nPtBins,BitMask);
  //Next, define what we want to calculate. This is done in terms of correlator configurations that are defined as "[POI_R1] REF_R1 [ | Overlap_R1 ] {Harmonics_R1} [ [POI_R2] REF_R2 [ | Overlap_R2] {Harmonics_R2} ]", where terms in [] brackets are optional, see README for more info. If we want to calculate <2, 2, -2, -2> for single event and 2-subevets (Positive side -- Negative side), the configurations will look as following:
  CorrConfigs configs{}; //First make a vector where we store all the configs
  configs.push_back(fGFW->GetCorrelatorConfig("FullReg {2 2 -2 -2}","FR2222",false)); //For har=2, 4PC, full event
  configs.push_back(fGFW->GetCorrelatorConfig("NegSide {2 2} PosSide {-2 -2}","SE2222",false)); //For har=2, 4PC, 2 subevent
  //In general, if no overlap and no POI are specified, then the reference region will be used to calculate the overlap (when phi1=phi2). However, if only one particle is taken from a regions ("PosSide {2}"), then there are no auto-correlations and thus no overlap
  //Next, let's make a configuration for PID particle which is also pT differential:
  configs.push_back(fGFW->GetCorrelatorConfig("PIDNeg NegSide | PIDNeg {2 2} PosSide {-2 -2}","PID2222",true));
  //In this case, on the left-hand side we have "PIDNeg NegSide | PIDNeg {H1 H2}". The first two terms indicate that first harmonic H1 is taken from PIDNeg region and the second one from NegSide. There might be some overlap between PIDNeg and NegSide, which in our case is exactly PIDNeg, and this is specified via a | symbol. The last argument "true" indicates that the configuration is pt-dependent, which will be storred in the configurator and later be used to decide whether we want to loop over pt bins or not.
  //After GFW is configured, let's initialize it. In principle there are multiple checks to make sure that it's initialised runtime, but it's also healthy practice to do it ourselves (in which case, the bool checks can be removed from Fill() and Calculate() methods to save a little bit of time)
  //IMPORTANT! If one uses the preferred AddRegion method (line 43), CreateRegions MUST be called after all the GetCorrelatorConfigs() calls. This is because GFW _needs_ to know which sets of harmonics will be considered.
  fGFW->CreateRegions();
  //Then let's define some storage container to keep track of event-by-event values. Normally one would use something like TProfile or sth, but here we just keep a simple vector of pairs (value and normalization)
  vector<pair<long double, long double> > storage = {make_pair(0,0),make_pair(0,0),make_pair(0,0),make_pair(0,0),make_pair(0,0)};
  //Finally, start calculations:
  for(auto PhiAngles:lEvents) { //Event loop
      fGFW->Clear(); //First, reset all the Q-vectors
      for(auto phi:PhiAngles) {
        //Lets generate a random eta, random pt, and random probability that it's a proton:
        double eta = (rand() % 20)*0.1 -1; //Eta between -1 and 1
        int ptInd  = rand()%3; //pt bin between 0 and 2
        bool isPID = (rand() % 100) < 30; //30% chance that this is a proton
        int BitMask= 1;
        if(isPID) BitMask+=2; //If this is a PID particle, then also need to fill the relevant Q vectors, see line 60
        double weight=1; //Some weight that could/would come from eg NUA/NUE
        fGFW->Fill(eta,ptInd,phi,weight,BitMask);
      }
      //Now that GFW has been filled, calculate the correlations defined in lines 61, 62, and 65:
      //First let's calculate FR2222:
      pair<double, double> fr22 = CalculateSingleConfig(fGFW,configs[0],0);
      storage[0].first += fr22.first*fr22.second; storage[0].second+=fr22.second;
      //Second, let's do SE2222
      pair<double, double> se22 = CalculateSingleConfig(fGFW,configs[1],0);
      storage[1].first += se22.first*se22.second; storage[1].second+=se22.second;
      //Then calculate the three pT bins for PID:
      for(int i=0;i<3;i++) {
        pair<double, double> pid22 = CalculateSingleConfig(fGFW,configs[2],i);
        storage[2+i].first += pid22.first*pid22.second; storage[2+i].second+=pid22.second;
      };
  }
  //Perform the normalization:
  for(auto &pr: storage) pr.first/=pr.second;
  //Print out calculated values:
  printf("%s: %Lf\n",configs[0].Head.c_str(),storage[0].first);
  printf("%s: %Lf\n",configs[1].Head.c_str(),storage[1].first);
  //And also the pT-dif PID:
  for(int i=0;i<3;i++) printf("%s_pt%i: %Lf\n",configs[2].Head.c_str(),i,storage[2+i].first);
  return 0;
}

C++ Implementation of Generic Framework for N-particle correlations
Author: Vytautas Vislavicius
Extention of https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.

Disclaimer:
The framework was originally written by V.V. and implemented to AliPhysics. With the move to O2, there has been a (somewhat) outdated port to O2, which dropped AliPhysics dependences, but ROOT dependences remained. Before I leave, I wanted to finalise all the unfinished ideas for the GFW and make it AliPhysics/O2/ROOT independent. The current version thus has no ROOT dependences, which means that it can be both easily moved to O2 (essentially, copy-paste -- I left all the "old" functions for compatibility with existing tasks), rivet, or any other framework that is build upon c++. In addition, it can also be used stand-alone without the necessity of e.g. ROOT.
Disclaimer 2: There is also a branch ROOT_Types, where instead of standard c++ prescription, variables have ROOT types (Double_t, Int_t, etc.). However, that version lacks few latest features, in particular, automatic power calculation that is now included into the c++ version.

Introduction:
My vision of GFW is to have a reliable, fast, and easy-to-use framework to calculate any order of (genuine) particle correlations with as many harmonics, subevents, POIs, etc. The framework uses recursive formula, but also keeps track of the degeneracy of the calculation in order to ensure that same terms, for a given configuration, are calculated only once. A guide how to use it:

-- Initialization:
   -- Create an object as any other object in C++:
   AliGFW *fGFW = new AliGFW();

   -- GFW works in terms of "regions". Each region corresponds to 3-dimensional array (pt, harmonic, and power) of Q-vectors, covering a particular eta region (from EtaMin to EtaMax) and having a pre-set bit mask (BitMask). When calling GFW::Fill(eta, ..., bitmask) method, only values with EtaMin < eta < EtaMax and (bitmask&BitMask)!=0 will be added to the corresponding Q-vector. Regions are created as:

   fGFW->AddRegion(string refName, double lEtaMin, double lEtaMax, int lNpT, int BitMask)

      -- For legacy reasons I left couple of old methods, where powers of the Q-vector have to be explicitly specified. However, this is a preferred method, as it calculates the required powers, so you don't need to spend time figuring it out yourself and possibly making errors
      -- Note that lNpT _has_ to be larger than 0, as it is one of the dimensions of Q vector. If we don't care about the pT, then lNpT=1.
      -- refName is a string that will only be used to find the corresponding region when building correlator configurations (next on the list)
      -- An example of this would be to add 2 regions in positive ("pos") and negative ("neg") eta to calculate e.g. 4PC with 2 subevents (gap 1)
      fGFW->AddRegion("pos",0.5,0.8,1,1);
      fGFW->AddRegion("neg",-0.8,-0.5,1,1);

   -- The calculation of correlations is done using correlator configurations. These are light-weight structures that "know" which Q-vectors need to be multiplied by which other Q-vectors, etc. To get this configuration, after adding all the regions, we call:

   fGFW->GetCorrelatorConfig(string config, string head = "", bool ptdif=false)

      -- The first argument is the most important and defined what we want to calculate. Let's say we have added positive and negative regions (see example few lines above). To calculate 4 particle correlation of e.g. mixed harmonics, ( <2, 3 | -2, -3>) from the two subevents, we make a configuration:

      AliGFW::CorrConfig my_config1 = fGFW->GetCorrelatorConfig("pos {2 3} neg {-2 -3}","SomeName",false);

      -- In principle, one can add as many subevents as one pleases. For instance, let's add another Q-vector at midrapidity (|eta|<0.5) and define 4-particle correlation with 3-subevents <3 | -3, -2 | 2>:

      fGFW->AddRegion("mid",-0.5,0.5,1,1);
      AliGFW::CorrConfig my_config2 = fGFW->GetCorrelatorConfig("pos {3} mid {-3 -2} neg {2}","SomeOtherName",false);

      -- The generic prescription of a single subevent is "[POI] Ref [ | Overlap] {H1[, H2, H3...]}" where terms in [] are optional. If we want to calculate a differential correlation, we need to specify the Q-vector that contains our particles of interest (POI). First, add another region that will contain POIs (let's say, some PID particle in the negative eta):

      fGFW->AddRegion("PID",-0.8,-0.5,1,2);
      and then calculate e.g. <2', 2 | -2 -2>, where 2' denotes second harmonic of PID particle:
      AliGFW::CorrConfig my_configPID = fGFW->GetCorrelatorConfig("PID neg | PID {2 2} pos {-2 -2}","PIDCorr",false);

      -- Note that here for the first subevent we have POI(="PID"), Ref(="neg"), and Overlap(="PID"). The overlap here has to be explicitly specified if there is any overlap between "PID" and "neg" regions. In this particular case, "PID" is a proper subset of "neg" (i.e. all "PID" particles are also part of the "neg" region), and so the overlap is also "PID". This is not always the case though, think e.g. "neg" particles are in pT-range [0.2-5] GeV/c and "PID" particles are in pT-range [3-7] GeV/c. Then, only "PID" particles in range [3-5] GeV/c will be in the overlap. In such case, one has to define a separate region to keep track of the "PID"-"neg" overlap. On the other hand, if there is _no_ overlap between "PID" and "neg" (e.g. "PID" is in pT-range [5-7] GeV/c), then _no_ overlap should be specified.

      -- In general, for any subevent with _more_ than 1 harmonic, the autocorrelation (overlap) is removed in 3 cases:
         1) Only reference region is provided (ie only one region for subevent and no POI's/overlaps)
         2) POI and ref are provided, and they _are_ the same regions (e.g. "neg neg {2 2} pos {-2 -2}")
         3) POI, ref, and overlap are provided

  -- Once all regions have been added _and_ all the correlator configurations have been fetched, we need to initialize the Q-vectors. This is done by:
    fGFW->CreateRegions();
    -- It is _extremely_ important that this method is called _only_ after all the regions have been created and all the correlator configurations are fetched. This is because when calling GetCorrelatorConfig(), the configuration of harmonics is stored internally in GFW, and then all the configurations are used to calculate the relevant power arrays for all the calculations. If one builds a new correlator config _after_ calling CreateRegions(), chances are you will have some powers of Q vector that are not available.
    -- I have also removed the checks on initialization from Fill() and Calculate() methods, because they take time and it _has_ to be users responsibility to initialize the GFW _before_ running calculations!

-- Running
  -- In the event loop (for each event), there are few things we need to care about:
    -- Clearing up all the vectors. This has to be done before filling up the GFW with new tracks
    fGFW->Clear();

    -- Fill up the vector in the track loop:
    fGFW->Fill(eta,ptInd,phi,weight,mask,secondWeight)
      -- eta is the pseudorapidity of track, only relevant regions get filled with this track
      -- ptInd defines which pt bin of the GFW should be filled (starting from 0). If the region has only 1 pt bin, then this number doesn't matter and will internally be overriden by 0
      -- phi is azimuthal angle
      -- weight is whatever weight (NUA, NUE, etc.) we want to add
      -- mask is a bitmask. Only regions with overlapping masks will be filled. This is relevant e.g. if we want to do PID: if the particle is PID, then we fill it with mask=3 (so it goes to both ref & POI vectors, since 3&1=1 and 3&2=2), otherwise mask=1 (so only ref region is filled).
      -- secondWeight is an advanced feature that will be irrelevant for most. It's role is to do a proper weight counting when e.g. a PID particle goes with one weight to reference region, and a different weight to POI. I will add more info on this later on, but rule of a thumb is, if it's a PID particle, it should have the same weight when going to POI and ref, and then secondWeight should have its default value, -1

    -- After finishing the track loop, Q-vectors are all filled and we can now calculate the N-particle correlations that we have defined in the correlation configurations. This is done by:

    fGFW->Calculate(CorrConfig corconf, int ptbin, bool SetHarmsToZero);
      -- The return value is a complex number, complex<double>
      -- ptbin is the pT bin for pT-differential regions (e.g. POI), and 0 otherwise (for regions with only 1 pT bin, it's automatically overridden)
      -- SetHarmsToZero is a bool argument that controls, as the name suggests, whether the harmonics should be set to zero. The total correlation function is calculated with SetHarmsToZero=false (default value), and it's normalization is calculated with SetHarmsToZero=true
      -- So to calculate e.g. <2, 3 | -2, -3>, as we had in one of the previous examples, we do:

      double Norm  = fGFW->Calculate(my_config1,0,true).real(); //Normalization
      double Value = fGFW->Calculate(my_config1,0,false).real(); //Correlation function
      Value/=Norm; //Normalize the correlation to number combinations. Also, one should check that Norm>0

    -- It is up to the user to decide where and how one wants to store/average these values; a typical choice in ROOT would be a TProfile, while for more simple approaches, one can just keep a vector of doubles. However, note that the event-averaging should be done by weighting each value by the number of combinations, that is:
    <<2, 3 | -2, -3>> = Sum(Value*Norm)/Sum(Norm)


Closing remarks:
-- I also include a Test.C file with an example of GFW in action. I added a ton of comments there, so you can look through that; the macro also compiles and runs (just type "make" and then "./Test")
-- You might ask why do "head" and "ptdif" arguments when making a correlator configuration. These have been added for simplicity when calling GFW::Calculate(...) function. In particular, if you have a whole array of CorrConfigs, you can fill a respective bin in e.g. TProfile that is called the same as "head", and you can also check whether the configuration is pT-differential (so you can have another loop over all the pT bins) or not, without writing explicit cases for each configuration.
-- There is also a new feature of specifying which pT bin should be used for each region. This is specified in parenthesis in the configurator as e.g. "PID (1) PID (2) {2 2} pos {-2 -2}", to correlate two PID particles from 2 different pT bins with reference. This can be useful when e.g. calculating vn-square bracket. I have not tested the feature excessively yet though.

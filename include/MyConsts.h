#ifndef MyConsts_h
#define MyConsts_h

#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// JK 2020
// spacial dimensions:
const int kNspace = 3;

// JK 2019

const double kGeV = 1e-3;
const double kEpsilon = 1e-5;
const double kInfty = 999e9;
const double kmW = 80.385;
const double kmZ = 91.1876;
const double kmt = 172.5;

const double kbtagFakeEff = 0.01;

/*
 * was:
const double kTau21WCut2 = 0.5;
//const double kTau32WCut1 = 0.5;
const double kmWcut2 = 125.;
const double kmWcut1 = 55.;
const double kTau32tCut = 0.6; // 0.75
const double kTau21tCut1 = 0.30; // 0.75
const double kTau21tCut2 = 0.70; // 0.75
const double kmtcut1 = 140.;
const double kmtcut2 = 215.;
*/

const double kTau21WCut1 = 0.10;
const double kTau21WCut2 = 0.60;
const double kTau32WCut1 = 0.50;
const double kTau32WCut2 = 0.85;
const double kmWcut1 = 66.; // was 50
const double kmWcut2 = 105.; // was 110

const double kTau21tCut1 = 0.30;
const double kTau21tCut2 = 0.70;
const double kTau32tCut1 = 0.30;
const double kTau32tCut2 = 0.80;
const double kmtcut1 = 130.;
const double kmtcut2 = 195.;


const int kntmax = 2; // 4
const int knwmax = 2; // 4; w or semiboosted??
const int knbmax = 4; //

// to add also inclusive cases like any semiboosted??
/*
const TString kSelections[] = {"4B0S", "3B1S", "2B2S", "1B3S", "0B4S",  // 4 B/S tops found
                               "3B0S", "2B1S", "1B2S", "0B3S",          // 3 B/S tops found
                               "2B0S", "1B1S", "0B2S",                  // 2 B/S tops found
                               "1B0S", "0B1S",                          // 1 B/S tops found
                               "0B0S"                                   // 0 B/S tops found
                              };
*/

const TString kSelections[] = {
                                // 2 B/S tops found
                                "2B0S",  // 2B incl., i.e. >=2B
                                "1B1S",  // 1B excl., 1S excl!
                                "0B2S",  // 0B excl., 2S inclusive == >=2!
                              };

enum kJetTypes {kallJets, kbJets, knonbJets, kTopTagJets, kWTagJets};

enum kLJetChoice { kResolved, // for historical resolved reasons
                   kDefault,
                   kSoftDroppedP4,
                   kPrunedP4,
                   kTrimmedP4};

enum kPseudotopType { kNonePseudo, kOldWhad, kStandard, kCloseMt, kSameMt, kTwoStep, kTwoStepII, kBestBsAndNu};

enum kLjets { kUndef, kEjets, kMujets};

// 27.5.220:
struct kGlobalVars {
    double met;
    double metphi;
    double HTj;
    double HTJ;
    double SumMJ;
    double SumMj;
    double MjVis; // visible mass from small jets
    double MJVis; // visible mass from large Jets
    // aplanarity and sphericity using small or large jets:
    double jApla;
    double jSpher;
    double JApla;
    double JSpher;
};

struct kPseudotop {

  TMyLorentzVector lepton;
  TMyLorentzVector neutrino;

  TMyLorentzVector pseudotophadron;
  TMyLorentzVector pseudotoplepton;

  TMyLorentzVector pseudoWhadron;
  TMyLorentzVector pseudoWlepton;

  TMyLorentzVector bhadron;
  TMyLorentzVector blepton;

  TMyLorentzVector j1;
  TMyLorentzVector j2;
    
  TMyLorentzVector pseudottbar;

  double Pout[2];
  double CosThetaStar[2];
  double deltattbar; // opening angle in lab
  double Yboost;
  double DeltaPhi;
  double Chittbar;

  // ratio variables
  double DiTopPoutRel[2];
  double TopPtRel[2]; // lep, had; or: 1, 2
  double WPtRel[2]; // lep, had; or: 1, 2
  double DiTopPtRel; // ratio to mtt
  double DiTopPtGeo; // ratio to sqrt(pTt1*pTt2)
  double DiTopMassGeo;
  double DiTopPoutGeo[2];
  double Rttbar;
  // signed discriminants^{1/6} of the quadratic equations for the neutino pz solution
  // based on the mlnu = mW or mthad = mtlep conditions
  double DmW;
  double Dmt;
  

};

    
#endif

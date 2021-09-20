//////////////////////////////////////////////////////////
// NOT ANY MORE:
// This class has been automatically generated on
// Thu Apr 27 16:00:22 2017 by ROOT version 5.34/03
// from TTree DelphesTree/Analysis tree
// found on file: /media/data/qitek/DelphesTreeOut/out_boosted_AtKt4_and_10_ljetsUPGRADE_run_06_zprime_850_alljets.root
// as Jiri Kvita moved from MakeClass to the Delphes TreeReader on 14th Oct 2020 (covid+ day of 8k in Czechia...)
//////////////////////////////////////////////////////////

/*
  DOCUMENTATION:
  https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
  [JK 13.2.2019]
 */

#ifndef DelphesTree_h
#define DelphesTree_h

// 14.10.2020
#include "DelphesClasses.h"
#include "ExRootTreeReader.h"

#include <fstream>

#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TF2.h"
#include "TChain.h"
#include "TFile.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TRandom3.h"

// Header file for the classes stored in the TChain if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// JK:
#include "TMyLorentzVector.h"
#include <iostream>
#include <vector>
#include "MyConsts.h"
#include "HistoMaker.h"

using namespace std;

class DelphesTree {
public :

    TChain          *fChain;   //!pointer to the analyzed TChain or TChain

    // JK:
    bool m_isSignal;
    TString m_sampleTag;
    bool m_CorrectLJES;
    bool m_CorrectSJES;
    bool m_FillJESclosure;
    TFile *m_JES_calibFile;
    TF2 *m_JES_LJets_response;
    TF2 *m_JES_SJets_response;

    TH2D *m_SelMigra;
    TH1D *m_hSelDet;
    TH1D *m_hSelPtcl;

    std::map<int,TH1D*> m_xPdfHistos;

    TH1D *m_InitialPartonsFracHisto;
    TH1D *m_InitialPartonsShat_gg;
    TH1D *m_InitialPartonsShat_qg;
    TH1D *m_InitialPartonsShat_qq;

    TH1D *m_CutValsHisto;

    TH1D *m_CutFlowHistoDet;
    TH1D *m_CutFlowHistoPtcl;

    TH1D *m_TopTagttruthMatchedPt_ptcl;
    TH1D *m_WTagWtruthMatchedPt_ptcl;
    TH1D *m_TopTagnontMatchedPt_ptcl;
    TH1D *m_WTagnonWMatchedPt_ptcl;
    TH1D *m_notTopTagttruthMatchedPt_ptcl;
    TH1D *m_notWTagWtruthMatchedPt_ptcl;
    TH1D *m_notTopTagnontMatchedPt_ptcl;
    TH1D *m_notWTagnonWMatchedPt_ptcl;

    TH1D *m_TopTagttruthMatchedPt_reco;
    TH1D *m_WTagWtruthMatchedPt_reco;
    TH1D *m_TopTagnontMatchedPt_reco;
    TH1D *m_WTagnonWMatchedPt_reco;
    TH1D *m_notTopTagttruthMatchedPt_reco;
    TH1D *m_notWTagWtruthMatchedPt_reco;
    TH1D *m_notTopTagnontMatchedPt_reco;
    TH1D *m_notWTagnonWMatchedPt_reco;

    void MakeHistoMakers(int nReplicas);
    void ResetSelectionDictionaries();
    std::map<TString, bool> m_SelectDict;

    TString m_AnyPassedTag;
    TString m_NothingPassedTag;

    // passed ptcl, binned in ptcl
    HistoMaker *m_hmaker_ptcl;
    // passed det, binned in det
    HistoMaker *m_hmaker_det;
    // passed reco cuts, binned in parton
    HistoMaker *m_hmaker_det_parton;
    // passed det && ptcl, binned in ptcl
    HistoMaker *m_hmaker_detptcl_ptcl;
    // passed det && ptcl, binned in det
    HistoMaker *m_hmaker_detptcl_det;

    // parton
    HistoMaker *m_hmaker_parton;

    // migrations particle-detector
    HistoMaker *m_hmaker_migra_ptcl;
    // migrations parton-detector
    HistoMaker *m_hmaker_migra_parton;
    // migrations parton-particle
    HistoMaker *m_hmaker_migra_partonptcl;

    int m_ndetbTags;
    int m_ndetWTags;
    int m_ndettopTags;

    int m_nptclbTags;
    int m_nptclWTags;
    int m_nptcltopTags;

    TString m_detPassedTag;
    TString m_ptclPassedTag;

    TRandom3 *m_rand;

    bool SkipTopologies(TString sel);
    TString MakePassedTag(int nbTags, int nWTags, int ntopTags);

    double CorrectMetForMuons(double &metx, double &mety, double pTmuMin);
    double ComputeMetPhi(double metx, double mety);

    vector<TMyLorentzVector*> RemoveParticleJetsOverlapWithLeadLepton(kLjets ljets, vector<TMyLorentzVector*>& jets, double ptcut, bool debug = false);
    bool DoNotOverlap(double eta, double phi, vector<TMyLorentzVector*>& jets, bool debug = false);
    void MakePout(kPseudotop &pst);
    void MakeGlobalVars(bool ptcl, kGlobalVars &globalVars, vector<TMyLorentzVector*> &jets, vector<TMyLorentzVector*> &ljets);
    int CheckPstMatching( kPseudotop pptcl, kPseudotop pdet);
    bool IsBtaggedPtclLevel(double jpt, double jeta, double jphi, double jm);
    TMyLorentzVector DressByPhotons(TMyLorentzVector lep);
    vector<TMyLorentzVector*> MakeJets(bool DetLevel, bool LJets, double ptcut, double etacut,
                                       bool CorrectJES,
                                       bool runFakeBtag = true,
                                       kLJetChoice LJetChoice = kDefault);

    vector<TMyLorentzVector*> MakeGenWs(double ptcut, double etacut, int status = 52);
    vector<TMyLorentzVector*> MakeGenTops(double ptcut, double etacut, int status = 52);
    vector<TMyLorentzVector*> MakeGenDiTops(vector<TMyLorentzVector*> GenTops);

    vector<TMyLorentzVector*> MakeTops(vector<TMyLorentzVector*> &jets, vector<TMyLorentzVector*> &ljets, double mindrbjet, int debug = 0);
    vector<TMyLorentzVector*> MakeDiTops(vector<TMyLorentzVector*> &topss, int debug = 0);
    vector<TMyLorentzVector*> MakeFourTops(vector<TMyLorentzVector*> tops, int debug = 0);

    void InitCutVals(TString LjetTypeStr);

    // info/dump
    ofstream *m_asciifile = 0;
    double m_EventWeight;
    double m_xsectWeight;
    bool m_ApplyXsectWeights;
    int m_debug;
    bool m_dumpASCII;
    bool m_doPdfAnalysis;
    
    // JET CUTS!
    double m_Ljetacut;
    double m_jetacut;
    double m_Ljptcut;
    double m_jptcut;
    // top and W tagging study
    double m_drtcut;
    double m_drWcut;
    double m_mindrbjet;

    bool m_runFakeBtag;
    kLJetChoice m_LJetChoice;

    // to veto an isolated lepton:
    double m_leptonEtaCut;
    double m_leptonPtCut;

    void FillTagEffHistos(HistoMaker *m_hmaker,
                          vector<TMyLorentzVector*> &detLjets,
                          vector<TMyLorentzVector*>& partonTops, vector<TMyLorentzVector*> &partonWs);
    void DumpAscii(HistoMaker *m_hmaker,
                   vector<TMyLorentzVector*>& detjets, vector<TMyLorentzVector*> &detLjets, vector<TMyLorentzVector*> &partonTops, vector<TMyLorentzVector*> &partonWs);

    bool MakeWjets(bool DetLevel, vector<TMyLorentzVector*>& jets);
    bool MakeBs(bool DetLevel, TMyLorentzVector lepton, vector<TMyLorentzVector*>& jets);
    
    bool MakeLepton(kLjets ljets, bool DetLevel, vector<TMyLorentzVector*>& jets);
    bool MakeNeutrino(bool DetLevel, bool SwapBs = false);
    bool MakePseudotops(bool DetLevel, vector<TMyLorentzVector*>& jets);

    kLjets m_ljets;
    
    // pseudotops:
    kPseudotopType m_PseudotopType;
    kPseudotop m_pseudo_det;
    kPseudotop m_pseudo_ptcl;
    kPseudotop m_pseudo_parton;

    kGlobalVars m_globalVars_ptcl;
    kGlobalVars m_globalVars_det;

    // some bools
    bool m_MadePstDetDiTops;
    bool m_MadePstPtclDiTops;
    bool m_MadePstPartonDiTops;

    DelphesTree(TString path, TString pattern, TString weightStr);
    virtual ~DelphesTree();
    virtual void     Init(TChain *tree);

    // JK
    void LoopBoosted(TString outtag, TString LjetTypeStr, int nReplicas);
    virtual void LoopResolved(TString outtag, TString channel, TString pseudo);


    // Branches for Delphes Tree Reader
    // JK 14.10.2020
    ExRootTreeReader *m_treeReader;

    TClonesArray *m_branchEvent;

    //TClonesArray *m_branchParticle;
    TClonesArray *m_branchElectron;
    TClonesArray *m_branchPhoton;
    TClonesArray *m_branchMuon;
    TClonesArray *m_branchJetJES;
    TClonesArray *m_branchLJet;
    TClonesArray *m_branchMet;

    TClonesArray *m_branchGenJet;
    TClonesArray *m_branchGenLJet;
    TClonesArray *m_branchGenMet;
    TClonesArray *m_branchGenElectron;
    TClonesArray *m_branchGenMuon;

    TClonesArray *m_branchGenPhoton;
    TClonesArray *m_branchGenBhadrons;
    TClonesArray *m_branchGenTop;
    TClonesArray *m_branchGenW;


};

#endif


#include "HistoHolder.h"

#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TMyLorentzVector.h"

#include "KinemUtils.h"

#include "MyConsts.h"

#include <map>
#include <vector>

using std::map;
using std::vector;

#ifndef HistoMaker_HPP
#define HistoMaker_HPP

const int kMaxJetsToFill = 4;

class varProps {
public:
  varProps(TString tit, int i, TString un, TString fstr) {
    this -> title = tit;
    this -> id = i;
    this -> unit = un;
    this -> fillstr = fstr;
  };
  varProps() {};
  ~varProps() {};
  TString title;
  int id;
  TString unit;
  TString fillstr;
};


class HistoMaker
{

    // class to hold histograms for each of the reco, particle, reco && particle (binned in reco or particle:) levels

private:
    map<TString,HistoHolder*> m_holders; // map of histo holders for each cut level
    TString m_level;
    // for migrations:
    TString m_level1;
    TString m_level2;

    std::map<TString,varProps> m_vars;
    map<TString, double*> m_PhysObjBins;
    map<TString, int> m_nPhysObjBins;

    //double *m_ptbins;
    int m_nptbins;
    //double *m_mbins;
    int m_nmbins;

    int m_nTauBins;
    double m_Taumin;
    double m_Taumax;
    
    int m_nReplicas;
    double *m_replicaWeights;


public:

    void SetReplicaWeights(double *weights);

    HistoMaker(TString level, int nReplicas = 0);
    ~HistoMaker();

    TString GetLevel();

    void AddCutLevel(TString cutname, bool jetsOnly = false);

    bool PassedJetType(kJetTypes jettype, TMyLorentzVector *jet);
    void MakeSingleObjHistos(TString ObjName, HistoHolder *hold);

    void MakeDijetHistos(TString cutname);
    void MakeControlHistos(bool makeLjet = false, bool makeDijet = false);
    void MakeAllJetHistos(bool makeReplicas);
    void MakeGlobalHistos(HistoHolder *hold);
    void MakePartonHistos();
    void MakeJetHistos(TString JetType, HistoHolder *hold, bool MakeSubs, bool makeReplicas = false);
    void MakeGlobalMigrations(TString cutname, TString objname, TString level1, TString level2);

    void MakeJetMigrations(TString cutname, TString JetType, TString level1, TString level2,
                           bool FillSub = true, int AddMore = 4);

    void FillAlljetMigrations(TString cutname,
                              vector<TMyLorentzVector*> ptclLjets, vector<TMyLorentzVector*> ptcltops, vector<TMyLorentzVector*> ptclDiTops, vector<TMyLorentzVector*> ptclFourTops,
                              vector<TMyLorentzVector*> detLjets, vector<TMyLorentzVector*> dettops, vector<TMyLorentzVector*> detDiTops, vector<TMyLorentzVector*> detFourTops,
                              double weight);

    void MakeAlljetMigrations(TString cutname, TString level1, TString level2);

    void MakeJESHistos(TString cutname, bool MakeClosure, bool MakeJES);
    void FillJESHistos(TString cutname, vector<TMyLorentzVector*> ptcljets, vector<TMyLorentzVector*> detjets, bool IsClosure, bool isLJets, double weight);

    void MakeSpecialHistos(HistoHolder *hold, TString jname, TString jtag, bool makeReplicas);

    void FillSingleObjectHistos(TString cutname, TString JetType, TMyLorentzVector& obj, double weight);

    void FillJetHistos(TString JetType, TString cutname, vector<TMyLorentzVector*> jets, double weight, kJetTypes jettype, bool FillSubs = false);
    void FillJetMigrations(TString cutname, TString jettag, vector<TMyLorentzVector*> ptcljets,
                           vector<TMyLorentzVector*> detjets,
                           double weight, int AddMore, kJetTypes jettype = kallJets, bool FillSubs = false);
    void FillMigrations(TString cutname, vector<TMyLorentzVector*> ptclLjets, vector<TMyLorentzVector*> ptcljets, vector<TMyLorentzVector*> ptclDijets,
                        vector<TMyLorentzVector*> detLjets, vector<TMyLorentzVector*> detjets, vector<TMyLorentzVector*> detDijets,
                        double weight);

    void MakeLjetHistos(TString cutname);
    void FillGlobalMigrations(TString cutname, TString objname, kGlobalVars globalVars_ptcl, kGlobalVars globalVars_det, double weight);
    void FillSpecialMigrations(TString cutname, TString objname, kPseudotop pseudo_ptcl, kPseudotop pseudo_det, double weight);
    void MakeLjetMigrations(TString cutname, TString level1, TString level2);
    void FillTTbarSpecialHistos(TString cutname, TString DiTopName, kPseudotop pseudo, kGlobalVars GlobalVars, double weight, bool fillReplicas);
    void FillGlobalHistos(TString cutname, kGlobalVars GlobalVars, double weight);
    void FillLjetHistos(TString cutname, kPseudotop pseudo, kGlobalVars GlobalVars, double weight, bool fillReplicas);
    void FillLjetMigrations(TString cutname, kPseudotop pseudo_ptcl, kPseudotop pseudo_det, double weight);
    void FillLjetMigrationsParton(TString cutname, kPseudotop pseudo_parton, kPseudotop pseudo, double weight);

    void MakeSingleDRHistos(TString cutname);
    void MakeSingleTagKinemHistos(TString cutname, TString objname);
    void MakeTaggingHistos();
    void FillSingleTagKinemHistos(TString objname, TMyLorentzVector *jet, double weight);
    void FillSingleTagDRHistos(mindrdetadphi drt, mindrdetadphi drw, double pt, double weight);


};

#endif

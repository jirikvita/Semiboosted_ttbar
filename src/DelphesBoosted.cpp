#include "DelphesTree.h"
#include "HistoMaker.h"

#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "MyConsts.h"
#include "KinemUtils.h"

// jk Jan 2019
// __________________________________________________________________

bool DelphesTree::SkipTopologies(TString sel)
{
    return sel == "0B0S" || sel == "0B1S" || sel ==  "1B0S";
}
// __________________________________________________________________

void DelphesTree::InitCutVals(TString LjetTypeStr)
{
    // STEERING!!!

    // info/dump
    m_debug = 0;
    m_dumpASCII = false;
    m_doPdfAnalysis = false;
    
    // Large-R jet cuts:
    m_Ljetacut = 2.0;
    m_Ljptcut = 80.;
    // small-R jet cuts:
    m_jetacut = 2.5;
    m_jptcut = 25.;

    // to veto an isolated lepton:
    m_leptonEtaCut = 2.5;
    m_leptonPtCut = 25;

    // allow mistags!
    // fr the rate, see MyConsts.h
    m_runFakeBtag = true;

    map<kLJetChoice, TString> LJetChoiceDict;
    LJetChoiceDict[kDefault]       = "Default";
    LJetChoiceDict[kSoftDroppedP4] = "SoftDropped";
    LJetChoiceDict[kPrunedP4]      = "Pruned";
    LJetChoiceDict[kTrimmedP4]     = "Trimmed";
    m_LJetChoice = kTrimmedP4; // default LJet collection!
    bool changed = false;
    for (auto item : LJetChoiceDict) { // note: here item is std::pair<key,val> ;-)
        auto key = item.first;
        auto val = LJetChoiceDict[key];
        if (LjetTypeStr == val) {
            cout << "OK, using the user-defined " << val.Data() << " LJet collection!" << endl;
            m_LJetChoice = key;
            changed = true;
            break;
        }
    }
    if (changed)
        cout << "Using the on-demand LJet collection " << LJetChoiceDict[m_LJetChoice] << endl;
    else
        cout << "Using the default LJet collection " << LJetChoiceDict[m_LJetChoice] << endl;


    // top and W tagging study
    m_drtcut = 0.15; // it's for R=1 jets, but let's try to be more precise;)
    m_drWcut = 0.12; // it's for R=1 jets, but let's try to be more precise;)

    // for adding to W candidate a bjet candidate:
    m_mindrbjet = 1.;

    // JES steering!
    // if false, jes correction histos are filled for signal sample
    // if true, closure plots are fill instead unless forbidden on demand below
    m_CorrectLJES = true;
    m_CorrectSJES = true;
    m_FillJESclosure = false; // !!!
    m_xsectWeight = -1;
    map<TString, pair<double,double> > xsectMap;
    // sample: (xsect, Nevens)
    // only the three baseline ttbar samples supported:
    TString refSample = "pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200";
    // numbers from running e.g.
    // ./python/XsectStack.py lists/list_y0_1000GeV.txt | grep Weigh
    // 30.5.2021:
    // merged cross section from running cd  /home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/
    // ./GetXsectMergedAver.sh 
    // merged NLO+PS = mg5+Pythia8 cross sections and generated matched/weighted number of events:
    xsectMap[refSample]                                      =  make_pair<double,double>(137.756388, 772427.5);  // 437089.5); // unmerged: 271.090909
    xsectMap["pp_2tj_allhad_NLO_ptj1j2min200"]               =  make_pair<double,double>(4.606223,   920929.5); // 507675.0); // unmerged: 7.2811
    xsectMap["pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200"] =  make_pair<double,double>(19.484641,  948314.5 );  // 523247.0); // unmerged: 29.9365
    // 13.7.2021, updated 7.9.2021:
    xsectMap["pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200"]   =  make_pair<double,double>(21773.217852, 11499884.0);
    xsectMap["pp_2b2j_LO_matched_ptj1j2min200"]               =  make_pair<double,double>(604.47387,  2275190);
    xsectMap["pp_2b2j_LO_matched_ptj1min200_ptj2min60max200"] =  make_pair<double,double>(2659.769536, 2482645.5);

    if (m_ApplyXsectWeights || !m_CorrectLJES || !m_CorrectSJES || m_FillJESclosure) {
        // weight the pT spectra of jets by the sample xsection weights
        // this is either on demand, or compulsury when deriving JES, one needs smooth spectra
        // JK 26.10.2020
        double refLum = xsectMap[refSample].second / xsectMap[refSample].first; // lumi of the fist sample
        if (xsectMap.find(m_sampleTag) != xsectMap.end()) {
            auto xsectData = xsectMap[m_sampleTag];
            m_xsectWeight = refLum * xsectData.first / xsectData.second;
            cout << "KIND NOTIFICATION: Weighting the sample " << m_sampleTag.Data() << " with weight L*xs/N = "
                 << refLum << " * " << xsectData.first << " / " << xsectData.second << " = " << m_xsectWeight << endl;
        } else
            cerr << "ERROR getting the on-JES demand xsection information for sample tag " << m_sampleTag.Data() << " !!!" << endl;
    }

    // lower the pT cuts for JES derivation:
    if (!m_CorrectLJES)
        m_Ljptcut = 60.;
    if (!m_CorrectSJES)
        m_jptcut = 20.;
}

// __________________________________________________________________
void DelphesTree::MakeHistoMakers(int nReplicas)

{


    gDirectory -> mkdir("pdfinfo");
    gDirectory -> cd("pdfinfo");
    // gg/gq/qq production fractions
    m_InitialPartonsFracHisto = new TH1D("InitPartonsFracs", "InitPartonsFracs", 3, -0.5, 2.5);
    // next are meant in TeV, so 10 GeV bins
    //  a1 = [0.01 * x for x in range(0,21)]

    double shatBins[] = {0.0, 0.005,
                         0.01, 0.015, 0.02, 0.025,
                         0.03, 0.04, 0.05,
                         0.06, 0.07, 0.08, 0.09, 0.1,
                         0.11, 0.12, 0.13, 0.14, 0.15,
                         0.16, 0.17, 0.18, 0.19, 0.2,
                         0.22, 0.24, 0.26, 0.28, 0.30,
                         0.32, 0.35, 0.4,
                         0.45, 0.5, 0.55, 0.6, 0.7,
                         0.8, 1., 2., 3., 4., 7., 10., 14.};
    int nshat = sizeof(shatBins) / sizeof(double);
    m_InitialPartonsShat_gg = new TH1D("InitPartonsFracsShat_gg", "InitPartonsFracsShat_gg;#hat{s} [TeV]", nshat-1, (double*)shatBins);
    m_InitialPartonsShat_qg = new TH1D("InitPartonsFracsShat_qg", "InitPartonsFracsShat_qg;#hat{s} [TeV]", nshat-1, (double*)shatBins);
    m_InitialPartonsShat_qq = new TH1D("InitPartonsFracsShat_qq", "InitPartonsFracsShat_qq;#hat{s} [TeV]", nshat-1, (double*)shatBins);
    int pdgids[] = {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21};
    double xbins[] = {0., 0.005, 0.01,
                      0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5,
                      0.55, 0.60, 0.65, 0.7, 0.8, 1.};
    int nx = sizeof(xbins) / sizeof(double);
    for (auto pdgid : pdgids) {
        TString xname = "h_xPdfHisto_";
        xname += pdgid;
        m_xPdfHistos[pdgid] = new TH1D(xname, xname + ";x", nx-1, (double*)xbins);
    }
    gDirectory -> cd("../");

    // store values of basic cuts
    m_CutValsHisto = new TH1D("CutVals", "CutVals", 13, 0, 13);

    int nSels = sizeof(kSelections) / sizeof(TString);
    m_SelMigra = new TH2D("SelectionMigrations", "SelectionMigrations", nSels, 0, nSels, nSels, 0, nSels);
    m_hSelDet = new TH1D("DetectorPassedSelection", "DetectorPassedSelection", nSels, 0, nSels);
    m_hSelPtcl = new TH1D("ParticlePassedSelection", "ParticlePassedSelection", nSels, 0, nSels);

    // cutflow histos
    TString cutnames[] = {"Unweighted", "Initial", "nLjets", "ElVeto", "MuVeto", "LeptVeto", "nBtags", "AnySel"};
    int ncuts = sizeof(cutnames) / sizeof(TString);
    m_CutFlowHistoDet = new TH1D("CutFlowDet", "CutFlowDet", ncuts, 0, ncuts);
    m_CutFlowHistoPtcl = new TH1D("CutFlowPtcl", "CutFlowPtcl", ncuts, 0, ncuts);
    int ii = 1;
    for (auto cutname : cutnames) {
        m_CutFlowHistoDet -> GetXaxis() -> SetBinLabel(ii, cutname);
        m_CutFlowHistoPtcl -> GetXaxis() -> SetBinLabel(ii, cutname);
        ii++;
    }

    // inclusive in passying any of the selections
    m_AnyPassedTag = "AnySel";

    // not passing any
    m_NothingPassedTag = this->MakePassedTag(0,0,0);

    // passed ptcl, binned in ptcl
    m_hmaker_ptcl = new HistoMaker("Particle");
    // passed det, binned in det
    m_hmaker_det = new HistoMaker("Detector", nReplicas); // replicas option!!!
    // passed reco cuts, binned in parton
    /// m_hmaker_det_parton = new HistoMaker("Detector_parton");
    // passed det && ptcl, binned in ptcl
    m_hmaker_detptcl_ptcl = new HistoMaker("DetectorAndParticle_ptcl");
    // passed det && ptcl, binned in det
    m_hmaker_detptcl_det = new HistoMaker("DetectorAndParticle_det");

    HistoMaker* hmakers[] = {m_hmaker_ptcl,
                             m_hmaker_det,
                             //                     m_hmaker_det_parton,
                             //m_hmaker_detptcl_ptcl,
                             //m_hmaker_detptcl_det
                            };

    for (auto hmaker : hmakers) {
        cout << "Working on histo maker " << hmaker -> GetLevel() << endl;
        cout << "  ...adding cut level NoCuts" << endl;
        bool jetsOnly = true;
        hmaker -> AddCutLevel("NoCuts", jetsOnly);
        hmaker -> AddCutLevel("btagCuts", jetsOnly);
        jetsOnly = false;
        for (auto selection : kSelections) {
            cout << "  ...adding cut level " << selection.Data() << endl;
            hmaker -> AddCutLevel(selection, jetsOnly);
        }
        cout << "  ...adding cut level AnySel" << endl;
        hmaker -> AddCutLevel(m_AnyPassedTag, jetsOnly);
        bool MakeReplicas = nReplicas && (hmaker == m_hmaker_det); // make replicas only for the detector level!
        hmaker -> MakeAllJetHistos(MakeReplicas);
    }

    // this is suitable for migrations in diagonal topologies in det and parti levels!
    // Migrations particle-detector!
    m_hmaker_migra_ptcl = new HistoMaker("Migrations_ptcl_det");

    bool jetsOnly = false;
    for (auto selection : kSelections) {
      m_hmaker_migra_ptcl -> AddCutLevel(selection, jetsOnly);
      m_hmaker_migra_ptcl -> MakeAlljetMigrations(selection, "Particle", "Detector");
    }
    m_hmaker_migra_ptcl -> AddCutLevel(m_AnyPassedTag, jetsOnly);
    m_hmaker_migra_ptcl -> MakeAlljetMigrations(m_AnyPassedTag, "Particle", "Detector");
    if (m_isSignal) // && ( m_FillJESclosure || !m_CorrectLJES || !m_CorrectSJES) )
        m_hmaker_migra_ptcl -> MakeJESHistos(m_AnyPassedTag, m_FillJESclosure, !m_CorrectLJES || !m_CorrectSJES);

    // tagging histos for efficiency:
    m_hmaker_ptcl -> MakeTaggingHistos();
    m_hmaker_det -> MakeTaggingHistos();

    // parton, Apr2019
    m_hmaker_parton = new HistoMaker("Parton");
    // so far fill only no selection parton stuff
    jetsOnly = true;
    m_hmaker_parton -> AddCutLevel("NoCuts", jetsOnly);
    m_hmaker_parton -> MakePartonHistos();

    /*
    // migrations parton-detector
    m_hmaker_migra_parton = new HistoMaker("Migrations_parton_det");
    m_hmaker_migra_parton -> AddCutLevel(m_AnyPassedTag);
    m_hmaker_migra_parton -> MakeMigrations(m_AnyPassedTag, "Parton", "Detector");

    // migrations parton-particle
    m_hmaker_migra_partonptcl = new HistoMaker("Migrations_parton_ptcl");
    m_hmaker_migra_partonptcl -> AddCutLevel(m_AnyPassedTag);
    m_hmaker_migra_partonptcl -> MakeMigrations(m_AnyPassedTag, "Parton", "Particle");
    */

}


// __________________________________________________________________

void DelphesTree::ResetSelectionDictionaries()
{

    TString levels[] = {"Particle", "Detector"};
    // B=boosted, S=Semiboosted
    // the following does not mean resolved tops are neither are not found!
    for (auto level : levels) {
        for (auto selection : kSelections) {
            m_SelectDict [selection + "_" + level] = false;
        }
    }
}

// __________________________________________________________________

TString DelphesTree::MakePassedTag(int nbTags, int nWTags, int ntopTags)
{
    if (ntopTags > kntmax)
        ntopTags = kntmax;
    if (nbTags > knbmax)
        nbTags = knbmax;
    if (nWTags > knwmax)
        nWTags = knwmax;

    // capping the total!
    while (ntopTags + nWTags > kntmax ) {
        if (nWTags > 0) nWTags--;
        else ntopTags--;
    }



    // return a string like "2B1S"
    return Form("%iB%iS", ntopTags, nWTags);
}




// __________________________________________________________________
void DelphesTree::FillTagEffHistos(HistoMaker *m_hmaker, vector<TMyLorentzVector*> &detLjets,
                                   vector<TMyLorentzVector*>& partonTops, vector<TMyLorentzVector*> &partonWs)
{

    for (auto ljet : detLjets) {

        // fill histograms needed for tagging efficiency, see python/PlotTagEffs.py ;-)
        mindrdetadphi mindrt;
        mindrdetadphi mindrw;
        double pt = ljet->Pt();
        bool dumm1 = MatchedWithinDR(ljet, partonTops, mindrt, m_drtcut);
        bool dumm2 = MatchedWithinDR(ljet, partonWs, mindrw, m_drWcut);
        m_hmaker ->FillSingleTagDRHistos(mindrt, mindrw, pt, m_EventWeight);
        if (ljet -> Toptag()) {
            m_hmaker -> FillSingleTagKinemHistos("TopTagLjets", ljet, m_EventWeight);
            if (MatchedWithinDR(ljet, partonTops, mindrt, m_drtcut)) {
                m_hmaker -> FillSingleTagKinemHistos("TopTagttruthMatchedLjets", ljet, m_EventWeight);
            } else {
                m_hmaker -> FillSingleTagKinemHistos("TopTagnontMatchedLjets", ljet, m_EventWeight);
            }
        } else {
            if (MatchedWithinDR(ljet, partonTops, mindrt, m_drtcut)) { // FIXED TO mindrt!! 31.7.2020
                m_hmaker -> FillSingleTagKinemHistos("notTopTagttruthMatchedLjets", ljet, m_EventWeight);
            } else {
                m_hmaker -> FillSingleTagKinemHistos("notTopTagnontMatchedLjets", ljet, m_EventWeight);
            }
        }
        if (ljet -> Wtag()) {
            m_hmaker -> FillSingleTagKinemHistos("WTagLjets", ljet, m_EventWeight);
            if (MatchedWithinDR(ljet, partonWs, mindrw, m_drWcut)) {
                m_hmaker -> FillSingleTagKinemHistos("WTagWtruthMatchedLjets", ljet, m_EventWeight);
            } else {
                m_hmaker -> FillSingleTagKinemHistos("WTagnonWMatchedLjets", ljet, m_EventWeight);
            }
        } else {
            if (MatchedWithinDR(ljet, partonWs, mindrw, m_drWcut)) {
                m_hmaker -> FillSingleTagKinemHistos("notWTagWtruthMatchedLjets", ljet, m_EventWeight);
            } else {
                m_hmaker -> FillSingleTagKinemHistos("notWTagnonWMatchedLjets", ljet, m_EventWeight);
            }
        }

    }

}

// __________________________________________________________________


void DelphesTree::DumpAscii(HistoMaker *m_hmaker,
                            vector<TMyLorentzVector*>& detjets, vector<TMyLorentzVector*> &detLjets,
                            vector<TMyLorentzVector*>& partonTops, vector<TMyLorentzVector*> &partonWs)
{

    int ijet = -1;
    mindrdetadphi mindrt;
    mindrdetadphi mindrw;

    (*m_asciifile) << detLjets.size() << "|" << detjets.size() << "|";

    // Feb2019
    // detLjet loop to dump ASCII info for mathematicians;)
    // and fill also dR information between truth parton tops/Ws and detector LJets
    for (auto ljet : detLjets) {
        ijet++;
        // fill some important stuff; and dump ascii info for input to ML learning guys at KMA:)
        TString match = "";
        if (MatchedWithinDR(ljet, partonTops, mindrt, m_drtcut)) {
            match = "t";
            m_hmaker -> FillSingleTagKinemHistos("ttruthMatchedLjets", ljet, m_EventWeight);
        }
        if (MatchedWithinDR(ljet, partonWs, mindrw, m_drWcut)) {
            match = match + "W";
            m_hmaker -> FillSingleTagKinemHistos("WtruthMatchedLjets", ljet, m_EventWeight);
        }
        if (match == "")
            match = "light";

        (*m_asciifile) << ijet
                       << " " << mindrw.mindr
                       << " " << mindrt.mindr
                       << " " << ljet->Pt()
                       << " " << ljet->Eta()
                       << " " << ljet->Phi()
                       << " " << ljet->Tau32()
                       << " " << ljet->Tau21()
                       << " " << ljet->M() << "; ";
    } // Ljets
    (*m_asciifile) << endl;

}

// __________________________________________________________________
// Loop! ;-)

void DelphesTree::LoopBoosted(TString outtag, TString LjetTypeStr, int nReplicas)
{

    this -> InitCutVals(LjetTypeStr);

    if (nReplicas > 0 && (!m_CorrectLJES || !m_CorrectSJES || m_FillJESclosure) ) {
        cout << "KINDLY WARNING! Asked to derive L/S JES or fill closure: ("
             << !m_CorrectLJES << "," << !m_CorrectSJES << "," << m_FillJESclosure
             << ") while asked to fill " << nReplicas
             << " replicas. This ain't good nor has sense, switching to NOT filling the replicas!" << endl;
        nReplicas = 0;
    }

    if (nReplicas > 0) {
        cout << "OK, asked to fill " << nReplicas << " replicas! Hope you're ready for the ride!:)" << endl;
    } else {
        cout << "OK, not filling the replicas!" << endl;
    }


    if (fChain == 0) return;
    Long64_t nentries = m_treeReader->GetEntries();

    if (m_ApplyXsectWeights)
        outtag += "_weighted";

    TString outfilename = "analyzed_histos_" + outtag + ".root";
    TFile *outfile = new TFile(outfilename, "recreate");
    this -> MakeHistoMakers(nReplicas);
    double *repweights = 0;
    if (nReplicas > 0) {
        repweights = new double[nReplicas];
    }

    int verbose = 10000;
    int ToRunOver = nentries;

    Long64_t nbytes = 0, nb = 0;

    cout << "Some settings" << endl
         << " debug     : " << m_debug << endl
         << " dumpASCII : " << m_dumpASCII << endl
         << " jetacut   : " << m_jetacut << endl
         << " Ljetacut  : " << m_Ljetacut << endl
         << " Ljptcut   : " << m_Ljptcut << endl
         << " jptcut    : " << m_jptcut << endl
         << " drtcut    : " << m_drtcut << endl
         << " drWcut    : " << m_drWcut << endl
         << " mindrbjet : " << m_mindrbjet << endl
         << " CorrectLJES        : " << m_CorrectLJES << endl
         << " CorrectSJES        : " << m_CorrectSJES << endl
         << " FillJESclosure     : " << m_FillJESclosure << endl
         << " nReplicas          : " << nReplicas << endl
         << " ApplyXsectWeights  : " << m_ApplyXsectWeights << endl;
         //<< " " << endl

    cout << "Will loop over " << ToRunOver << " events." << endl;

    m_asciifile = 0;
    if (m_dumpASCII) {
        TString asciiname = outfilename;
        asciiname.ReplaceAll("analyzed_histos", "ascii");
        asciiname.ReplaceAll(".root", ".txt");
        m_asciifile = new ofstream(asciiname.Data());
    }


    // loop over events in the TTree
    for (Long64_t jentry = 0; jentry < ToRunOver; jentry++) {

        bool passedDetector = false;
        bool passedParticle = false;

        this -> ResetSelectionDictionaries();
        m_EventWeight = 1.;
        if (m_xsectWeight > 0)
            m_EventWeight *= m_xsectWeight;

        if (jentry % verbose == 0) {
            cout << "Processing " << jentry << "/" << nentries << " [" << Form("%2.1f", jentry/(1.*nentries) * 100) << "%]" << endl;
        }
        if (m_debug > 0) cout << "Event " << jentry << endl;

        m_treeReader->ReadEntry(jentry);

        int icutdet = 0;
        int icutptcl = 0;
        m_CutFlowHistoDet -> Fill(icutdet++);
        m_CutFlowHistoPtcl -> Fill(icutptcl++);

        // skip zero-weight events from mg5
        // which are zero due to MLM matching
        HepMCEvent *event = (HepMCEvent*)m_branchEvent->At(0);
        double Event_Weight = event->Weight;
        if (Event_Weight < 1.e-20)
            continue;


	if (m_doPdfAnalysis) {
	// access collided partons flavours:
	int pid1 = event->ID1;
	int pid2 = event->ID2;
	float px1 = event->X1;
	float px2 = event->X2;
    double shat = 14. * px1 * px2; // hardcoded the 14 TeV here!
    /*
    cout << " Weight=" << Event_Weight
	 << " pid1=" << pid1
	 << " pid2=" << pid2
	 << " px1=" << px1
	 << " px2=" << px2
	 << endl;
    */
    if (pid1 == 21 && pid2 == 21) {
        // gg-initiated
        m_InitialPartonsFracHisto -> Fill(0., Event_Weight);
        m_InitialPartonsShat_gg -> Fill(shat, Event_Weight);
    } else if ( (pid1 == 21 || pid2 == 21) && (abs(pid1) < 6 || abs(pid2) < 6 )) {
        // qg-initiated
        m_InitialPartonsFracHisto -> Fill(1., Event_Weight);
        m_InitialPartonsShat_qg -> Fill(shat, Event_Weight);
    } else if ( abs(pid1) < 6 && abs(pid2) < 6 ) {
        // qq initiated
        m_InitialPartonsFracHisto -> Fill(2., Event_Weight);
        m_InitialPartonsShat_qq -> Fill(shat, Event_Weight);
    } else {
        cerr << "Unclear initial parton fractions! This should not happen!" << endl;
        //m_InitialPartonsFracHisto -> Fill(-1., Event_Weight);
    }
    if (pid1 != 0) {
        m_xPdfHistos[pid1]->Fill(px1, Event_Weight);
    }
    if (pid2 != 0) {
        m_xPdfHistos[pid2]->Fill(px2, Event_Weight);
    }

        // print some metadata:
        /*
        if (m_debug > 1) {
            cout << " Event_CrossSection: " << Event_CrossSection[0] << endl
                 << " Event_Weight: " << Event_Weight[0] << endl
                 << " Event_X1: " << Event_X1[0] << endl
                 << " Event_X2: " << Event_X2[0] << endl
                 << " Event_MPI: " << Event_MPI[0] << endl
                 << " Event_ID1: " << Event_ID1[0] << endl
                 << " Event_ID2: " << Event_ID2[0] << endl
                 << " Event_AlphaQED: " << Event_AlphaQED[0] << endl
                 << " Event_AlphaQCD: " << Event_AlphaQCD[0] << endl;
        }
        */

	}
        /* JK 9.12.2018
         * Need to keep special regard to the possible definition/overlap of topologies
         * So let's count first boosted tops
         * which should be exclusive to semiboosted
         * and one should treat properly overlap between jets!
         * */

        // +------------------------+
        // |  Parton level cuts     |
        // +------------------------+
        // none for the moment;-)
        // just make the parton tops and histograms of their kinematics

        if (m_debug > 1) cout << "  making parton tops" << endl;
        vector<TMyLorentzVector*> partonTops = MakeGenTops(0., 999., 52);
        if (m_debug > 1) cout << "     ...got " << partonTops.size() << " parton tops" << endl;
        m_hmaker_parton -> FillJetHistos("NoCuts", "Top", partonTops, m_EventWeight, kallJets, 0);

        if (m_debug > 1) cout << "  making parton Ws" << endl;
        vector<TMyLorentzVector*> partonWs = MakeGenWs(0., 999., 52);
        //if (m_debug > 1) cout << "     ...got " << partonWs.size() << " parton Ws" << endl;
        m_hmaker_parton -> FillJetHistos("NoCuts", "W", partonWs, m_EventWeight, kallJets, 0);

        // warning, so far only the first pair is returned!!
        if (m_debug > 1) cout << "  ok, making parton ditops" << endl;
        vector<TMyLorentzVector*> partonDiTops = MakeGenDiTops(partonTops);
        m_hmaker_parton -> FillJetHistos("NoCuts", "DiTop", partonDiTops, m_EventWeight, kallJets, 0);


        // +------------------------+
        // |  Detector level cuts:  |
        // +------------------------+

        // JK Jan 2019:
        // just make the Ljets and
        // then count number of boosted, semiboosted,
        // and fill the migration matrix between regimes

        m_CutFlowHistoDet -> Fill(icutdet++);

        // no pT cut, detector jets:
        /*
        vector<TMyLorentzVector*>  alldetLjets = MakeJets(1, 1, 0., etacut, CorrectJES);
        vector<TMyLorentzVector*>  alldetjets  = MakeJets(1, 0, 0., etacut, CorrectJES);
        if (m_debug > 1) cout << "  got " << alldetLjets.size() << " jets." << endl;
        m_hmaker_det -> FillJetHistos("NoCuts", "LJet", alldetLjets, m_EventWeight, kallJets, 1);
        m_hmaker_det -> FillJetHistos("NoCuts", "Jet", alldetjets, m_EventWeight, kallJets, 0);
        */
        // pt cuts on all jets, get the collection of those passing:
        vector<TMyLorentzVector*> detLjets    = MakeJets(1, 1, m_Ljptcut, m_Ljetacut, m_CorrectLJES,
                                                         m_runFakeBtag, m_LJetChoice);
        vector<TMyLorentzVector*> detjets     = MakeJets(1, 0, m_jptcut, m_jetacut, m_CorrectSJES, m_runFakeBtag);
        if (m_debug > 1) cout << "  got " << detLjets.size() << " large det jets in pt eta acceptance." << endl;

        // make leptons and veto on them:
        bool hasElDet = this -> MakeLepton(kEjets, true, detjets);
        bool hasMuDet = this -> MakeLepton(kMujets, true, detjets);

        if (m_dumpASCII && detLjets.size() > 0)
            this -> DumpAscii(m_hmaker_det, detjets, detLjets, partonTops, partonWs);

        this -> FillTagEffHistos(m_hmaker_det, detLjets, partonTops, partonWs);

        // vector<TMyLorentzVector*> detDijets;

        // fill small jet histograms before cuts (except jet pT and eta cut
        // we always look at jets in pT and eta acceptance;-)
        if (detjets.size() > 0) {
            m_hmaker_det -> FillJetHistos("NoCuts", "Jet", detjets, m_EventWeight, kallJets, 0);
            m_hmaker_det -> FillJetHistos("NoCuts", "bJet", detjets, m_EventWeight, kbJets, 0);
            m_hmaker_det -> FillJetHistos("NoCuts", "nonbJet", detjets, m_EventWeight, knonbJets, 0);
        }


        vector<TMyLorentzVector*> dettops = MakeTops(detjets, detLjets, m_mindrbjet);
        vector<TMyLorentzVector*> detDiTops = MakeDiTops(dettops);
        if (dettops.size() > 1 && detDiTops.size() < 1)
            cerr << "WARNING: det: got >=1 tops but zero ditops!" << endl;
        vector<TMyLorentzVector*> detFourTops = MakeFourTops(dettops);
        // store also in the pseudotop struct to compute cos theta* etc:
        kPseudotop pst_detDiTops;
        m_MadePstDetDiTops = false;
        if (dettops.size() > 1) {
            // cout << " Making Detector Pout..." << endl;
            pst_detDiTops.pseudotoplepton = TMyLorentzVector(*dettops[0]);
            pst_detDiTops.pseudotophadron = TMyLorentzVector(*dettops[1]);
            pst_detDiTops.pseudottbar = pst_detDiTops.pseudotophadron + pst_detDiTops.pseudotoplepton;
            this -> MakePout(pst_detDiTops);
            m_MadePstDetDiTops = true;
        }
        this -> MakeGlobalVars(false, m_globalVars_det, detjets, detLjets);


        if (nReplicas > 0) {
            for (unsigned int irep = 0; irep < nReplicas; ++irep) {
                repweights[irep] = m_rand->Poisson(1.);
            }
            m_hmaker_det -> SetReplicaWeights(repweights);
        }

        if (detjets.size() > 0) {
            m_hmaker_det -> FillGlobalHistos("NoCuts", m_globalVars_det, m_EventWeight);
        }

        if (detLjets.size() > 0) {

            m_CutFlowHistoDet -> Fill(icutdet++);


            if (!hasElDet) {
                m_CutFlowHistoDet -> Fill(icutdet++);
            } else icutdet++;
            if (!hasMuDet) {
                m_CutFlowHistoDet -> Fill(icutdet++);
            } else icutdet++;

            if (hasElDet || hasMuDet) {

            } else {

                m_CutFlowHistoDet -> Fill(icutdet++);

                if (m_debug > 1)
                    cout << "Reconstructed" <<
                            " dettops: " << dettops.size() <<
                            " detDiTops: " << detDiTops.size() <<
                            " detFourTops: " << detFourTops.size() << endl;

                if (m_debug > 1) cout << "  first jet pt: " << detLjets[0]->Pt() << endl;
                if (m_debug > 1) cout << "  ok, passed njets cut" << endl;

                m_hmaker_det -> FillJetHistos("NoCuts", "LJet", detLjets, m_EventWeight, kallJets, 1);

                // count how many tags found
                m_ndetbTags = CountBtags(detjets);
                m_ndetWTags = CountWtags(detLjets);
                m_ndettopTags = Counttoptags(detLjets);

                // 31.3.2020, 24.4.2020
                // number of top and W tags is not sufficient
                // it may not be possible to find a small b-jet to add it to W jet
                // and this should be counted as inefficiency...and event should me marked
                // as NOT passign 1B1S etc.
                // This is ensured below by requiring that bot tops were reconsructed and ttbar related
                // special quantities are reconstructed in pseudotop structure

                // created a string tag encoding which topology was passed
                m_detPassedTag = this->MakePassedTag(m_ndetbTags, m_ndetWTags, m_ndettopTags);
                if (m_debug > 0) {
                    cout << "det passed tag: " << m_detPassedTag.Data()
                         << " m_ndetbTags: " << m_ndetbTags
                         << " m_ndetWTags: " << m_ndetWTags
                         << " m_ndettopTags: " << m_ndettopTags << endl;
                }

                // require >= 2b tags!
                if (m_ndetbTags < 2) {
                    // important! ;-)
                    m_detPassedTag = m_NothingPassedTag;
                } else {
                    m_CutFlowHistoDet -> Fill(icutdet++);

                    if (detjets.size() > 0) {
                        m_hmaker_det -> FillJetHistos("btagCuts", "Jet", detjets, m_EventWeight, kallJets, 0);
                        m_hmaker_det -> FillJetHistos("btagCuts", "bJet", detjets, m_EventWeight, kbJets, 0);
                        m_hmaker_det -> FillJetHistos("btagCuts", "nonbJet", detjets, m_EventWeight, knonbJets, 0);
                    }
                }

                // this is effectivelly also after 2b req.;-)
                if (m_detPassedTag != m_NothingPassedTag && ! this->SkipTopologies(m_detPassedTag)) {
                    m_CutFlowHistoDet -> Fill(icutdet++);

                    passedDetector = true;
                    m_hmaker_det -> FillJetHistos(m_detPassedTag, "LJet", detLjets, m_EventWeight, kallJets, 1);
                    m_hmaker_det -> FillJetHistos(m_detPassedTag, "Jet", detjets, m_EventWeight, kallJets, 0);
                    m_hmaker_det -> FillJetHistos(m_detPassedTag, "bJet", detjets, m_EventWeight, kbJets, 0);
                    m_hmaker_det -> FillJetHistos(m_detPassedTag, "nonbJet", detjets, m_EventWeight, knonbJets, 0);

                    // AnySel dir. fill:
                    m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "LJet", detLjets, m_EventWeight, kallJets, 1);
                    m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "Jet", detjets, m_EventWeight, kallJets, 0);
                    m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "bJet", detjets, m_EventWeight, kbJets, 0);
                    m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "nonbJet", detjets, m_EventWeight, knonbJets, 0);
                    m_hmaker_det -> FillGlobalHistos(m_AnyPassedTag, m_globalVars_det, m_EventWeight);

                    //            m_hmaker_det_parton -> FillJetHistos(m_AnyPassedTag, "Dijet", partonDijets, m_EventWeight, kallJets, 0);

                    // fill W histos, no need to check if (m_ndetWTags > 0)
                    m_hmaker_det -> FillJetHistos(m_detPassedTag, "W", detLjets, m_EventWeight, kWTagJets, 1);
                    // fill top histos, no need to check if (m_ndettopTags > 0)
                    m_hmaker_det -> FillJetHistos(m_detPassedTag, "Top", dettops, m_EventWeight, kTopTagJets, 1);
                    // fill ditop histos JK May2019, check >=2 reco tops via pst struct 2020
                    if (m_MadePstDetDiTops) {
                        m_hmaker_det -> FillJetHistos(m_detPassedTag, "DiTop", detDiTops, m_EventWeight, kallJets, 0);
                        m_hmaker_det -> FillTTbarSpecialHistos(m_detPassedTag, "DiTop", pst_detDiTops, m_globalVars_det, m_EventWeight, nReplicas);
                    }
                    m_hmaker_det -> FillGlobalHistos(m_detPassedTag, m_globalVars_det, m_EventWeight);

                    // fill fourtophistos, no need to check if (m_ndettopTags > 3)
                    m_hmaker_det -> FillJetHistos(m_detPassedTag, "FourTop", detFourTops, m_EventWeight, kallJets, 0);

                    // bloody the same, but now fill AnySel directory, too;)
                    // fill W histos, no need to check if (m_ndetWTags > 0)
                    m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "W", detLjets, m_EventWeight, kWTagJets, 1);
                    // fill top histos, no need to check if (m_ndettopTags > 0)
                    m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "Top", dettops, m_EventWeight, kTopTagJets, 1);
                    // fill ditop histos JK May2019, check >=2 reco tops via pst struct 2020
                    if (m_MadePstDetDiTops) {
                        m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "DiTop", detDiTops, m_EventWeight, kallJets, 0);
                        m_hmaker_det -> FillTTbarSpecialHistos(m_AnyPassedTag, "DiTop", pst_detDiTops, m_globalVars_det, m_EventWeight, nReplicas);
                    }
                    // fill fourtophistos, no need to check if (m_ndettopTags > 3)
                    m_hmaker_det -> FillJetHistos(m_AnyPassedTag, "FourTop", detFourTops, m_EventWeight, kallJets, 0);

                } // pased one of the selections

            } // lepton veto

        } // Ljets size


        // +------------------------+
        // |  Particle level cuts:  |
        // +------------------------+

        bool ptclCorrectJES = false; // ALWAYS false!;)
        bool FillReplicas = false;
        m_CutFlowHistoPtcl -> Fill(icutptcl++);

        // no pT cut, ptcl jets:
/*
        vector<TMyLorentzVector*>  allptclLjets = MakeJets(0, 1, 0., etacut, ptclCorrectJES);
        vector<TMyLorentzVector*>  allptcljets  = MakeJets(0, 0, 0., etacut, ptclCorrectJES);
        if (m_debug > 1) cout << "  got ptcl " << allptclLjets.size() << " jets." << endl;
        m_hmaker_ptcl -> FillJetHistos("NoCuts", "LJet", allptclLjets, m_EventWeight, kallJets, 1);
        m_hmaker_ptcl -> FillJetHistos("NoCuts", "Jet", allptcljets, m_EventWeight, kallJets, 0);
*/

        // pt cuts on all jets, get the collection of those passing:
        vector<TMyLorentzVector*> ptclLjets    = MakeJets(0, 1, m_Ljptcut, m_Ljetacut, ptclCorrectJES,
                                                          m_runFakeBtag, m_LJetChoice);
        vector<TMyLorentzVector*> ptcljets     = MakeJets(0, 0, m_jptcut, m_jetacut, ptclCorrectJES, m_runFakeBtag);
        if (m_debug > 1) cout << "  got " << ptclLjets.size() << " large pctl jets in pt eta acceptance." << endl;

        // make leptons and veto on them:
        bool hasElPtcl = this -> MakeLepton(kEjets, false, ptcljets);
        bool hasMuPtcl = this -> MakeLepton(kMujets, false, ptcljets);

        this -> FillTagEffHistos(m_hmaker_ptcl, ptclLjets, partonTops, partonWs);

        vector<TMyLorentzVector*> ptcltops = MakeTops(ptcljets, ptclLjets, m_mindrbjet);
        vector<TMyLorentzVector*> ptclDiTops = MakeDiTops(ptcltops);
        if (ptcltops.size() > 1 && ptclDiTops.size() < 1)
            cerr << "WARNING: ptcl: got >=1 tops but zero ditops!" << endl;
        vector<TMyLorentzVector*> ptclFourTops = MakeFourTops(ptcltops);
        // store also in the pseudotop struct to compute cos theta* etc:
        kPseudotop pst_ptclDiTops;
        m_MadePstPtclDiTops = false;
        if (ptcltops.size() > 1) {
            //cout << " Making Particle Pout..." << endl;
            pst_ptclDiTops.pseudotoplepton = TMyLorentzVector(*ptcltops[0]);
            pst_ptclDiTops.pseudotophadron = TMyLorentzVector(*ptcltops[1]);
            pst_ptclDiTops.pseudottbar = pst_ptclDiTops.pseudotophadron + pst_ptclDiTops.pseudotoplepton;
            this -> MakePout(pst_ptclDiTops);
            m_MadePstPtclDiTops = true;
        }
        this -> MakeGlobalVars(true, m_globalVars_ptcl, ptcljets, ptclLjets);

        if (detjets.size() > 0) {
            m_hmaker_ptcl -> FillGlobalHistos("NoCuts", m_globalVars_ptcl, m_EventWeight);
        }


        if (ptcljets.size() > 0) {
            m_hmaker_ptcl -> FillJetHistos("NoCuts", "Jet",  ptcljets, m_EventWeight, kallJets, 0);
            m_hmaker_ptcl -> FillJetHistos("NoCuts", "bJet", ptcljets, m_EventWeight, kbJets, 0);
            m_hmaker_ptcl -> FillJetHistos("NoCuts", "nonbJet", ptcljets, m_EventWeight, knonbJets, 0);
        }

        if (ptclLjets.size() > 0) {

            m_CutFlowHistoPtcl -> Fill(icutptcl++);

            // TO FINISH, ALSO FOR PTCL LEVEL!!!

            if (!hasElPtcl) {
                m_CutFlowHistoPtcl -> Fill(icutptcl++);
            } else icutptcl++;

            if (!hasMuPtcl) {
                m_CutFlowHistoPtcl -> Fill(icutptcl++);
            } else icutptcl++;

            if (hasElPtcl || hasMuPtcl) {

            } else {
                m_CutFlowHistoPtcl -> Fill(icutptcl++);
                if (m_debug > 1)
                    cout << "Reconstructed" <<
                            " ptcltops: " << ptcltops.size() <<
                            " ptclDiTops: " << ptclDiTops.size() <<
                            " ptclFourTops: " << ptclFourTops.size() << endl;

                m_hmaker_ptcl -> FillJetHistos("NoCuts", "LJet", ptclLjets, m_EventWeight, kallJets, 1);

                m_nptclbTags = CountBtags(ptcljets);
                m_nptclWTags = CountWtags(ptclLjets);
                m_nptcltopTags = Counttoptags(ptclLjets);

                m_ptclPassedTag = this->MakePassedTag(m_nptclbTags, m_nptclWTags, m_nptcltopTags);
                if (m_debug > 0)
                    cout << "ptcl passed tag: " << m_ptclPassedTag.Data() << endl;

                // but do require >= 2b tags!
                if (m_nptclbTags < 2) {
                    m_ptclPassedTag = m_NothingPassedTag;
                } else {

                    m_CutFlowHistoPtcl -> Fill(icutptcl++);
                    if (ptcljets.size() > 0) {
                        m_hmaker_ptcl -> FillJetHistos("btagCuts", "Jet",  ptcljets, m_EventWeight, kallJets, 0);
                        m_hmaker_ptcl -> FillJetHistos("btagCuts", "bJet", ptcljets, m_EventWeight, kbJets, 0);
                        m_hmaker_ptcl -> FillJetHistos("btagCuts", "nonbJet", ptcljets, m_EventWeight, knonbJets, 0);
                    }
                }
                if (m_ptclPassedTag != m_NothingPassedTag && ! this->SkipTopologies(m_ptclPassedTag)) {
                    m_CutFlowHistoPtcl -> Fill(icutptcl++);
                    passedParticle = true;
                    m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "LJet", ptclLjets, m_EventWeight, kallJets, 1);
                    m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "Jet", ptcljets, m_EventWeight, kallJets, 0);
                    m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "bJet", ptcljets, m_EventWeight, kbJets, 0);
                    m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "nonbJet", ptcljets, m_EventWeight, knonbJets, 0);
                    m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "LJet", ptclLjets, m_EventWeight, kallJets, 1);
                    m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "Jet", ptcljets, m_EventWeight, kallJets, 0);
                    m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "bJet", ptcljets, m_EventWeight, kbJets, 0);
                    m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "nonbJet", ptcljets, m_EventWeight, knonbJets, 0);
                    m_hmaker_ptcl -> FillGlobalHistos(m_AnyPassedTag, m_globalVars_ptcl, m_EventWeight);


                    //m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "Dijet", ptclDijets, m_EventWeight, kallJets, 0);
                    // fill migrations: parton-particle:
                    // m_hmaker_migra_partonptcl -> FillMigrations(m_AnyPassedTag, partonTops, dummyjets, partonDijets, ptclLjets, ptcljets, ptclDijets, m_EventWeight);

                    // fill W histos, no need to check if (m_nptclWTags > 0)
                    m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "W", ptclLjets, m_EventWeight, kWTagJets, 1);
                    // fill top histos, no need to check if (m_nptcltopTags > 0)
                    m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "Top", ptcltops, m_EventWeight, kTopTagJets, 1);
                    // fill ditop histos JK May2019, check >=2 reco tops via pst struct 2020
                    if (m_MadePstPtclDiTops) {
                        m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "DiTop", ptclDiTops, m_EventWeight, kallJets, 0);
                        m_hmaker_ptcl -> FillTTbarSpecialHistos(m_ptclPassedTag, "DiTop", pst_ptclDiTops, m_globalVars_ptcl, m_EventWeight, FillReplicas);
                    }
                    m_hmaker_ptcl -> FillGlobalHistos(m_ptclPassedTag, m_globalVars_ptcl, m_EventWeight);
                    // fill fourtophistos, no need to check if (m_nptcltopTags > 3)
                    m_hmaker_ptcl -> FillJetHistos(m_ptclPassedTag, "FourTop", ptclFourTops, m_EventWeight, kallJets, 0);

                    // bloody the same, but now fill AnySel directory, too;)
                    // fill W histos, no need to check if (m_nptclWTags > 0)
                    m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "W", ptclLjets, m_EventWeight, kWTagJets, 1);
                    // fill top histos, no need to check if (m_nptcltopTags > 0)
                    m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "Top", ptcltops, m_EventWeight, kTopTagJets, 1);
                    // fill ditop histos JK May2019, check >=2 reco tops via pst struct 2020
                    if (m_MadePstPtclDiTops) {
                        m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "DiTop", ptclDiTops, m_EventWeight, kallJets, 0);
                        m_hmaker_ptcl -> FillTTbarSpecialHistos(m_AnyPassedTag, "DiTop", pst_ptclDiTops, m_globalVars_ptcl, m_EventWeight, FillReplicas);
                    }
                    // fill fourtophistos, no need to check if (m_nptcltopTags > 3)
                    m_hmaker_ptcl -> FillJetHistos(m_AnyPassedTag, "FourTop", ptclFourTops, m_EventWeight, kallJets, 0);


                } // passed one of the selections

            } // lepton veto

        } // ptclLjets size


        int i = 0;
        int iptcl = -1;
        int idet = -1;
        for (auto selection : kSelections) {
            if (passedParticle && m_ptclPassedTag == selection)
                iptcl = i;
            if (passedDetector && m_detPassedTag == selection)
                idet = i;
            i++;
        }
        if (passedDetector)
            m_hSelDet -> Fill(idet);
        if (passedParticle)
            m_hSelPtcl -> Fill(iptcl);

        if (m_isSignal && ( m_FillJESclosure || !m_CorrectLJES || !m_CorrectSJES) ) {
            // now we have both particle and detector jets and we can do some JES derivation or closure;)
            m_hmaker_migra_ptcl -> FillJESHistos(m_AnyPassedTag, ptclLjets, detLjets, m_CorrectLJES, true, m_EventWeight);
            m_hmaker_migra_ptcl -> FillJESHistos(m_AnyPassedTag, ptcljets, detjets, m_CorrectSJES, false, m_EventWeight);
        }


        // +-----------------------------+
        // | passed detector && particle |
        // +-----------------------------+

        // Migration matrices and eff and acc ingredients
        if (passedDetector && passedParticle) {

            if (idet >=0 && iptcl >= 0) {
                if (m_debug > 1)
                    cout << "filling sel. migra at bins " << iptcl << " " << idet << endl;
                m_SelMigra -> Fill(iptcl, idet);
            }

            // Passed both particle and detector levels:
            // fill migrations: particle-detector:
            // Note: Ljets contain Wjets, their tag is checked in filling;-)
            // so far only support for diagonality in selections at det and ptcl levels!
            if (m_ptclPassedTag == m_detPassedTag && m_MadePstDetDiTops && m_MadePstPtclDiTops) {

                // fill these only for events diagonal in topologies between particle and detector levels! ;-)

                TString selection = m_detPassedTag;
                m_hmaker_migra_ptcl -> FillAlljetMigrations(selection,
                                                            ptclLjets, ptcltops, ptclDiTops, ptclFourTops,
                                                            detLjets, dettops, detDiTops, detFourTops, m_EventWeight);
                // fill special migrations:
                m_hmaker_migra_ptcl -> FillSpecialMigrations(selection, "DiTop", pst_ptclDiTops, pst_detDiTops, m_EventWeight);
                // migrations of global variables
                m_hmaker_migra_ptcl -> FillGlobalMigrations(selection, "Global", m_globalVars_ptcl, m_globalVars_det, m_EventWeight);


                m_hmaker_migra_ptcl -> FillAlljetMigrations(m_AnyPassedTag,
                                                            ptclLjets, ptcltops, ptclDiTops, ptclFourTops,
                                                            detLjets, dettops, detDiTops, detFourTops, m_EventWeight);
                // special migrations!
                m_hmaker_migra_ptcl -> FillSpecialMigrations(m_AnyPassedTag, "DiTop", pst_ptclDiTops, pst_detDiTops, m_EventWeight);
                // migrations of global variables
                m_hmaker_migra_ptcl -> FillGlobalMigrations(m_AnyPassedTag, "Global", m_globalVars_ptcl, m_globalVars_det, m_EventWeight);

                // histograms for unfolfing corrections:
                // but swithcing off, and should not be needed for unfoding corrs anyway, either,
                // as corrs are nmade form histos massign ptc or det, and from projections of migrations!
                // jk 3.12.2020
                /*
                // fill 1D histos binned in Detector:
                m_hmaker_detptcl_det -> FillJetHistos(m_AnyPassedTag, "LJet", detLjets, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_det -> FillJetHistos(m_AnyPassedTag, "Jet", detjets, m_EventWeight, kallJets, 0);

                m_hmaker_detptcl_det -> FillJetHistos(m_AnyPassedTag, "W", detLjets, m_EventWeight, kWTagJets, 1);
                m_hmaker_detptcl_det -> FillJetHistos(m_AnyPassedTag, "Top", dettops, m_EventWeight, kTopTagJets, 1);
                m_hmaker_detptcl_det -> FillJetHistos(m_AnyPassedTag, "DiTop", detDiTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_det -> FillJetHistos(m_AnyPassedTag, "FourTop", detFourTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_det -> FillGlobalHistos(m_AnyPassedTag, m_globalVars_det, m_EventWeight);
                //m_hmaker_detptcl_det -> FillTTbarSpecialHistos(m_AnyPassedTag, "DiTop", pst_detDiTops, m_EventWeight);

                m_hmaker_detptcl_det -> FillJetHistos(m_detPassedTag, "W", detLjets, m_EventWeight, kWTagJets, 1);
                m_hmaker_detptcl_det -> FillJetHistos(m_detPassedTag, "Top", dettops, m_EventWeight, kTopTagJets, 1);
                m_hmaker_detptcl_det -> FillJetHistos(m_detPassedTag, "DiTop", detDiTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_det -> FillJetHistos(m_detPassedTag, "FourTop", detFourTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_det -> FillGlobalHistos(m_detPassedTag, m_globalVars_det, m_EventWeight);
                //m_hmaker_detptcl_det -> FillTTbarSpecialHistos(m_detPassedTag, "DiTop", pst_detDiTops, m_EventWeight);

                // fill 1D histos binned in Particle:
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_AnyPassedTag, "LJet", ptclLjets, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_AnyPassedTag, "Jet", ptcljets, m_EventWeight, kallJets, 0);

                m_hmaker_detptcl_ptcl -> FillJetHistos(m_AnyPassedTag, "W", ptclLjets, m_EventWeight, kWTagJets, 1);
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_AnyPassedTag, "Top", ptcltops, m_EventWeight, kTopTagJets, 1);
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_AnyPassedTag, "DiTop", ptclDiTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_AnyPassedTag, "FourTop", ptclFourTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_ptcl -> FillGlobalHistos(m_AnyPassedTag, m_globalVars_ptcl, m_EventWeight);
                //m_hmaker_detptcl_ptcl -> FillTTbarSpecialHistos(m_AnyPassedTag, "DiTop", pst_ptclDiTops, m_EventWeight);

                m_hmaker_detptcl_ptcl -> FillJetHistos(m_ptclPassedTag, "W", ptclLjets, m_EventWeight, kWTagJets, 1);
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_ptclPassedTag, "Top", ptcltops, m_EventWeight, kTopTagJets, 1);
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_ptclPassedTag, "DiTop", ptclDiTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_ptcl -> FillJetHistos(m_ptclPassedTag, "FourTop", ptclFourTops, m_EventWeight, kallJets, 0);
                m_hmaker_detptcl_ptcl -> FillGlobalHistos(m_ptclPassedTag, m_globalVars_ptcl, m_EventWeight);
                //m_hmaker_detptcl_ptcl -> FillTTbarSpecialHistos(m_ptclPassedTag, "DiTop", pst_ptclDiTops, m_EventWeight);
               */

            } // same selection


        } // det && ptcl


        if (m_debug > 1)
            cout << "Cleaning!" << endl;

        // clean up:
        for (auto jet : detjets)   { delete jet; }
        for (auto jet : detLjets)  { delete jet; }

        for (auto jet : dettops)      { if (jet) delete jet; }
        for (auto jet : detDiTops)    { if (jet) delete jet; }
        for (auto jet : detFourTops)  { if (jet) delete jet; }

        for (auto jet : ptcljets)  { delete jet; }
        for (auto jet : ptclLjets) { delete jet; }

        for (auto jet : ptcltops)      { if (jet) delete jet; }
        for (auto jet : ptclDiTops)    { if (jet) delete jet; }
        for (auto jet : ptclFourTops)  { if (jet) delete jet; }

        //for (auto jet : alldetjets) {
        //    delete jet;
        //}
        for (auto jet : partonTops) { delete jet; }
        for (auto jet : partonWs)   { delete jet; }
        for (auto jet : partonDiTops) { delete jet; }

    }
    // +--------------------+
    // |  End of the loop!  |
    // +--------------------+

    cout << "End of event loop!" << endl;
    if (m_asciifile) m_asciifile->close();

    cout << "Some histo labels postprocessing..." << endl;

    // postprocessing

    m_InitialPartonsFracHisto -> GetXaxis() -> SetBinLabel(1, "gg");
    m_InitialPartonsFracHisto -> GetXaxis() -> SetBinLabel(2, "qg");
    m_InitialPartonsFracHisto -> GetXaxis() -> SetBinLabel(3, "qq");


    int i = 1;
    for (auto selection : kSelections) {
        m_SelMigra -> GetXaxis() -> SetBinLabel(i, selection);
        m_SelMigra -> GetYaxis() -> SetBinLabel(i, selection);
        m_hSelDet -> GetXaxis() -> SetBinLabel(i, selection);
        m_hSelPtcl -> GetXaxis() -> SetBinLabel(i, selection);
        i++;
    }

    int ibin = 1;
    // BAD behaviour in hadd...:
    m_CutValsHisto -> Fill(ibin, m_drtcut);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "Match #DeltaR_{t} cut");
    m_CutValsHisto -> Fill(ibin, m_drWcut);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "Match #DeltaR_{W} cut");
    m_CutValsHisto -> Fill(ibin, m_jetacut);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "#jeta cut");
    m_CutValsHisto -> Fill(ibin, m_Ljetacut);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "#Ljeta cut");
    m_CutValsHisto -> Fill(ibin, m_Ljptcut);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "Ljet p_{T} cut");

    m_CutValsHisto -> Fill(ibin, kTau21WCut2);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "W #tau_{21} cut2");
    m_CutValsHisto -> Fill(ibin, kmWcut1);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "m_{W} cut1");
    m_CutValsHisto -> Fill(ibin, kmWcut2);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "m_{W} cut2");

    m_CutValsHisto -> Fill(ibin, kTau32tCut1);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "W #tau_{32} cut2");
    m_CutValsHisto -> Fill(ibin, kTau21tCut1);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "W #tau_{21} cut2");
    m_CutValsHisto -> Fill(ibin, kTau21tCut2);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "W #tau_{21} cut2");
    m_CutValsHisto -> Fill(ibin, kmtcut1);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "m_{T} cut1");
    m_CutValsHisto -> Fill(ibin, kmtcut2);
    m_CutValsHisto -> GetXaxis() -> SetBinLabel(ibin++, "m_{T} cut2");

    cout << "Writing output..." << endl;
    outfile->Write();

    cout << "DONE! Enjoy the output!;-)" << endl;

}

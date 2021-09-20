#include "HistoMaker.h"
// JK 2017, 2018, 2020


// __________________________________________________
double* AllocateArray(int n, double *arr)
{
  //cout << "   initializing by " << endl;;
  double *vals = new double[n+1]; // yes, n+1 !!!
    for (unsigned int i = 0; i <= n; ++i) // !!! note the <= here, as n is already number of bins, not edges!;-)
    {
        //cout << "  " << arr[i] << ", ";
        vals[i] = arr[i];
    }
    //cout << endl;
    return vals;
}

// __________________________________________________
double* MakeUniformBins(int n, double x1, double x2)
{
    double *vals = new double[n+1]; // yes, n+1 !!!
    for (unsigned int i = 0; i <= n; ++i) // !!! note the <= here, as n is already number of bins, not edges!;-)
    {
        vals[i] = x1 + i*(x2 - x1) / n;
    }
    //cout << endl;
    return vals;

}

// __________________________________________________
int ComputeMaxJetsToFill(TString JetType)
{

    int nMaxJetsToFill = kMaxJetsToFill;
    if (JetType == "Top") {
        nMaxJetsToFill = 2;
    } else if (JetType == "FourTop") {
            nMaxJetsToFill = 0;
    } else if (JetType == "DiTop") {
        nMaxJetsToFill = 1;
    } else if (JetType == "W") {
        nMaxJetsToFill = 2;
    }

    return nMaxJetsToFill;

}

// __________________________________________________

void HistoMaker::SetReplicaWeights(double *weights) {
    for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
        m_replicaWeights[irep] = weights[irep];
    }

}


// __________________________________________________

HistoMaker::HistoMaker(TString m_level, int nReplicas) : m_level(m_level), m_nReplicas(nReplicas)
{

    // 2.9.2020

    // see also matching TestGenCorrCode!!! jk 30.11.2020

    m_replicaWeights = new double[m_nReplicas];
    for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
       m_replicaWeights[irep] = 1.;
    }

    TString GeV = " [GeV]";
    int ivar = 0;
    m_vars["DiTopPout"] = varProps("|p_{out}|", ivar++, GeV, "pseudo.Pout[0]");
    m_vars["DiTopDeltaPhi"] = varProps("#Delta#phi^{t#bar{t}}", ivar++, "", "pseudo.DeltaPhi");
    m_vars["DiTopMass"] = varProps("m^{t#bar{t}}", ivar++, GeV, "pseudo.pseudottbar.M()");
    m_vars["DiTopPt"] = varProps("p^{t#bar{t}}_{T}", ivar++, GeV, "pseudo.pseudottbar.Pt()");
    m_vars["TopPt"] = varProps("p_{T}^{t}", ivar++, GeV, "pseudo.pseudotophadron.Pt()");
    m_vars["CosThetaStar"] = varProps("|cos#theta*|", ivar++, "", "pseudo.CosThetaStar[0]");
    m_vars["Delta"] = varProps("#delta^{t#bar{t}}", ivar++, "", "pseudo.deltattbar");
    m_vars["Chittbar"] = varProps("#chi^{t#bar{t}}", ivar++, "", "pseudo.Chittbar");
    m_vars["Yboost"] = varProps("y_{boost}^{t#bar{t}}", ivar++, "", "pseudo.Yboost");
    m_vars["Rttbar"] = varProps("R^{t2,t1}", ivar++, "", "pseudo.Rttbar");

    m_vars["jApla"] = varProps("jets Aplanarity", ivar++, "", "GlobalVars.jApla");
    m_vars["jSpher"] = varProps("jets Sphericity", ivar++, "", "GlobalVars.jSpher");
    //m_vars["JApla"] = varProps("Jets Aplanarity", ivar++, "", "GlobalVars.JApla");
    //m_vars["JSpher"] = varProps("Jets Sphericity", ivar++, "", "GlobalVars.JSpher");
    m_vars["HTj"] = varProps("H_{T}^{j}", ivar++, GeV, "GlobalVars.HTj");
    //m_vars["SumMj"] = varProps("#sum m^{j}", ivar++, GeV, "GlobalVars.SumMj");
    //m_vars["HTJ"] = varProps("H_{T}^{J}", ivar++, GeV, "GlobalVars.HTJ");
    m_vars["HTjPlusMet"] = varProps("H_{T}^{j} + E^{miss}_{T}", ivar++, GeV, "GlobalVars.HTj + GlobalVars.met");
    //m_vars["HTJPlusMet"] = varProps("H_{T}^{J} + E^{miss}_{T}", ivar++, GeV, "GlobalVars.HTJ + GlobalVars.met");
    m_vars["SumMJ"] = varProps("#sum m^{J}", ivar++, GeV, "GlobalVars.SumMJ");

    m_vars["Met"] = varProps("E^{miss}_{T}", ivar++, GeV, "GlobalVars.met");
    //m_vars["MjVis"] = varProps("m_{sum j}^{vis}", ivar++, GeV, "GlobalVars.MjVis");
    m_vars["MJVis"] = varProps("m_{sum J}^{vis}", ivar++, GeV, "GlobalVars.MJVis");

    m_vars["DiTopPtRel"] = varProps("p_{T}^{t#bar{t}} / m^{t#bar{t}}", ivar++, "", "pseudo.DiTopPtRel");
    m_vars["DiTopPtGeo"] = varProps("p_{T}^{t#bar{t}} / #sqrt{p_{T}^{t1} p_{T}^{t2}}", ivar++, "", "pseudo.DiTopPtGeo");
    m_vars["DiTopMassGeo"] = varProps("m^{t#bar{t}} / #sqrt{p_{T}^{t1} p_{T}^{t2}}", ivar++, "", "pseudo.DiTopMassGeo");
    m_vars["DiTopPoutRel"] = varProps("|p_{out}| / m^{t#bar{t}}", ivar++, "", "pseudo.DiTopPoutRel[0]");
    m_vars["DiTopPoutGeo"] = varProps("|p_{out}| / #sqrt{p_{T}^{t1} p_{T}^{t2}}", ivar++, "", "pseudo.DiTopPoutGeo[0]");
    m_vars["TopPtRel"] = varProps("p_{T}^{t} / m^{t#bar{t}}", ivar++, "", "pseudo.TopPtRel[0]");

    cout << "Will study correlations among " << ivar << " variables!" << endl;

    m_nTauBins = 150;
    m_Taumin = 0;
    m_Taumax = 1.1;


    // template:
    //    m_nPhysObjBins[""] = ;
    //    m_PhysObjBins[""] = MakeUniformBins(m_nPhysObjBins[""], );

    m_nPhysObjBins["CosThetaStar"] = 12;
    m_PhysObjBins["CosThetaStar"] = MakeUniformBins(m_nPhysObjBins["CosThetaStar"], 0., 1.);

    m_nPhysObjBins["Delta"] = 10.;
    m_PhysObjBins["Delta"] = MakeUniformBins(m_nPhysObjBins["CosThetaStar"], 0., TMath::Pi());

    m_nPhysObjBins["Chittbar"] = 10;
    m_PhysObjBins["Chittbar"] = MakeUniformBins(m_nPhysObjBins["Chittbar"], 1., 20.);

    m_nPhysObjBins["Yboost"] = 10;
    m_PhysObjBins["Yboost"] = MakeUniformBins(m_nPhysObjBins["Yboost"], 0., 2.);
    // leading over subleading
    // double rttbins[] = {1., 1.1, 1.2, 1.4, 1.6, 1.9, 2.5, 3.5, 5.};
    // subleading oved leading
    double rttbins[] = {0., 0.2, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.90, 0.95, 1.001};
    m_nPhysObjBins["Rttbar"] = sizeof(rttbins)/ sizeof(double) - 1;
    m_PhysObjBins["Rttbar"] = AllocateArray(m_nPhysObjBins["Rttbar"], rttbins); //MakeUniformBins(m_nPhysObjBins["Rttbar"], 0.5, 3.);


    // non-uniform bins:
    double dphibins[] = {0., 1., 2., 2.34, 2.54, 2.69, 2.84, 2.94, 3.04, TMath::Pi()};
    m_nPhysObjBins["DiTopDeltaPhi"] = sizeof(dphibins) / sizeof(double) - 1;
    m_PhysObjBins["DiTopDeltaPhi"] = AllocateArray(m_nPhysObjBins["DiTopDeltaPhi"], dphibins);

    double poutbins[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600,  800}; // 550; 700
    m_nPhysObjBins["DiTopPout"] = sizeof(poutbins) / sizeof(double) - 1;
    m_PhysObjBins["DiTopPout"] = AllocateArray(m_nPhysObjBins["DiTopPout"], poutbins);

    double rapiditybins[] = {-2.5, -1.8, -1.6, -1.4, -1.2, -1., -0.8, -0.6, -0.4, -0.2,
                             0.,
                             0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.5};
    m_nPhysObjBins["Rapidity"] = sizeof(rapiditybins) / sizeof(double) - 1;
    m_PhysObjBins["Rapidity"] = AllocateArray(m_nPhysObjBins["Rapidity"], rapiditybins);


    // large-R jet pT and mass bins:
    double LJetPtbins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                           120, 140, 160, 180, 200,
                           233, 266, 300, 333, 366,
                           400, 433, 466, 500, 533, 566, 600, 650, 700, 800, 900, 1000};
    m_nPhysObjBins["LJetPt"] = sizeof(LJetPtbins) / sizeof(double) - 1;
    m_PhysObjBins["LJetPt"] = AllocateArray(m_nPhysObjBins["LJetPt"], LJetPtbins);

    double LJetMbins[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115,
                          120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230,
                          240, 250, 260, 270, 280, 290, 300,
                          350, 400, 500};
    m_nPhysObjBins["LJetMass"] = sizeof(LJetMbins) / sizeof(double) - 1;
    m_PhysObjBins["LJetMass"] = AllocateArray(m_nPhysObjBins["LJetMass"], LJetMbins);

    double ditopmassbins_0B2S[] = {340.0, 500.0, // 400; 450; 500
                                   //550.0, 600.0, 650.0, 700.0,
                                   //750.0, 800.0, 850.0, 900.0,
                                   //950.0, 1000.0, 1050.0, 1100.0, 1150.0,
                                   600., 700.,
                                   800.0,  900.0,
                                   1000.0,  1100.0,
                                   1200.0,// 1300.0,
                                   1400.,
                                   1600.,
                                   // 1500.0, 1750.,
                                   2000.};
    double ditopmassbins_1B1S[] = {500,
                                   //600.,
                                   700.,
                                   800.0,  900.0,
                                   1000.0,  1100.0,
                                   1200.0, // 1300.0,
                                   1400.,
                                   1600.,
                                   //1500.0, 1750.,
                                   2000.};
    double ditopmassbins_2B0S[] = {700.0,
                                   // 750.0, 800.0, 850.0, 900.0,
                                   // 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
                                   800.0,  900.0,
                                   1000.0,  1100.0,
                                   1200.0, //1300.0,
                                   1400.,
                                   1600.,
                                   //1500.0, 1750.,
                                   2000.};
    double ditopmassbins_AnySel[] = {340., 500, 600, 700, 
                                   800.0,  900.0,
                                   1000.0,  1100.0,
                                   1200.0, //1300.0,
                                   1400.,
                                   1600.,
                                   //1500.0, 1750.,
                                   2000.};

    double mvisbins[] = {340., 500, 600, 700, 
                                   800.0,  900.0,
                                   1000.0,  1100.0,
                                   1200.0, //1300.0,
                                   1400.,
                                   1600.,
                                   //1500.0, 1750.,
                                   2000.,3000., 4000., 5000., 8000.};

    m_nPhysObjBins["DiTopMass_0B2S"] = sizeof(ditopmassbins_0B2S) / sizeof(double) - 1;
    m_PhysObjBins["DiTopMass_0B2S"] = AllocateArray(m_nPhysObjBins["DiTopMass_0B2S"], ditopmassbins_0B2S);
    m_nPhysObjBins["DiTopMass_1B1S"] = sizeof(ditopmassbins_1B1S) / sizeof(double) - 1;
    m_PhysObjBins["DiTopMass_1B1S"] = AllocateArray(m_nPhysObjBins["DiTopMass_1B1S"], ditopmassbins_1B1S);
    m_nPhysObjBins["DiTopMass_2B0S"] = sizeof(ditopmassbins_2B0S) / sizeof(double) - 1;
    m_PhysObjBins["DiTopMass_2B0S"] = AllocateArray(m_nPhysObjBins["DiTopMass_2B0S"], ditopmassbins_2B0S);
    m_nPhysObjBins["DiTopMass_AnySel"] = sizeof(ditopmassbins_AnySel) / sizeof(double) - 1;
    m_PhysObjBins["DiTopMass_AnySel"] = AllocateArray(m_nPhysObjBins["DiTopMass_AnySel"], ditopmassbins_AnySel);
    // 1.12.2020
    m_nPhysObjBins["MVis"] = sizeof(mvisbins) / sizeof(double) - 1;
    m_PhysObjBins["MVis"] = AllocateArray(m_nPhysObjBins["MVis"], mvisbins);


    double topptbins_0B2S[] = {0., 100.,
                               150.0, 200,
                               250.0,  300.0,
                               350.0, 400.0,
                               450.0, 500.0,
                               550.0, 600.0 //, 700.0, 1000.0
                              };
    double topptbins_1B1S[] = {100., 210.,
                               250.0, 300.0,
                               350.0, 400.0,
                               450.0, 500,
                               550.0, 600.0,
                               700.0, 1000.0};
    double topptbins_2B0S[] = {   240.0,
                                  350.0,
                                  400, 450.0, //480.0, 510.0, 540.0, 570.0, 600.0,
                                  500, 600,
                                  700.0, 800.0, 1000.0};
    double topptbins_AnySel[] = {0., 100.,
                               150.0, 200,
                               250.0,  300.0,
                               350.0, 400.0,
                               450.0, 500.0,
                               550.0, 600.0, 700.0, 800., 1000.0
                              };

    m_nPhysObjBins["TopPt_0B2S"] = sizeof(topptbins_0B2S) / sizeof(double) - 1;
    m_PhysObjBins["TopPt_0B2S"] = AllocateArray(m_nPhysObjBins["TopPt_0B2S"], topptbins_0B2S);
    m_nPhysObjBins["TopPt_2B0S"] = sizeof(topptbins_2B0S) / sizeof(double) - 1;
    m_PhysObjBins["TopPt_2B0S"] = AllocateArray(m_nPhysObjBins["TopPt_2B0S"], topptbins_2B0S);
    m_nPhysObjBins["TopPt_1B1S"] = sizeof(topptbins_1B1S) / sizeof(double) - 1;
    m_PhysObjBins["TopPt_1B1S"] = AllocateArray(m_nPhysObjBins["TopPt_1B1S"], topptbins_1B1S);
    m_nPhysObjBins["TopPt_AnySel"] = sizeof(topptbins_AnySel) / sizeof(double) - 1;
    m_PhysObjBins["TopPt_AnySel"] = AllocateArray(m_nPhysObjBins["TopPt_AnySel"], topptbins_AnySel);

    //                    so far general for all topologies:

    double ditopptbins[] = {0., 30., 60.0, 90.,
                            120.0, 160.0, 200.0,
                            250.0, 300.0,
                            350.0, 450.0,  600.0}; // 400.0,  ; 500.0,
    m_nPhysObjBins["DiTopPt"] = sizeof(ditopptbins) / sizeof(double) - 1;
    m_PhysObjBins["DiTopPt"] = AllocateArray(m_nPhysObjBins["DiTopPt"], ditopptbins);

    double ditoppoutrel[] = {0., 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 1.};
    m_nPhysObjBins["DiTopPoutRel"] = sizeof(ditoppoutrel) / sizeof(double) - 1;
    m_PhysObjBins["DiTopPoutRel"] = AllocateArray(m_nPhysObjBins["DiTopPoutRel"], ditoppoutrel);

    double ditoppoutgeo [] = {0., 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.5};
    m_nPhysObjBins["DiTopPoutGeo"] = sizeof(ditoppoutgeo) / sizeof(double) - 1;
    m_PhysObjBins["DiTopPoutGeo"] = AllocateArray(m_nPhysObjBins["DiTopPoutGeo"], ditoppoutgeo);

    double ditopptrel[] = {0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.5};
    m_nPhysObjBins["DiTopPtRel"] = sizeof(ditopptrel) / sizeof(double) - 1;
    m_PhysObjBins["DiTopPtRel"] = AllocateArray(m_nPhysObjBins["DiTopPtRel"], ditopptrel);

    double ditopptgeo[] = {0.,  0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.};
    m_nPhysObjBins["DiTopPtGeo"] = sizeof(ditopptgeo) / sizeof(double) - 1;
    m_PhysObjBins["DiTopPtGeo"] = AllocateArray(m_nPhysObjBins["DiTopPtGeo"], ditopptgeo);

    double ditopmassgeo[] = {1.5, 2.5, 3., 3.5, 4., 4.5, 5., 7., 10., 20.};
    m_nPhysObjBins["DiTopMassGeo"] = sizeof(ditopmassgeo) / sizeof(double) - 1;
    m_PhysObjBins["DiTopMassGeo"] = AllocateArray(m_nPhysObjBins["DiTopMassGeo"], ditopmassgeo);

    double topptrel[] = {0.,  0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.75, 1.};
    m_nPhysObjBins["TopPtRel"] = sizeof(topptrel) / sizeof(double) - 1;
    m_PhysObjBins["TopPtRel"] = AllocateArray(m_nPhysObjBins["TopPtRel"], topptrel);


    double spher[] = {0.,  0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.75, 1.};
    m_nPhysObjBins["Spher"] = sizeof(spher) / sizeof(double) - 1;
    m_PhysObjBins["Spher"] = AllocateArray(m_nPhysObjBins["Spher"], spher);

    double apla[] = {0.,  0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.2, 0.3, 0.5};
    m_nPhysObjBins["Apla"] = sizeof(apla) / sizeof(double) - 1;
    m_PhysObjBins["Apla"] = AllocateArray(m_nPhysObjBins["Apla"], apla);

    double metbins[] = {0., 10, 20.0, 30, 40.0, 50, 60.0, 80.0, 100.0,
                        120.0, 140.0, 160.0, 180.0, 200.0, 500.};
    m_nPhysObjBins["Met"] = sizeof(metbins) / sizeof(double) - 1;
    m_PhysObjBins["Met"] = AllocateArray(m_nPhysObjBins["Met"], metbins);


    double htbins[] = {0.,  50.0,  100.0,
                       150.0,  200.0,
                       250.0, 300.0,
                        350.0,  400.0,
                        450.0,  500.0,
                        550.0, 600.0,
                        650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 1000.0, 1100, 1200, 1300, 1500, 2000, 3000};
    m_nPhysObjBins["HT"] = sizeof(htbins) / sizeof(double) - 1;
    m_PhysObjBins["HT"] = AllocateArray(m_nPhysObjBins["HT"], htbins);

    double SumMjbins[] = {0.,  50.0,  100.0,
                          150.0,  200.0,
                          250.0, 300.0,
                           350.0,  400.0,
                           450.0,  500.0,
                           550.0, 600.0,
                           650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 1000.0};
    m_nPhysObjBins["SumM"] = sizeof(SumMjbins) / sizeof(double) - 1;
    m_PhysObjBins["SumM"] = AllocateArray(m_nPhysObjBins["SumM"], htbins);

    TString sels[] = {"0B2S", "1B1S", "2B0S", "AnySel"};

    // and now some postprocessing, and if missing the topology-specific versions of binning, make them!
    // 6.8.2020
    for (auto& [var, varProp]: m_vars) {
        for (auto sel : sels) {
            TString bvar = var; // remove j a J aka small and LargeJets specific tags in var names
            bvar.ReplaceAll("J","");
            bvar.ReplaceAll("j","");
            if (bvar.Contains("PlusMet"))
                continue;
            if (m_PhysObjBins.find(bvar + "_" + sel) == m_PhysObjBins.end()) {
                //cout << "Creating bins for " << (bvar + "_" + sel).Data() << endl;
                m_PhysObjBins[bvar + "_" + sel] = m_PhysObjBins[bvar];
                m_nPhysObjBins[bvar + "_" + sel] = m_nPhysObjBins[bvar];
            }
        }
    }

    /*
    for (auto sel : sels) {
        TString vars[] = {"DiTopMass", "TopPt"};
        for (auto var : vars) {
            TString key = var + "_" + sel;
            cout << "Using " << m_nPhysObjBins[key] << " " << key.Data() << " bins:" << endl;
            for (int i = 0; i <= m_nPhysObjBins[key]; ++i)
                cout << " " << m_PhysObjBins[key][i] << ", ";
            cout << endl;
        }
    }
    cout << "Using " << m_nPhysObjBins["DiTopPt"] << " DiTopPt bins:" << endl;
    for (int i = 0; i <= m_nPhysObjBins["DiTopPt"]; ++i)
        cout << " " << m_PhysObjBins["DiTopPt"][i] << ", ";
    cout << endl;
    */
}

// __________________________________________________

HistoMaker::~HistoMaker()
{
    if (m_nReplicas > 0) {
        delete [] m_replicaWeights;
    }
}

// __________________________________________________

TString HistoMaker::GetLevel() {
    return m_level;
}


// __________________________________________________
	
void HistoMaker::AddCutLevel(TString cutname, bool jetsOnly)
{
  if ( m_holders.find(cutname) == m_holders.end() ) {
    m_holders[cutname] = new HistoHolder(cutname, jetsOnly);
  }
  if (!gDirectory->GetDirectory(cutname))
      gDirectory->mkdir(cutname);
}
// __________________________________________________

bool HistoMaker::PassedJetType(kJetTypes jettype, TMyLorentzVector *jet)
{
    if (jettype == kbJets && ! jet->Btag() )
        return false;
    if (jettype == knonbJets && jet->Btag() )
        return false;
    if (jettype == kTopTagJets && ! jet->Toptag() )
        return false;
    if (jettype == kWTagJets && ! jet->Wtag() )
        return false;
    return true;

}

// __________________________________________________

void HistoMaker::MakeSpecialHistos(HistoHolder *hold, TString jname, TString jtag, bool makeReplicas)
{

    TString GeV = " [GeV]";

    //    cout << "In MakeSpecialHistos, jname=" << jname.Data() << " jtag=" << jtag.Data() << endl;

    //   double pOmax = 1000;
    // hold -> AddTH1D(jname + "Pout", jname + "Pout;" + jtag + " |p_{out}|" + GeV, nPtBins, 0, pOmax);
    hold -> AddTH1D(jname + "Pout", jname + "Pout;" + jtag + " |p_{out}|" + GeV, m_nPhysObjBins["DiTopPout"], m_PhysObjBins["DiTopPout"]);
    hold -> AddTH1D(jname + "DeltaPhi", jname + "DeltaPhi;" + jtag + " #Delta#phi^{t#bar{t}}", m_nPhysObjBins["DiTopDeltaPhi"], m_PhysObjBins["DiTopDeltaPhi"]);
    hold -> AddTH1D(jname + "CosThetaStar", jname + "CosThetaStar;" + jtag + " |cos#theta*|", m_nPhysObjBins["CosThetaStar"], m_PhysObjBins["CosThetaStar"]); // abs?!
    hold -> AddTH1D(jname + "Delta", jname + "Delta;" + jtag + " #delta^{t#bar{t}}", m_nPhysObjBins["Delta"], m_PhysObjBins["Delta"]);
    hold -> AddTH1D(jname + "Chittbar", jname + "Chittbar;" + jtag + " #chi^{t#bar{t}}", m_nPhysObjBins["Chittbar"], m_PhysObjBins["Chittbar"]);
    hold -> AddTH1D(jname + "Yboost", jname + "Yboost;" + jtag + " y_{boost}^{t#bar{t}}", m_nPhysObjBins["Yboost"], m_PhysObjBins["Yboost"]);
    hold -> AddTH1D(jname + "Rttbar", jname + "Rttbar;" + jtag + " R^{t2,t1}", m_nPhysObjBins["Rttbar"], m_PhysObjBins["Rttbar"]);


    if (makeReplicas) {
        TString repdir = "replicas";
        if (!gDirectory->GetDirectory(repdir)) {
            gDirectory->mkdir(repdir);
        }
        gDirectory->cd(repdir);
        for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
            TString reptag = Form("_rep%i",irep);
            // hold -> AddTH1D(jname + "Pout", jname + "Pout;" + jtag + " |p_{out}|" + GeV, nPtBins, 0, pOmax);
            hold -> AddTH1D(jname + "Pout" + reptag, jname + "Pout;" + jtag + " |p_{out}|" + GeV, m_nPhysObjBins["DiTopPout"], m_PhysObjBins["DiTopPout"]);
            hold -> AddTH1D(jname + "DeltaPhi" + reptag, jname + "DeltaPhi;" + jtag + " #Delta#phi^{t#bar{t}}", m_nPhysObjBins["DiTopDeltaPhi"], m_PhysObjBins["DiTopDeltaPhi"]);
            hold -> AddTH1D(jname + "CosThetaStar" + reptag, jname + "CosThetaStar;" + jtag + " |cos#theta*|", m_nPhysObjBins["CosThetaStar"], m_PhysObjBins["CosThetaStar"]); // abs?!
            hold -> AddTH1D(jname + "Delta" + reptag, jname + "Delta;" + jtag + " #delta^{t#bar{t}}", m_nPhysObjBins["Delta"], m_PhysObjBins["Delta"]);
            hold -> AddTH1D(jname + "Chittbar" + reptag, jname + "Chittbar;" + jtag + " #chi^{t#bar{t}}", m_nPhysObjBins["Chittbar"], m_PhysObjBins["Chittbar"]);
            hold -> AddTH1D(jname + "Yboost" + reptag, jname + "Yboost;" + jtag + " y_{boost}^{t#bar{t}}", m_nPhysObjBins["Yboost"], m_PhysObjBins["Yboost"]);
            hold -> AddTH1D(jname + "Rttbar" + reptag, jname + "Rttbar;" + jtag + " R^{t1,t2}", m_nPhysObjBins["Rttbar"], m_PhysObjBins["Rttbar"]);

        }
        gDirectory->cd("../");
    }

    // additional variables:

    // new ratio variables:
    hold -> AddTH1D(jname + "PoutRel", jname + "PoutRel;" + jtag + " |p_{out}| / m^{t#bar{t}}", m_nPhysObjBins["DiTopPoutRel"], m_PhysObjBins["DiTopPoutRel"]);
    hold -> AddTH1D(jname + "PoutGeo", jname + "PoutGeo;" + jtag + " |p_{out}| / #sqrt{p_{T}^{t1} p_{T}^{t2}}", m_nPhysObjBins["DiTopPoutGeo"], m_PhysObjBins["DiTopPoutGeo"]);
    hold -> AddTH1D(jname + "PtRel", jname + "PtRel;" + jtag + " p_{T}^{t#bar{t}} / m^{t#bar{t}}", m_nPhysObjBins["DiTopPtRel"], m_PhysObjBins["DiTopPtRel"]);
    hold -> AddTH1D(jname + "PtGeo", jname + "PtGeo;" + jtag + " p_{T}^{t#bar{t}} / #sqrt{p_{T}^{t1} p_{T}^{t2}}", m_nPhysObjBins["DiTopPtGeo"], m_PhysObjBins["DiTopPtGeo"]);
    hold -> AddTH1D(jname + "MassGeo", jname + "MassGeo;" + jtag + " m^{t#bar{t}} / #sqrt{p_{T}^{t1} p_{T}^{t2}}", m_nPhysObjBins["DiTopMassGeo"], m_PhysObjBins["DiTopMassGeo"]);

    // create top pT geo histos in jet histos section?
    // also other jet pTs? and pT themselves in TeV??

    // now Rel and Geo for Top1 and Top2 pT, by replacing DiTop:
    jname.ReplaceAll("DiTop", "Top1");
    jtag.ReplaceAll("DiTop", "Top1");
    hold -> AddTH1D(jname + "PtRel", jname + "PtRel;" + jtag + " p_{T}^{t1} / m^{t#bar{t}}", m_nPhysObjBins["TopPtRel"], m_PhysObjBins["TopPtRel"]);
    jname.ReplaceAll("Top1", "Top2");
    jtag.ReplaceAll("Top1", "Top2");
    hold -> AddTH1D(jname + "PtRel", jname + "PtRel;" + jtag + " p_{T}^{t2} / m^{t#bar{t}}", m_nPhysObjBins["TopPtRel"], m_PhysObjBins["TopPtRel"]);
    jname.ReplaceAll("Top2", "Top");
    jtag.ReplaceAll("Top2", "Top");
    hold -> AddTH1D(jname + "PtRel", jname + "PtRel;" + jtag + " p_{T}^{t} / m^{t#bar{t}}", m_nPhysObjBins["TopPtRel"], m_PhysObjBins["TopPtRel"]);


}

// __________________________________________________

void HistoMaker::MakeGlobalHistos(HistoHolder *hold)
{

    if (!gDirectory->GetDirectory(hold->GetName())) {
        gDirectory->mkdir(hold->GetName());
    }
    gDirectory->cd(hold->GetName());

    hold -> AddTH1D(m_level + "jAplanarity", m_level + "jAplanarity;" + m_level + " jets Aplanarity",  m_nPhysObjBins["Apla"], m_PhysObjBins["Apla"]);
    hold -> AddTH1D(m_level + "jSphericity", m_level + "jSphericity;" + m_level + " jets Sphericity",  m_nPhysObjBins["Spher"], m_PhysObjBins["Spher"]);
    hold -> AddTH1D(m_level + "JAplanarity", m_level + "JAplanarity;" + m_level + " LJets Aplanarity",  m_nPhysObjBins["Apla"], m_PhysObjBins["Apla"]);
    hold -> AddTH1D(m_level + "JSphericity", m_level + "JSphericity;" + m_level + " LJets Sphericity",  m_nPhysObjBins["Spher"], m_PhysObjBins["Spher"]);

    hold -> AddTH1D(m_level + "Met", m_level + "Met;" + m_level + " E^{miss}_{T}", m_nPhysObjBins["Met"], m_PhysObjBins["Met"]);
    hold -> AddTH1D(m_level + "MetPhi", m_level + "MetPhi;" + m_level + " #phi_{E^{miss}_{T}}", 100, -TMath::Pi(), TMath::Pi());
    hold -> AddTH1D(m_level + "HTj", m_level + "HTj;" + m_level + " H_{T}^{j}", m_nPhysObjBins["HT"], m_PhysObjBins["HT"]);
    hold -> AddTH1D(m_level + "HTJ", m_level + "HTJ;" + m_level + " H_{T}^{J}", m_nPhysObjBins["HT"], m_PhysObjBins["HT"]);
    hold -> AddTH1D(m_level + "HTjPlusMet", m_level + "HTjPlusMet;" + m_level + " H_{T}^{j} + MET", m_nPhysObjBins["HT"], m_PhysObjBins["HT"]);
    hold -> AddTH1D(m_level + "HTJPlusMet", m_level + "HTJPlusMet;" + m_level + " H_{T}^{J} + MET", m_nPhysObjBins["HT"], m_PhysObjBins["HT"]);
    hold -> AddTH1D(m_level + "SumMj", m_level + "SumMj;" + m_level + " #sum m^{j}", m_nPhysObjBins["SumM"], m_PhysObjBins["SumM"]);
    hold -> AddTH1D(m_level + "SumMJ", m_level + "SumMJ;" + m_level + " #sum m^{J}", m_nPhysObjBins["SumM"], m_PhysObjBins["SumM"]);

    hold -> AddTH1D(m_level + "MJVis", m_level + "MJVis;" + m_level + " Visible mass J", m_nPhysObjBins["MVis"], m_PhysObjBins["MVis"]);
    hold -> AddTH1D(m_level + "MjVis", m_level + "MjVis;" + m_level + " Visible mass j", m_nPhysObjBins["MVis"], m_PhysObjBins["MVis"]);

    gDirectory->cd("../");
}

// __________________________________________________

void HistoMaker::MakeSingleObjHistos(TString ObjName, HistoHolder *hold)
{
    if (!gDirectory->GetDirectory(hold->GetName())) {
        gDirectory->mkdir(hold->GetName());
    }
    gDirectory->cd(hold->GetName());


	cout << "Making kinematic histos for " << ObjName.Data() << endl;

    TString jname = m_level + "_" + ObjName;
    TString jtag = m_level + " " + ObjName;
	TString GeV = " [GeV]";

	// HERE	
	int nPtBins = 100;
	double pTmin = 0;
	double pTmax = 500;


	if (ObjName.Contains("pseudoW")) {
		pTmax = 500;
	}
    if (ObjName.Contains("ttbar") || ObjName.Contains("DiTop") ) {
       this -> MakeSpecialHistos(hold, jname, jtag, false); // do NOT create replicas
    }
	hold -> AddTH1D(jname + "Pt", jname + "Pt;" + jtag + " p_{T}" + GeV, nPtBins, pTmin, pTmax);
	
	int nRapidityBins = 40;
	double rapmin = -2.5;
	double rapmax = 2.5;

	
	int nMassBins = 100;
	double Massmin = 0;
	double Massmax = 1400;
    if (ObjName.Contains("ttbar")) {
	  nRapidityBins = 40;
	  Massmin = 300;
	  Massmax = 1300;
	  nMassBins = 250;
	}

	if (ObjName.Contains("pseudotop")) {
	  nRapidityBins = 40;
	  Massmin = 100;
	  Massmax = 300;
	  nMassBins = 100;

	}
	if (ObjName.Contains("pseudoW")) {
	  nRapidityBins = 40;
	  Massmax = 160;
	  nMassBins = 100;
	}
	hold -> AddTH1D(jname + "Rapidity", jname + "Rapidity;" + jtag + " y", nRapidityBins, rapmin, rapmax);
	hold -> AddTH1D(jname + "Mass", jname + "Mass;" + jtag + " m" + GeV, nMassBins, Massmin, Massmax);
  
	gDirectory->cd("../");
}



// __________________________________________________

void HistoMaker::MakeLjetHistos(TString cutname)
{
  HistoHolder* hold = 0;
  try {
    hold = m_holders[cutname];
  }  catch (exception& e) {
    cout << "ERROR in HistoMaker::MakeLjetHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
   }

  this -> MakeSingleObjHistos("pseudolepton", hold);
  this -> MakeSingleObjHistos("pseudoneutrino", hold);
  this -> MakeSingleObjHistos("pseudoWlepton", hold);
  this -> MakeSingleObjHistos("pseudoWhadron", hold);
  this -> MakeSingleObjHistos("pseudotoplepton", hold);
  this -> MakeSingleObjHistos("pseudotophadron", hold);
  this -> MakeSingleObjHistos("pseudottbar", hold);    
}

// __________________________________________________

void HistoMaker::MakeLjetMigrations(TString cutname, TString m_level1, TString m_level2)
{

    auto iter = m_holders.find(cutname);
    if (iter != m_holders.end()) {
        auto hold = iter->second;
        if (!gDirectory->GetDirectory(hold->GetName()))
            gDirectory->mkdir(hold->GetName());
        gDirectory->cd(hold->GetName());
        TString migradir = "migrations";
        if (!gDirectory->GetDirectory(migradir)) {
            gDirectory->mkdir(migradir);
        }
        gDirectory->cd(migradir);

        this -> MakeJetMigrations(cutname, "pseudolepton", m_level1, m_level2, false, 0);
        this -> MakeJetMigrations(cutname, "pseudoneutrino", m_level1, m_level2, false, 0);
        this -> MakeJetMigrations(cutname, "pseudoWlepton", m_level1, m_level2, false, 0);
        this -> MakeJetMigrations(cutname, "pseudoWhadron", m_level1, m_level2, false, 0);
        this -> MakeJetMigrations(cutname, "pseudotoplepton", m_level1, m_level2, false, 0);
        this -> MakeJetMigrations(cutname, "pseudotophadron", m_level1, m_level2, false, 0);
        this -> MakeJetMigrations(cutname, "pseudottbar", m_level1, m_level2, false, 0);

        gDirectory->cd("../");
        gDirectory->cd("../");

    } else
    cerr << "ERROR HistoMaker::MakeLjetMigrations: making dijet histos!" << endl;
	
  return;
  
}


// __________________________________________________

void HistoMaker::FillTTbarSpecialHistos(TString cutname, TString DiTopName, kPseudotop pseudo, kGlobalVars GlobalVars,
                                        double weight, bool fillReplicas)
{

    HistoHolder* hold = 0;
    try {
      hold = m_holders[cutname];
    }  catch (exception& e) {
      cout << "ERROR in HistoMaker::MakeLjetHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
      return;
     }
    TString addtag = "";
    if (DiTopName.Contains("pseudo"))
        addtag = "_";

    //cout << "Filling level " << m_level.Data() << endl;
    //cout << "Filling with weights of 0.5 pouts of " << pseudo.Pout[0] << " and " << pseudo.Pout[1] << endl;
    hold -> FillTH1D(m_level + addtag + DiTopName + "Pout", fabs(pseudo.Pout[0]), weight*0.5);
    hold -> FillTH1D(m_level + addtag + DiTopName + "Pout", fabs(pseudo.Pout[1]), weight*0.5);
    hold -> FillTH1D(m_level + addtag + DiTopName + "DeltaPhi", pseudo.DeltaPhi, weight);
    hold -> FillTH1D(m_level + addtag + DiTopName + "Yboost", pseudo.Yboost, weight);
    hold -> FillTH1D(m_level + addtag + DiTopName + "Chittbar", pseudo.Chittbar, weight);
    if (pseudo.Rttbar > 0.) {
        hold -> FillTH1D(m_level + addtag + DiTopName + "Rttbar", pseudo.Rttbar, weight);
    }
    hold -> FillTH1D(m_level + addtag + DiTopName + "Delta", pseudo.deltattbar, weight);
    //cout << "Filling with weights of 0.5 cos theta*'s of " << pseudo.CosThetaStar[0] << " and " << pseudo.CosThetaStar[1] << endl;
    hold -> FillTH1D(m_level + addtag + DiTopName + "CosThetaStar", fabs(pseudo.CosThetaStar[0]), weight*0.5);
    hold -> FillTH1D(m_level + addtag + DiTopName + "CosThetaStar", fabs(pseudo.CosThetaStar[1]), weight*0.5);

    // fill replicas here!
    if (m_nReplicas) {
        for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
            TString reptag = Form("_rep%i", irep);

            hold -> FillTH1D(m_level + addtag + DiTopName + "Pout" + reptag, fabs(pseudo.Pout[0]), weight*0.5*m_replicaWeights[irep]);
            hold -> FillTH1D(m_level + addtag + DiTopName + "Pout" + reptag, fabs(pseudo.Pout[1]), weight*0.5*m_replicaWeights[irep]);
            hold -> FillTH1D(m_level + addtag + DiTopName + "DeltaPhi" + reptag, pseudo.DeltaPhi, weight*m_replicaWeights[irep]);
            hold -> FillTH1D(m_level + addtag + DiTopName + "Yboost" + reptag, pseudo.Yboost, weight*m_replicaWeights[irep]);
            hold -> FillTH1D(m_level + addtag + DiTopName + "Chittbar" + reptag, pseudo.Chittbar, weight*m_replicaWeights[irep]);
            if (pseudo.Rttbar > 0) {
                hold -> FillTH1D(m_level + addtag + DiTopName + "Rttbar" + reptag, pseudo.Rttbar, weight*m_replicaWeights[irep]);
            }
            hold -> FillTH1D(m_level + addtag + DiTopName + "Delta" + reptag, pseudo.deltattbar, weight*m_replicaWeights[irep]);
            hold -> FillTH1D(m_level + addtag + DiTopName + "CosThetaStar" + reptag, fabs(pseudo.CosThetaStar[0]), weight*0.5*m_replicaWeights[irep]);
            hold -> FillTH1D(m_level + addtag + DiTopName + "CosThetaStar" + reptag, fabs(pseudo.CosThetaStar[1]), weight*0.5*m_replicaWeights[irep]);

        }
    }


    // special relative 1D histos
    hold -> FillTH1D(m_level + addtag + DiTopName + "PtRel", pseudo.DiTopPtRel, weight);
    hold -> FillTH1D(m_level + addtag + DiTopName + "PtGeo", pseudo.DiTopPtGeo, weight);
    hold -> FillTH1D(m_level + addtag + DiTopName + "MassGeo", pseudo.DiTopMassGeo, weight);
    hold -> FillTH1D(m_level + addtag + DiTopName + "PoutRel", pseudo.DiTopPoutRel[0], weight*0.5);
    hold -> FillTH1D(m_level + addtag + DiTopName + "PoutGeo", pseudo.DiTopPoutGeo[0], weight*0.5);
    hold -> FillTH1D(m_level + addtag + DiTopName + "PoutRel", pseudo.DiTopPoutRel[1], weight*0.5);
    hold -> FillTH1D(m_level + addtag + DiTopName + "PoutGeo", pseudo.DiTopPoutGeo[1], weight*0.5);

    // special correlation histos
    // code generated by ./bin/TestGenCorrCode_x
    // with by-hand split to (not)filling Rtt


    // for top related rel and geo histos
    TString tname = m_level + "Top";

    // make and fill this only for Detector level!!!
    if (m_level.Contains("Detector")) {

        // fill replicas or just the standard case
        // fixed filling the 0th replica, 26.8.2021
        for (int irep = -1; irep < max (0, m_nReplicas); ++irep) {
            TString reptag = Form("_rep%i",irep);
            double rweight = weight;
            if (irep >= 0) {
                rweight *= m_replicaWeights[irep];
            } else {
                reptag = "";
            }

            // non-Rtt filling code:
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.Yboost, rweight);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , GlobalVars.jApla, rweight);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , GlobalVars.jSpher, rweight);
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.Chittbar, rweight*0.5);
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.Chittbar, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.deltattbar, rweight*0.5);
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.deltattbar, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.DiTopMassGeo, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.DiTopMassGeo, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.DiTopPoutGeo[0], rweight*0.25);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.DiTopPoutGeo[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.DiTopPoutGeo[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.CosThetaStar[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.DiTopPoutRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.DiTopPoutRel[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.DiTopPoutRel[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.CosThetaStar[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.DiTopPtGeo, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.DiTopPtGeo, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.DiTopPtRel, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.DiTopPtRel, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , GlobalVars.HTj, rweight*0.5);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , GlobalVars.HTj, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , GlobalVars.HTj + GlobalVars.met, rweight*0.5);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , GlobalVars.HTj + GlobalVars.met, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , GlobalVars.MJVis, rweight*0.5);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , GlobalVars.MJVis, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , GlobalVars.met, rweight*0.5);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , GlobalVars.met, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , GlobalVars.SumMJ, rweight*0.5);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , GlobalVars.SumMJ, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.TopPtRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.TopPtRel[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.TopPtRel[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.CosThetaStar[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.Yboost, rweight*0.5);
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.Yboost, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , GlobalVars.jApla, rweight*0.5);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , GlobalVars.jApla, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , GlobalVars.jSpher, rweight*0.5);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , GlobalVars.jSpher, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.Chittbar, rweight);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.Yboost, rweight);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , GlobalVars.jApla, rweight);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , GlobalVars.jSpher, rweight);
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.Chittbar, rweight);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.CosThetaStar[0], rweight*0.5);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.CosThetaStar[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.deltattbar, rweight);
            hold -> FillTH2D(m_level + "DiTopMass" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.pseudottbar.M(), rweight);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPt" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.pseudottbar.Pt(), rweight);
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.pseudotophadron.Pt(), rweight*0.5);
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.pseudotoplepton.Pt(), rweight*0.5); // add1
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.Yboost, rweight);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , GlobalVars.jApla, rweight);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , GlobalVars.jSpher, rweight);
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.Chittbar, rweight);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.CosThetaStar[0], rweight*0.5);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.CosThetaStar[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.deltattbar, rweight);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPt" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.pseudottbar.Pt(), rweight);
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.pseudotophadron.Pt(), rweight*0.5);
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.pseudotoplepton.Pt(), rweight*0.5); // add1
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.Yboost, rweight);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , GlobalVars.jApla, rweight);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , GlobalVars.jSpher, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopMassGeo" + reptag, pseudo.DiTopMassGeo , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopMassGeo" + reptag, pseudo.DiTopMassGeo , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopMassGeo" + reptag, pseudo.DiTopMassGeo , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopMassGeo" + reptag, pseudo.DiTopMassGeo , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopMassGeo" + reptag, pseudo.DiTopMassGeo , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopMassGeo" + reptag, pseudo.DiTopMassGeo , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.Chittbar, rweight*0.5);
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.Chittbar, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.CosThetaStar[0], rweight*0.25);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.CosThetaStar[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.CosThetaStar[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.Pout[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.deltattbar, rweight*0.5);
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.deltattbar, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopDeltaPhi" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DeltaPhi, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopDeltaPhi" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.DeltaPhi, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopMass" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.pseudottbar.M(), rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopMass" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.pseudottbar.M(), rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DiTopMassGeo, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.DiTopMassGeo, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DiTopPoutGeo[0], rweight*0.25);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.DiTopPoutGeo[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DiTopPoutGeo[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.Pout[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DiTopPoutRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.DiTopPoutRel[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DiTopPoutRel[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.Pout[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "DiTopPt" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.pseudottbar.Pt(), rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPt" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.pseudottbar.Pt(), rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DiTopPtGeo, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.DiTopPtGeo, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.DiTopPtRel, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.DiTopPtRel, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , GlobalVars.HTj, rweight*0.5);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , GlobalVars.HTj, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , GlobalVars.HTj + GlobalVars.met, rweight*0.5);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , GlobalVars.HTj + GlobalVars.met, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , GlobalVars.MJVis, rweight*0.5);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , GlobalVars.MJVis, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , GlobalVars.met, rweight*0.5);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , GlobalVars.met, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , GlobalVars.SumMJ, rweight*0.5);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , GlobalVars.SumMJ, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.pseudotophadron.Pt(), rweight*0.25);
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.pseudotoplepton.Pt(), rweight*0.25); // add1
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.pseudotophadron.Pt(), rweight*0.25); // add2
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.pseudotoplepton.Pt(), rweight*0.25); // add2
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.TopPtRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.TopPtRel[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.TopPtRel[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.Pout[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.Yboost, rweight*0.5);
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.Yboost, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , GlobalVars.jApla, rweight*0.5);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , GlobalVars.jApla, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , GlobalVars.jSpher, rweight*0.5);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , GlobalVars.jSpher, rweight*0.5); // add2
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutGeo" + reptag, pseudo.DiTopPoutGeo[0] , pseudo.TopPtRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutGeo" + reptag, pseudo.DiTopPoutGeo[1] , pseudo.TopPtRel[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutGeo" + reptag, pseudo.DiTopPoutGeo[0] , pseudo.TopPtRel[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutGeo" + reptag, pseudo.DiTopPoutGeo[1] , pseudo.DiTopPoutGeo[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[0] , pseudo.DiTopPoutGeo[0], rweight*0.25);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[1] , pseudo.DiTopPoutGeo[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[0] , pseudo.DiTopPoutGeo[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[1] , pseudo.DiTopPoutRel[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[0] , pseudo.TopPtRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[1] , pseudo.TopPtRel[1], rweight*0.25); // add4
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[0] , pseudo.TopPtRel[1], rweight*0.25); // add5
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPoutRel" + reptag, pseudo.DiTopPoutRel[1] , pseudo.DiTopPoutRel[0], rweight*0.25); // add6
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.Chittbar, rweight);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.CosThetaStar[0], rweight*0.5);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.CosThetaStar[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.deltattbar, rweight);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.pseudotophadron.Pt(), rweight*0.5);
            hold -> FillTH2D(m_level + "TopPt" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.pseudotoplepton.Pt(), rweight*0.5); // add1
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.Yboost, rweight);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , GlobalVars.jApla, rweight);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , GlobalVars.jSpher, rweight);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "DiTopPtGeo" + reptag, pseudo.DiTopPtGeo , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPtGeo" + reptag, pseudo.DiTopPtGeo , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPtGeo" + reptag, pseudo.DiTopPtGeo , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPtGeo" + reptag, pseudo.DiTopPtGeo , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPtGeo" + reptag, pseudo.DiTopPtGeo , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPtGeo" + reptag, pseudo.DiTopPtGeo , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPtGeo" + reptag, pseudo.DiTopPtGeo , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "DiTopPtRel" + reptag, pseudo.DiTopPtRel , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "HTj" + reptag, GlobalVars.HTj , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "HTjPlusMet" + reptag, GlobalVars.HTj + GlobalVars.met , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "MJVis" + reptag, GlobalVars.MJVis , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Met" + reptag, GlobalVars.met , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "SumMJ" + reptag, GlobalVars.SumMJ , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.Chittbar, rweight*0.5);
            hold -> FillTH2D(m_level + "Chittbar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.Chittbar, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.CosThetaStar[0], rweight*0.25);
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.CosThetaStar[0], rweight*0.25); // add1
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.CosThetaStar[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "CosThetaStar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.CosThetaStar[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.deltattbar, rweight*0.5);
            hold -> FillTH2D(m_level + "Delta" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.deltattbar, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.DiTopMassGeo, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.DiTopMassGeo, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.DiTopPoutGeo[0], rweight*0.25);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.DiTopPoutGeo[0], rweight*0.25); // add1
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.DiTopPoutGeo[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.DiTopPoutGeo[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.DiTopPoutRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.DiTopPoutRel[0], rweight*0.25); // add1
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.DiTopPoutRel[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.DiTopPoutRel[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.DiTopPtGeo, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.DiTopPtGeo, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.DiTopPtRel, rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.DiTopPtRel, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , GlobalVars.HTj, rweight*0.5);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , GlobalVars.HTj, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , GlobalVars.HTj + GlobalVars.met, rweight*0.5);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , GlobalVars.HTj + GlobalVars.met, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , GlobalVars.MJVis, rweight*0.5);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , GlobalVars.MJVis, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , GlobalVars.met, rweight*0.5);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , GlobalVars.met, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , GlobalVars.SumMJ, rweight*0.5);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , GlobalVars.SumMJ, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.TopPtRel[0], rweight*0.25);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.TopPtRel[0], rweight*0.25); // add1
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.TopPtRel[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.TopPtRel[1], rweight*0.25); // add3
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.Yboost, rweight*0.5);
            hold -> FillTH2D(m_level + "Yboost" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.Yboost, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , GlobalVars.jApla, rweight*0.5);
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , GlobalVars.jApla, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , GlobalVars.jSpher, rweight*0.5);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , GlobalVars.jSpher, rweight*0.5); // add1
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , GlobalVars.jApla, rweight);
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , GlobalVars.jSpher, rweight);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , pseudo.TopPtRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "jApla" + reptag, GlobalVars.jApla , GlobalVars.jSpher, rweight);
            hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.DiTopMassGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.DiTopPoutGeo[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.DiTopPoutRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
            hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.DiTopPtGeo, rweight);
            hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.DiTopPtRel, rweight);
            hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , GlobalVars.HTj, rweight);
            hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , GlobalVars.HTj + GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , GlobalVars.MJVis, rweight);
            hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , GlobalVars.met, rweight);
            hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , GlobalVars.SumMJ, rweight);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.TopPtRel[0], rweight*0.5);
            hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "jSpher" + reptag, GlobalVars.jSpher , pseudo.TopPtRel[1], rweight*0.5); // add3


            // Rttbar filling code:
            if (pseudo.Rttbar > 0.) {
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "Chittbar" + reptag, pseudo.Chittbar , pseudo.Rttbar, rweight);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[0] , pseudo.Rttbar, rweight*0.5);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "CosThetaStar" + reptag, pseudo.CosThetaStar[1] , pseudo.Rttbar, rweight*0.5); // add2
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "Delta" + reptag, pseudo.deltattbar , pseudo.Rttbar, rweight);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "DiTopDeltaPhi" + reptag, pseudo.DeltaPhi , pseudo.Rttbar, rweight);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "DiTopMass" + reptag, pseudo.pseudottbar.M() , pseudo.Rttbar, rweight);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[0] , pseudo.Rttbar, rweight*0.5);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "DiTopPout" + reptag, pseudo.Pout[1] , pseudo.Rttbar, rweight*0.5); // add2
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "DiTopPt" + reptag, pseudo.pseudottbar.Pt() , pseudo.Rttbar, rweight);
                hold -> FillTH2D(m_level + "DiTopMassGeo" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.DiTopMassGeo, rweight);
                hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.DiTopPoutGeo[0], rweight*0.5);
                hold -> FillTH2D(m_level + "DiTopPoutGeo" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.DiTopPoutGeo[1], rweight*0.5); // add3
                hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.DiTopPoutRel[0], rweight*0.5);
                hold -> FillTH2D(m_level + "DiTopPoutRel" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.DiTopPoutRel[1], rweight*0.5); // add3
                hold -> FillTH2D(m_level + "DiTopPtGeo" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.DiTopPtGeo, rweight);
                hold -> FillTH2D(m_level + "DiTopPtRel" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.DiTopPtRel, rweight);
                hold -> FillTH2D(m_level + "HTj" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , GlobalVars.HTj, rweight);
                hold -> FillTH2D(m_level + "HTjPlusMet" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , GlobalVars.HTj + GlobalVars.met, rweight);
                hold -> FillTH2D(m_level + "MJVis" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , GlobalVars.MJVis, rweight);
                hold -> FillTH2D(m_level + "Met" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , GlobalVars.met, rweight);
                hold -> FillTH2D(m_level + "SumMJ" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , GlobalVars.SumMJ, rweight);
                hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.TopPtRel[0], rweight*0.5);
                hold -> FillTH2D(m_level + "TopPtRel" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , pseudo.TopPtRel[1], rweight*0.5); // add3
                hold -> FillTH2D(m_level + "jApla" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , GlobalVars.jApla, rweight);
                hold -> FillTH2D(m_level + "jSpher" + "Vs" + m_level + "Rttbar" + reptag, pseudo.Rttbar , GlobalVars.jSpher, rweight);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotophadron.Pt() , pseudo.Rttbar, rweight*0.5);
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "TopPt" + reptag, pseudo.pseudotoplepton.Pt() , pseudo.Rttbar, rweight*0.5); // add1
                hold -> FillTH2D(m_level + "Rttbar" + "Vs" + m_level + "Yboost" + reptag, pseudo.Yboost , pseudo.Rttbar, rweight);
            }
        } // replicas or just the standard case

    }

    // special relative
    hold -> FillTH1D(addtag + tname + "PtRel", pseudo.TopPtRel[0], weight);
    hold -> FillTH1D(addtag + tname + "PtRel", pseudo.TopPtRel[1], weight);
    // Top1, Top2:
    tname.ReplaceAll("Top", "Top1");
    hold -> FillTH1D(addtag + tname + "PtRel", pseudo.TopPtRel[0], weight);
    tname.ReplaceAll("Top1", "Top2");
    hold -> FillTH1D(addtag + tname + "PtRel", pseudo.TopPtRel[1], weight);

}

// __________________________________________________


void HistoMaker::FillGlobalHistos(TString cutname, kGlobalVars GlobalVars, double weight)
{
    HistoHolder* hold = 0;
    try {
      hold = m_holders[cutname];
    }  catch (exception& e) {
      cout << "ERROR in HistoMaker::MakeLjetHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
      return;
     }

    hold -> FillTH1D(m_level + "jAplanarity", GlobalVars.jApla, weight);
    hold -> FillTH1D(m_level + "jSphericity", GlobalVars.jSpher, weight);
    hold -> FillTH1D(m_level + "JAplanarity", GlobalVars.JApla, weight);
    hold -> FillTH1D(m_level + "JSphericity", GlobalVars.JSpher, weight);
    hold -> FillTH1D(m_level + "Met", GlobalVars.met, weight);
    hold -> FillTH1D(m_level + "MetPhi", GlobalVars.metphi, weight);
    hold -> FillTH1D(m_level + "HTj", GlobalVars.HTj, weight);
    hold -> FillTH1D(m_level + "HTJ", GlobalVars.HTJ, weight);
    hold -> FillTH1D(m_level + "HTjPlusMet", GlobalVars.HTj + GlobalVars.met, weight);
    hold -> FillTH1D(m_level + "HTJPlusMet", GlobalVars.HTJ + GlobalVars.met, weight);
    hold -> FillTH1D(m_level + "SumMJ", GlobalVars.SumMJ, weight);
    hold -> FillTH1D(m_level + "SumMj", GlobalVars.SumMj, weight);
    hold -> FillTH1D(m_level + "MjVis", GlobalVars.MjVis, weight);
    hold -> FillTH1D(m_level + "MJVis", GlobalVars.MJVis, weight);

}

// __________________________________________________

void HistoMaker::FillLjetHistos(TString cutname, kPseudotop pseudo, kGlobalVars GlobalVars, double weight, bool fillReplicas)
{

    // not used in allhad studies

  HistoHolder* hold = 0;
  try {
    hold = m_holders[cutname];
  }  catch (exception& e) {
    cout << "ERROR in HistoMaker::MakeLjetHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
   }

  hold -> FillTH1D(m_level + "_pseudottbarPt", pseudo.pseudottbar.Pt(), weight);
  hold -> FillTH1D(m_level + "_pseudottbarRapidity", pseudo.pseudottbar.Rapidity(), weight);
  hold -> FillTH1D(m_level + "_pseudottbarMass", pseudo.pseudottbar.M(), weight);

  this -> FillTTbarSpecialHistos(cutname, "pseudottbar", pseudo, GlobalVars, weight, fillReplicas);

  hold -> FillTH1D(m_level + "_pseudoleptonPt", pseudo.lepton.Pt(), weight);
  hold -> FillTH1D(m_level + "_pseudoleptonRapidity", pseudo.lepton.Rapidity(), weight);
  hold -> FillTH1D(m_level + "_pseudoleptonMass", pseudo.lepton.M(), weight);

  hold -> FillTH1D(m_level + "_pseudoneutrinoPt", pseudo.neutrino.Pt(), weight);
  hold -> FillTH1D(m_level + "_pseudoneutrinoRapidity", pseudo.neutrino.Rapidity(), weight);
  hold -> FillTH1D(m_level + "_pseudoneutrinoMass", pseudo.neutrino.M(), weight);

  hold -> FillTH1D(m_level + "_pseudoWleptonPt", pseudo.pseudoWlepton.Pt(), weight);
  hold -> FillTH1D(m_level + "_pseudoWleptonRapidity", pseudo.pseudoWlepton.Rapidity(), weight);
  hold -> FillTH1D(m_level + "_pseudoWleptonMass", pseudo.pseudoWlepton.M(), weight);

  hold -> FillTH1D(m_level + "_pseudoWhadronPt", pseudo.pseudoWhadron.Pt(), weight);
  hold -> FillTH1D(m_level + "_pseudoWhadronRapidity", pseudo.pseudoWhadron.Rapidity(), weight);
  hold -> FillTH1D(m_level + "_pseudoWhadronMass",pseudo.pseudoWhadron.M(), weight);

  hold -> FillTH1D(m_level + "_pseudotopleptonPt", pseudo.pseudotoplepton.Pt(), weight);
  hold -> FillTH1D(m_level + "_pseudotopleptonRapidity", pseudo.pseudotoplepton.Rapidity(), weight);
  hold -> FillTH1D(m_level + "_pseudotopleptonMass", pseudo.pseudotoplepton.M(), weight);

  hold -> FillTH1D(m_level + "_pseudotophadronPt", pseudo.pseudotophadron.Pt(), weight);
  hold -> FillTH1D(m_level + "_pseudotophadronRapidity", pseudo.pseudotophadron.Rapidity(), weight);
  hold -> FillTH1D(m_level + "_pseudotophadronMass", pseudo.pseudotophadron.M(), weight);

}

// __________________________________________________

void HistoMaker::FillGlobalMigrations(TString cutname, TString objname, kGlobalVars globalVars_ptcl, kGlobalVars globalVars_det, double weight)
{

    HistoHolder* hold = 0;
    try {
      hold = m_holders[cutname];
    }  catch (exception& e) {
      cout << "ERROR in HistoMaker::FillGlobalMigrations: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
      return;
    }

    TString migratag = "_" + m_level1 + "_" + m_level2;

    hold -> FillTH2D(objname + "Met" + migratag,  globalVars_ptcl.met, globalVars_det.met, weight);
    hold -> FillTH2D(objname + "MetPhi" + migratag,  globalVars_ptcl.metphi, globalVars_det.metphi, weight);
    hold -> FillTH2D(objname + "jApla" + migratag,  globalVars_ptcl.jApla, globalVars_det.jApla, weight);
    hold -> FillTH2D(objname + "JApla" + migratag,  globalVars_ptcl.JApla, globalVars_det.JApla, weight);
    hold -> FillTH2D(objname + "jSpher" + migratag,  globalVars_ptcl.jSpher, globalVars_det.jSpher, weight);
    hold -> FillTH2D(objname + "JSpher" + migratag,  globalVars_ptcl.JSpher, globalVars_det.JSpher, weight);


    }

// __________________________________________________

void HistoMaker::FillSpecialMigrations(TString cutname, TString objname, kPseudotop pseudo_ptcl, kPseudotop pseudo_det, double weight)
{

    HistoHolder* hold = 0;
    try {
      hold = m_holders[cutname];
    }  catch (exception& e) {
      cout << "ERROR in HistoMaker::FillSpecialMigrations: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
      return;
    }

    TString migratag = "_" + m_level1 + "_" + m_level2;
/*
 *  this would be duplication ;-)
    hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_ptcl.pseudottbar.Pt(),       pseudo_det.pseudottbar.Pt(), weight);
    hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_ptcl.pseudottbar.Rapidity(), pseudo_det.pseudottbar.Rapidity(), weight);
    hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_ptcl.pseudottbar.M(),        pseudo_det.pseudottbar.M(), weight);
  */
    hold -> FillTH2D(objname + "Pout" + migratag,     fabs(pseudo_ptcl.Pout[0]),        fabs(pseudo_det.Pout[0]), weight);
    hold -> FillTH2D(objname + "Pout" + migratag,     fabs(pseudo_ptcl.Pout[1]),        fabs(pseudo_det.Pout[1]), weight);
    hold -> FillTH2D(objname + "DeltaPhi" + migratag,     pseudo_ptcl.DeltaPhi,        pseudo_det.DeltaPhi, weight);

    hold -> FillTH2D(objname + "Yboost" + migratag,     pseudo_ptcl.Yboost,        pseudo_det.Yboost, weight);
    if (pseudo_det.pseudotophadron.Pt() > 0 && pseudo_ptcl.pseudotophadron.Pt())
        hold -> FillTH2D(objname + "Rttbar" + migratag, pseudo_ptcl.pseudotoplepton.Pt() / pseudo_ptcl.pseudotophadron.Pt(), pseudo_det.pseudotoplepton.Pt() / pseudo_det.pseudotophadron.Pt(), weight);

    hold -> FillTH2D(objname + "Chittbar" + migratag,     pseudo_ptcl.Chittbar,        pseudo_det.Chittbar, weight);

    hold -> FillTH2D(objname + "Delta" + migratag,     pseudo_ptcl.deltattbar,        pseudo_det.deltattbar, weight);
    //hold -> FillTH2D(objname + "CosThetaStar" + migratag,     fabs(pseudo_ptcl.CosThetaStar[0]),        fabs(pseudo_det.CosThetaStar[0]), weight);
    // beware of swaps, esp. in multiple hadron topologies...
    hold -> FillTH2D(objname + "CosThetaStar" + migratag,     fabs(pseudo_ptcl.CosThetaStar[0]),        fabs(pseudo_det.CosThetaStar[0]), weight);
    hold -> FillTH2D(objname + "CosThetaStar" + migratag,     fabs(pseudo_ptcl.CosThetaStar[1]),        fabs(pseudo_det.CosThetaStar[1]), weight);


}
// __________________________________________________

void HistoMaker::FillLjetMigrations(TString cutname, kPseudotop pseudo_ptcl, kPseudotop pseudo_det, double weight)
{

  HistoHolder* hold = 0;
  try {
    hold = m_holders[cutname];
  }  catch (exception& e) {
    cout << "ERROR in HistoMaker::FillLjetMigrations: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
  }
  
  TString migratag = "_" + m_level1 + "_" + m_level2;
  
  TString objname = "DiscriminantDmW";
  hold -> FillTH2D(objname + migratag,       
		   pseudo_ptcl.DmW / fabs(pseudo_ptcl.DmW) * pow(fabs(pseudo_ptcl.DmW), 1/6.),
		   pseudo_det.DmW / fabs(pseudo_det.DmW) * pow(fabs(pseudo_det.DmW), 1/6.),
		   weight);
  objname = "DiscriminantDmt";
  hold -> FillTH2D(objname + migratag,       
		   pseudo_ptcl.Dmt / fabs(pseudo_ptcl.Dmt) * pow(fabs(pseudo_ptcl.Dmt), 1/6.),
		   pseudo_det.Dmt / fabs(pseudo_det.Dmt) * pow(fabs(pseudo_det.Dmt), 1/6.),
		   weight);


  objname = "pseudottbar";
  this -> FillSpecialMigrations(cutname, objname, pseudo_ptcl, pseudo_det, weight);

  objname = "pseudotophadron";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_ptcl.pseudotophadron.Pt(),       pseudo_det.pseudotophadron.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_ptcl.pseudotophadron.Rapidity(), pseudo_det.pseudotophadron.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_ptcl.pseudotophadron.M(),        pseudo_det.pseudotophadron.M(), weight);

  objname = "pseudotoplepton";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_ptcl.pseudotoplepton.Pt(),       pseudo_det.pseudotoplepton.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_ptcl.pseudotoplepton.Rapidity(), pseudo_det.pseudotoplepton.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_ptcl.pseudotoplepton.M(),        pseudo_det.pseudotoplepton.M(), weight);

  objname = "pseudotopleptonVspseudotophadron";
  hold -> FillTH2D(objname + "Mass_" + m_level2, pseudo_det.pseudotophadron.M(), pseudo_det.pseudotoplepton.M(), weight);
  hold -> FillTH2D(objname + "Mass_" + m_level1, pseudo_ptcl.pseudotophadron.M(), pseudo_ptcl.pseudotoplepton.M(), weight);
  hold -> FillTH2D(objname + "Pt_" + m_level2, pseudo_det.pseudotophadron.Pt(), pseudo_det.pseudotoplepton.Pt(), weight);
  hold -> FillTH2D(objname + "Pt_" + m_level1, pseudo_ptcl.pseudotophadron.Pt(), pseudo_ptcl.pseudotoplepton.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity_" + m_level2, pseudo_det.pseudotophadron.Rapidity(), pseudo_det.pseudotoplepton.Rapidity(), weight);
  hold -> FillTH2D(objname + "Rapidity_" + m_level1, pseudo_ptcl.pseudotophadron.Rapidity(), pseudo_ptcl.pseudotoplepton.Rapidity(), weight);
  
  objname = "pseudoWhadron";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_ptcl.pseudoWhadron.Pt(),       pseudo_det.pseudoWhadron.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_ptcl.pseudoWhadron.Rapidity(), pseudo_det.pseudoWhadron.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_ptcl.pseudoWhadron.M(),        pseudo_det.pseudoWhadron.M(), weight);

  objname = "pseudoWlepton";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_ptcl.pseudoWlepton.Pt(),       pseudo_det.pseudoWlepton.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_ptcl.pseudoWlepton.Rapidity(), pseudo_det.pseudoWlepton.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_ptcl.pseudoWlepton.M(),        pseudo_det.pseudoWlepton.M(), weight);

  objname = "pseudolepton";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_ptcl.lepton.Pt(),       pseudo_det.lepton.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_ptcl.lepton.Rapidity(), pseudo_det.lepton.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_ptcl.lepton.M(),        pseudo_det.lepton.M(), weight);
   
  objname = "pseudoneutrino";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_ptcl.neutrino.Pt(),       pseudo_det.neutrino.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_ptcl.neutrino.Rapidity(), pseudo_det.neutrino.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_ptcl.neutrino.M(),        pseudo_det.neutrino.M(), weight);
   

  return;
  
}


// __________________________________________________

void HistoMaker::FillLjetMigrationsParton(TString cutname, kPseudotop pseudo_parton, kPseudotop pseudo, double weight)
{

    // not maintained!!! jk 2020

  HistoHolder* hold = 0;
  try {
    hold = m_holders[cutname];
  }  catch (exception& e) {
    cout << "ERROR in HistoMaker::FillLjetMigrationsParton: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
  }
  
  TString migratag = "_" + m_level1 + "_" + m_level2;

  // double pseudo.DmW / fabs(pseudo.DmW);

  TString objname = "pseudottbar";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_parton.pseudottbar.Pt(),       pseudo.pseudottbar.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_parton.pseudottbar.Rapidity(), pseudo.pseudottbar.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudo_parton.pseudottbar.M(),        pseudo.pseudottbar.M(), weight);

  if (pseudo_parton.pseudotoplepton.DeltaR(pseudo.pseudotoplepton) < pseudo_parton.pseudotoplepton.DeltaR(pseudo.pseudotophadron)) {
    objname = "pseudotoplepton";
    hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_parton.pseudotoplepton.Pt(),       pseudo.pseudotoplepton.Pt(), weight);
    hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_parton.pseudotoplepton.Rapidity(), pseudo.pseudotoplepton.Rapidity(), weight);
    objname = "pseudotophadron";
    hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_parton.pseudotophadron.Pt(),       pseudo.pseudotophadron.Pt(), weight);
    hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_parton.pseudotophadron.Rapidity(), pseudo.pseudotophadron.Rapidity(), weight);
    objname = "pseudottbar";
    hold -> FillTH2D(objname + "Pout" + migratag,     pseudo_parton.Pout[0],        pseudo.Pout[0], weight);
    hold -> FillTH2D(objname + "Pout" + migratag,     pseudo_parton.Pout[1],        pseudo.Pout[1], weight);
  } else {
    // switch partonic top roles:
    objname = "pseudotoplepton";
    hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_parton.pseudotophadron.Pt(),       pseudo.pseudotoplepton.Pt(), weight);
    hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_parton.pseudotophadron.Rapidity(), pseudo.pseudotoplepton.Rapidity(), weight);
    objname = "pseudotophadron";
    hold -> FillTH2D(objname + "Pt" + migratag,       pseudo_parton.pseudotoplepton.Pt(),       pseudo.pseudotophadron.Pt(), weight);
    hold -> FillTH2D(objname + "Rapidity" + migratag, pseudo_parton.pseudotoplepton.Rapidity(), pseudo.pseudotophadron.Rapidity(), weight);
    objname = "pseudottbar";
    hold -> FillTH2D(objname + "Pout" + migratag,     pseudo_parton.Pout[1],        pseudo.Pout[0], weight);
    hold -> FillTH2D(objname + "Pout" + migratag,     pseudo_parton.Pout[0],        pseudo.Pout[1], weight);
  }
  
  // This are global variales and the lep/had identification does not matter: 
  objname = "pseudottbar";
  hold -> FillTH2D(objname + "DeltaPhi" + migratag,     pseudo_parton.DeltaPhi,        pseudo.DeltaPhi, weight);

  hold -> FillTH2D(objname + "Chittbar" + migratag,     pseudo_parton.Chittbar,        pseudo.Chittbar, weight);
  hold -> FillTH2D(objname + "Yboost" + migratag,     pseudo_parton.Yboost,        pseudo.Yboost, weight);
  
  hold -> FillTH2D(objname + "Delta" + migratag,     pseudo_parton.deltattbar,        pseudo.deltattbar, weight);
  hold -> FillTH2D(objname + "CosThetaStar" + migratag,     fabs(pseudo_parton.CosThetaStar[0]),        fabs(pseudo.CosThetaStar[0]), weight);
  hold -> FillTH2D(objname + "CosThetaStar" + migratag,     fabs(pseudo_parton.CosThetaStar[1]),        fabs(pseudo.CosThetaStar[1]), weight);

  /*  
  objname = "pseudotophadron";
  hold -> FillTH2D(objname + "Pt" + migratag,       pseudotophadron.Pt(),       pseudo.pseudotophadron.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, pseudotophadron.Rapidity(), pseudo.pseudotophadron.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     pseudotophadron.M(),        pseudo.pseudotophadron.M(), weight);

  objname = "pseudotoplepton";
  hold -> FillTH2D(objname + "Pt" + migratag,       partontoplepton.Pt(),       pseudo.pseudotoplepton.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, partontoplepton.Rapidity(), pseudo.pseudotoplepton.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     partontoplepton.M(),        pseudo.pseudotoplepton.M(), weight);
  */
  /*
  objname = "pseudoWhadron";
  hold -> FillTH2D(objname + "Pt" + migratag,       partonWhadron.Pt(),       pseudo.pseudoWhadron.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, partonWhadron.Rapidity(), pseudo.pseudoWhadron.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     partonWhadron.M(),        pseudo.pseudoWhadron.M(), weight);

  objname = "pseudoWlepton";
  hold -> FillTH2D(objname + "Pt" + migratag,       partonWlepton.Pt(),       pseudo.pseudoWlepton.Pt(), weight);
  hold -> FillTH2D(objname + "Rapidity" + migratag, partonWlepton.Rapidity(), pseudo.pseudoWlepton.Rapidity(), weight);
  hold -> FillTH2D(objname + "Mass" + migratag,     partonWlepton.M(),        pseudo.pseudoWlepton.M(), weight);
  */



  return;
  
}


// __________________________________________________

void HistoMaker::MakeJetHistos(TString JetType, HistoHolder *hold, bool MakeSubs, bool makeReplicas) // whether to make jet substructure histos
{
    if (!gDirectory->GetDirectory(hold->GetName())) {
        gDirectory->mkdir(hold->GetName());
    }
    gDirectory->cd(hold->GetName());
	
    int nMaxJetsToFill = ComputeMaxJetsToFill(JetType);


    bool IsSingleObject = JetType.Contains("Lepton") || JetType.Contains("lepton") ||
                              JetType.Contains("Neutrino") || JetType.Contains("MET") ||
                              JetType.Contains("neutrino");

    bool IsTop = (!JetType.Contains("DiTop")) && (JetType.Contains("Top") );
    bool IsDiTop = JetType.Contains("DiTop");

    TString jname = m_level + JetType;
    TString jtag = m_level + " " + JetType;
	TString GeV = " [GeV]";
	
    if (!IsSingleObject) {
      int nNBins = 16;
	  double nmin = -0.5;
      double nmax = 15.5;
      // cout << "making " << jname.Data() << "N" << endl;
      hold -> AddTH1D(jname + "N", jname + "N;N_{" + jname + "}", nNBins, nmin, nmax);
    }

    if (JetType.Contains("ttbar") || JetType.Contains("DiTop") ) {
       this -> MakeSpecialHistos(hold, jname, jtag, makeReplicas);
    }

	int nPtBins = 300;
	double pTmin = 0;
	double pTmax = 1500;
    if (JetType.Contains("Dijet") || JetType.Contains("DiTop")) {
        pTmax = 1000;
    } else if (JetType.Contains("Top")) {
        pTmax = 1200;
        pTmin = 0;
    }
    if (JetType.Contains("W")) {
        pTmax = 1000;
    }
    TString topotag = "_" + hold -> GetName();
    // cout << "topotag: " << topotag.Data() << endl;
    // for chosing some existing binning, mostly for parton histos:
    if (topotag == "_NoCuts") {
      topotag = "_AnySel";
    }

    TString repdir = "replicas";
    if (makeReplicas && !gDirectory->GetDirectory(repdir)) {
        gDirectory->mkdir(repdir);
    }
    if (IsTop) {
            hold -> AddTH1D(jname + "Pt", jname + "Pt;" + jtag + " p_{T}" + GeV, m_nPhysObjBins["TopPt" + topotag], m_PhysObjBins["TopPt"  + topotag]);
            if (makeReplicas) {
                for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
                    TString reptag = Form("_rep%i",irep);
                    gDirectory->cd(repdir);
                    hold -> AddTH1D(jname + "Pt" + reptag, jname + "Pt;" + jtag + " p_{T}" + GeV, m_nPhysObjBins["TopPt" + topotag], m_PhysObjBins["TopPt"  + topotag]);
                    gDirectory->cd("../");
                }
            }
    } else if (IsDiTop ) {
        hold -> AddTH1D(jname + "Pt", jname + "Pt;" + jtag + " p_{T}" + GeV, m_nPhysObjBins["DiTopPt"], m_PhysObjBins["DiTopPt"]);
        if (makeReplicas) {
            for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
                TString reptag = Form("_rep%i",irep);
                gDirectory->cd(repdir);
                hold -> AddTH1D(jname + "Pt" + reptag, jname + "Pt;" + jtag + " p_{T}" + GeV, m_nPhysObjBins["DiTopPt"], m_PhysObjBins["DiTopPt"]);
                gDirectory->cd("../");
            }
        }
        TString repdir = "replicas";
        if (makeReplicas && !gDirectory->GetDirectory(repdir)) {
            gDirectory->mkdir(repdir);
        }

        // correlations between variables:
        // special correlation histos!
        // 6.8.2020
        // c.f. also TestGenCorrCode.cxx
        // https://stackoverflow.com/questions/6963894/how-to-use-range-based-for-loop-with-stdmap
        // make and fill this only for Detector level!!!
        if (m_level.Contains("Detector")) {
            int nvars = 0;
            for (auto& [var1, varProp1]: m_vars) {
                for (auto& [var2, varProp2]: m_vars) {
                    if (var1 == var2)
                        continue;
                    if (varProp1.id > varProp2.id)
                        continue;
                    /* jk 30.4.2021
                    if (var1.Contains("j") && var2.Contains("J"))
                        continue;
                    if (var1.Contains("J") && var2.Contains("j"))
                        continue;
                    */
                    nvars++;
                    // fragile, but works to use for jApla the Apla binning etc
                    // jk 30.11.2020
                    TString binvar1 = var1;
                    binvar1.ReplaceAll("J","");
                    binvar1.ReplaceAll("j","");
                    binvar1.ReplaceAll("PlusMet","");
                    binvar1 += topotag;
                    TString binvar2 = var2;
                    binvar2.ReplaceAll("J","");
                    binvar2.ReplaceAll("j","");
                    binvar2.ReplaceAll("PlusMet","");
                    binvar2 += topotag;

                    TString name = m_level + var2 + "Vs" + m_level + var1;
                    //cout << "Making special " << name.Data() << " in topo " << topotag << endl;
                    //cout << "   nbins: based on bv1=" << binvar1.Data() << ",binvar2=" << binvar2.Data() << ": "
                    //     << m_nPhysObjBins[binvar1] << "," << m_nPhysObjBins[binvar2] << endl;
                    hold -> AddTH2D(name,
                                    name + ";" + varProp1.title + varProp1.unit + ";" + varProp2.title  + varProp2.unit + ";",
                                    m_nPhysObjBins[binvar1], m_PhysObjBins[binvar1],
                                    m_nPhysObjBins[binvar2], m_PhysObjBins[binvar2]);

                    if (makeReplicas) {
                        TString repdir = "replicas";
                        gDirectory->cd(repdir);
                        for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
                            TString reptag = Form("_rep%i",irep);
                            hold -> AddTH2D(name + reptag,
                                            name + ";" + varProp1.title + varProp1.unit + ";" + varProp2.title  + varProp2.unit + ";",
                                            m_nPhysObjBins[binvar1], m_PhysObjBins[binvar1],
                                            m_nPhysObjBins[binvar2], m_PhysObjBins[binvar2]);
                        }
                        gDirectory->cd("../");
                    } // make replicas


                } // var2
            } // var1
            cout << "Will study total of " << nvars << " 2D histos between variables!" << endl;
        } // Detector

    } else {
        hold -> AddTH1D(jname + "Pt", jname + "Pt;" + jtag + " p_{T}" + GeV, nPtBins, pTmin, pTmax);
        hold -> AddTH1D("LeadDi" + jname + "Mass", "LeadDi" + jname + "Mass;LeadinDi" + jtag + " Mass" + GeV, 250., 0., 1000.);
        hold -> AddTH1D("LeadDi" + jname + "dR", "LeadDi" + jname + "dR;LeadingDi" + jtag + " Mass" + GeV, 100, 0, 5.);
    }

    int nRapidityBins = 250;
    double rapmin = -2.5;
    double rapmax = 2.5;

    if (IsTop || IsDiTop) {
        hold -> AddTH1D(jname + "Rapidity", jname + "Rapidity;" + jtag + " y", m_nPhysObjBins["Rapidity"], m_PhysObjBins["Rapidity"]);
        if (makeReplicas) {
            for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
                TString reptag = Form("_rep%i",irep);
                gDirectory->cd(repdir);
                hold -> AddTH1D(jname + "Rapidity" + reptag, jname + "Rapidity;" + jtag + " y", m_nPhysObjBins["Rapidity"], m_PhysObjBins["Rapidity"]);
                gDirectory->cd("../");
            }
        }
    } else {
        hold -> AddTH1D(jname + "Rapidity", jname + "Rapidity;" + jtag + " y", nRapidityBins, rapmin, rapmax);
    }

	int nMassBins = 300;
	double Massmin = 0;
    double Massmax = 300;

    if (JetType.Contains("Dijet") || JetType.Contains("DiTop")) {
        Massmax = 3000;
        nMassBins = 300;
    } else if (JetType.Contains("Top")) {
        Massmax = 250;
        Massmin = 100;
        nMassBins = 150;

    } else if (JetType.Contains("W")) {
        Massmax = 150;
        Massmin = 50;
        nMassBins = 200;

    }

    if (IsDiTop ) {
        hold -> AddTH1D(jname + "Mass", jname + "Mass;" + jtag + " m" + GeV, m_nPhysObjBins["DiTopMass" + topotag], m_PhysObjBins["DiTopMass" + topotag]);
        if (makeReplicas) {
            for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
                TString reptag = Form("_rep%i",irep);
                gDirectory->cd(repdir);
                hold -> AddTH1D(jname + "Mass" + reptag, jname + "Mass;" + jtag + " m" + GeV, m_nPhysObjBins["DiTopMass" + topotag], m_PhysObjBins["DiTopMass" + topotag]);
                gDirectory->cd("../");


            }
        }
    } else {
        hold -> AddTH1D(jname + "Mass", jname + "Mass;" + jtag + " m" + GeV, nMassBins, Massmin, Massmax);
    }
    if (!IsSingleObject) {
        for (int ij = 1; ij <= ComputeMaxJetsToFill(JetType); ++ij) {
            TString sij = Form("%i", ij);
            if (IsTop) {
                hold -> AddTH1D(jname + sij + "Pt", jname + sij + "Pt;" + jtag + sij + " p_{T}" + GeV, m_nPhysObjBins["TopPt" + topotag], m_PhysObjBins["TopPt" + topotag]);
            } else if (IsDiTop ) {
                hold -> AddTH1D(jname + sij + "Pt", jname + sij + "Pt;" + jtag + sij + " p_{T}" + GeV, m_nPhysObjBins["DiTopPt"], m_PhysObjBins["DiTopPt"]);
            } else {
                hold -> AddTH1D(jname + sij + "Pt", jname + sij + "Pt;" + jtag + sij + " p_{T}" + GeV, nPtBins, pTmin, pTmax);
            }
            if (IsDiTop || IsTop) {
                hold -> AddTH1D(jname + sij + "Rapidity", jname + sij + "Rapidity;" + jtag + sij + " y", m_nPhysObjBins["Rapidity"], m_PhysObjBins["Rapidity"]);
            } else {
                hold -> AddTH1D(jname + sij + "Rapidity", jname + sij + "Rapidity;" + jtag + sij + " y", nRapidityBins, rapmin, rapmax);
            }
            // eventually modify this for 4t study and additional DiTop mass pair!!! jk 5.6.2020
            hold -> AddTH1D(jname + sij + "Mass", jname + sij + "Mass;" + jtag + sij + " m" + GeV, nMassBins, Massmin, Massmax);
        }
    }


    if (MakeSubs && !IsSingleObject) {
        hold -> AddTH1D(jname + "Tau32", jname + "Tau32;" + jtag + " #tau_{32}", m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH1D(jname + "Tau21", jname + "Tau21;" + jtag + " #tau_{21}", m_nTauBins, m_Taumin, m_Taumax);

        hold -> AddTH2D(jname + "Tau32VsMass", jname + "Tau32VsMass;" + jtag + " Mass;" + jtag + " #tau_{32}",
                        nMassBins, Massmin, Massmax, m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(jname + "Tau21VsMass", jname + "Tau21VsMass;" + jtag + " Mass;" + jtag + " #tau_{21}",
                        nMassBins, Massmin, Massmax, m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(jname + "Tau21VsTau32", jname + "Tau21VsTau32;" + jtag + " #tau_{32};" + jtag + " #tau_{21}",
                        m_nTauBins, m_Taumin, m_Taumax, m_nTauBins, m_Taumin, m_Taumax);
        // 2.7.2021:
        hold -> AddTH2D(jname + "Tau32VsPt", jname + "Tau32VsPt;" + jtag + " p_{T};" + jtag + " #tau_{32}",
                        nPtBins, pTmin, pTmax, m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(jname + "Tau21VsPt", jname + "Tau21VsPt;" + jtag + " p_{T};" + jtag + " #tau_{21}",
                        nPtBins, pTmin, pTmax, m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(jname + "PtVsMass", jname + "PtVsMass;" + jtag + " Mass;" + jtag + " p_{T}",
                        nMassBins, Massmin, Massmax, nPtBins, pTmin, pTmax);
        hold -> AddTH2D(jname + "PtVsMass", jname + "PtVsMass;" + jtag + " Mass;" + jtag + " p_{T}",
                        nMassBins, Massmin, Massmax, nPtBins, pTmin, pTmax);


        for (int ij = 1; ij <= nMaxJetsToFill; ++ij) {
            TString sij = Form("%i", ij);
            hold -> AddTH1D(jname + sij + "Tau32", jname + sij + "Tau32;" + jtag + sij + " #tau_{32}", m_nTauBins, m_Taumin, m_Taumax);
            hold -> AddTH1D(jname + sij + "Tau21", jname + sij + "Tau21;" + jtag + sij + " #tau_{21}", m_nTauBins, m_Taumin, m_Taumax);

        }
    }
    gDirectory->cd("../");
}


// __________________________________________________

void HistoMaker::MakeDijetHistos(TString cutname)
{
   auto iter = m_holders.find(cutname);
   if (iter != m_holders.end())
	 this -> MakeJetHistos("Dijet", iter->second, false);	
   else
     cerr << "ERROR HistoMaker::MakeDijetHistos: making dijet histos!" << endl;
}


// __________________________________________________

void HistoMaker::MakeJESHistos(TString cutname, bool MakeClosure, bool MakeJES)
{

  HistoHolder* hold = 0;
  try {
    hold = m_holders[cutname];
  }  catch (exception& e) {
    cout << "ERROR in HistoMaker::MakeJESHisto: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
  }

  gDirectory->mkdir("JES");
  gDirectory->cd("JES");

  TString jtags[] = {"LJets", "Jets"};
  for (auto jtag : jtags) {
      TString corrs[] = {"RPt", "RE", "RPtClosure", "REClosure",
                        //"DEta", "DPhi", "DEtaClosure", "DPhiClosure"
                        };
      std::map<TString,TString> zlabels;
      zlabels["RPt"] = "R_{pT}";
      zlabels["RE"] = "R_E";
      zlabels["RPtClosure"] = "R_{pT} Closure";
      zlabels["REClosure"] = "R_{E} Closure";

      zlabels["DEta"] = "#Delta#eta";
      zlabels["DEtaClosure"] = "#Delta#eta Closure";
      zlabels["DPhi"] = "#Delta#phi";
      zlabels["DPhiClosure"] = "#Delta#phi Closure";

      int nz = 400;
      int neta = 100;
      double zzmin = 0.;
      double zzmax = 2.;
      for (auto corr : corrs) {

          if (corr.Contains("Closure") && ! MakeClosure)
              continue;
          if (!corr.Contains("Closure") && ! MakeJES)
              continue;

          double zmin = zzmin;
          double zmax = zzmax;

          if (corr.Contains("Eta")) {
              zmin = -0.3;
              zmax = 0.3;
          }
          if (corr.Contains("Phi")) {
              zmin = 0.;
              zmax = 0.25;
          }

          double etam = 2.5;
          double ptmin = 20.;
          double ptmax = 1020.;
          int npt = 200;
          int ne = 300;
          double emin = 20.;
          double emax = 1520.;
          if (jtag.Contains("LJets")) {
              etam = 2.;
              ptmin = 60;
              ptmax = 1060;
              emin = 60.;
              emax = 1560.;
          }
          hold -> AddTH2D("JES" + corr + jtag + "Pt", "JES" + corr + jtag + "Pt;p_{t};"+zlabels[corr], npt, ptmin, ptmax, nz, zmin, zmax);
          hold -> AddTH2D("JES" + corr + jtag + "E", "JES" + corr + jtag + "E;E;"+zlabels[corr], ne, emin, emax, nz, zmin, zmax);
          hold -> AddTH2D("JES" + corr + jtag + "Eta", "JES" + corr + jtag + "Eta;#eta;"+zlabels[corr], neta, -etam, etam, nz, zmin, zmax);

          // 3D histo of eta pt JES
          if (corr.Contains("R")) {
              hold -> AddTH3D("JES" + corr + jtag + "EtaPt", "JES" + corr + jtag + "EtaPt;#eta;p_{T};"+zlabels[corr], neta, -etam, etam, npt, ptmin, ptmax, nz, zmin, zmax);
              hold -> AddTH3D("JES" + corr + jtag + "EtaE", "JES" + corr + jtag + "EtaE;#eta;E;"+zlabels[corr], neta, -etam, etam, ne, emin, emax, nz, zmin, zmax);
          }
      }
  }
  gDirectory->cd("../");

}

// __________________________________________________

void HistoMaker::MakeControlHistos(bool makeLjet, bool makeDijet)
{
	for (auto iter : m_holders) {
      MakeJetHistos("Jet", iter.second, false);
      MakeJetHistos("bJet", iter.second, false);
      MakeJetHistos("nonbJet", iter.second, false);
	  MakeJetHistos("Lepton", iter.second, false);
	  MakeJetHistos("Neutrino", iter.second, false);
	  if (makeLjet)
	    MakeJetHistos("LJet", iter.second, true);	
	  if (makeDijet)
	    MakeJetHistos("Dijet", iter.second, false);	
    }
}

// __________________________________________________

void HistoMaker::MakeAllJetHistos(bool makeReplicas)
{
    for (auto iter : m_holders) {
      MakeJetHistos("Jet", iter.second, false);
      MakeJetHistos("bJet", iter.second, false);
      MakeJetHistos("nonbJet", iter.second, false);
      MakeJetHistos("LJet", iter.second, true);
      MakeGlobalHistos(iter.second);
      
      if (! iter.second -> JetsOnly()) {
          MakeJetHistos("Top", iter.second, true, makeReplicas);
          MakeJetHistos("W", iter.second, true);
          MakeJetHistos("DiTop", iter.second, false, makeReplicas);
          // MakeJetHistos("TriTop", iter.second, false);
          MakeJetHistos("FourTop", iter.second, false);

          /* this, or use Top1 and Top2?
     * or keep Top1 and Top2 for BS in the order??
     * */
          //MakeJetHistos("LeadingTop", iter.second, true); // , makeReplicas ???
          //MakeJetHistos("SubleadingTop", iter.second, true); // , makeReplicas ???
      } // !jets only

    }
}
// __________________________________________________

void HistoMaker::MakePartonHistos()
{
    for (auto iter : m_holders) {
        MakeJetHistos("Top", iter.second, false);
        MakeJetHistos("W", iter.second, false);
        MakeJetHistos("DiTop", iter.second, false);
	/*
        MakeJetHistos("PartonTop", iter.second, false);
        MakeJetHistos("PartonW", iter.second, false);
        MakeJetHistos("PartonDiTop", iter.second, false);
	*/
    }
}


// __________________________________________________

void HistoMaker::MakeGlobalMigrations(TString cutname, TString objname, TString level1, TString level2)
{

    auto iter = m_holders.find(cutname);
    if (iter != m_holders.end()) {

        auto hold = iter->second;
       /*
        if (!gDirectory->GetDirectory(hold->GetName())) {
            gDirectory->mkdir(hold->GetName());
        }
        gDirectory->cd(hold->GetName());
        */
        m_level1 = level1;
        m_level2 = level2;

        TString migratag = "_" + m_level1 + "_" + m_level2;

        TString ptclname = m_level1;
        TString detname = m_level2;
        TString GeV = " [GeV]";

        int nMetBins = 30;
        double maxMet = 300;
        int nPhiBins = 20;
        double minPhi = -TMath::Pi();
        double maxPhi = TMath::Pi();

        // or: use m_nPhysObjBins["Met"], m_PhysObjBins["Met"]
        hold -> AddTH2D(objname + "Met" + migratag,
                        objname + "Met" + migratag + ";" + ptclname + " MET" + GeV + ";" + detname + " MET" + GeV,
                        nMetBins, 0, maxMet, nMetBins, 0, maxMet );

        hold -> AddTH2D(objname + "MetPhi" + migratag,
                        objname + "MetPhi"  + migratag + ";" + ptclname + " MET #phi;" + detname + " MET #phi",
                        nPhiBins, minPhi, maxPhi, nPhiBins, minPhi, maxPhi);

        hold -> AddTH2D(objname + "jApla" + migratag,
                        objname + "jApla" + migratag + ";" + ptclname + " jets Aplanarity;" + detname + " jets Aplanarity",
                        m_nPhysObjBins["Apla"], m_PhysObjBins["Apla"], m_nPhysObjBins["Apla"], m_PhysObjBins["Apla"]);
        hold -> AddTH2D(objname + "JApla" + migratag,
                        objname + "JApla" + migratag + ";" + ptclname + " LJets Aplanarity;" + detname + " LJets Aplanarity",
                        m_nPhysObjBins["Apla"], m_PhysObjBins["Apla"], m_nPhysObjBins["Apla"], m_PhysObjBins["Apla"]);
        hold -> AddTH2D(objname + "jSpher" + migratag,
                        objname + "jSpher" + migratag + ";" + ptclname + " jets Sphericity;" + detname + " jets Sphericity",
                        m_nPhysObjBins["Spher"], m_PhysObjBins["Spher"], m_nPhysObjBins["Spher"], m_PhysObjBins["Spher"] );
        hold -> AddTH2D(objname + "JSpher" + migratag,
                        objname + "JSpher" + migratag + ";" + ptclname + " LJets Sphericity;" + detname + " LJets Sphericity",
                        m_nPhysObjBins["Spher"], m_PhysObjBins["Spher"], m_nPhysObjBins["Spher"], m_PhysObjBins["Spher"] );
      //  gDirectory->cd("../");

    } // holders
}
// __________________________________________________

void HistoMaker::MakeJetMigrations(TString cutname, TString JetType, TString level1, TString level2,
                                   bool FillSub, int AddMore)
{

   auto iter = m_holders.find(cutname);
   
   if (iter != m_holders.end()) {
	   
    auto hold = iter->second;
    //if (!gDirectory->GetDirectory(hold->GetName())) {
    //    gDirectory->mkdir(hold->GetName());
    //}
    //gDirectory->cd(hold->GetName());
    //cout << "Creating migrations for " << cutname.Data() << " in directory " << hold->GetName() << endl;
    m_level1 = level1;
    m_level2 = level2;

    TString migratag = "_" + m_level1 + "_" + m_level2;

    TString ptclname = m_level1 + " " + JetType;
    TString detname = m_level2 + " " + JetType;
	TString GeV = " [GeV]";

	int nNBins = 11;
	double nmin = -0.5;
	double nmax = 10.5; 

	int nPtBins = 20;
	double pTmin = 0;
	double pTmax = 1500;

	int nRapidityBins = 20;
	double Rapiditymin = -2.5;
	double Rapiditymax = 2.5;

	int nMassBins = 20;
	double Massmin = 0;
	double Massmax = 200;

	if (JetType.Contains("pseudotop")) {
	  nPtBins = 20;
	  pTmin = 0;
	  pTmax = 500;
	}
	if (JetType.Contains("pseudoW")) {
	  pTmax = 500;
	  nPtBins = 20;
	  pTmin = 0;
	  pTmax = 500;
	}

    if (JetType.Contains("ttbar") || JetType.Contains("DiTop")) {
      nMassBins = 68;
      Massmin = 300;
      Massmax = 2000;
      nPtBins = 50;
      pTmin = 0;
      pTmax = 500;
	}
	if (JetType.Contains("pseudotop")) {
	  nMassBins = 20;
	  Massmin = 100;
	  Massmax = 300;
	}
	if (JetType.Contains("pseudoW")) {
	  nMassBins = 20;
	  Massmin = 0;
	  Massmax = 160;
	}

    if (JetType.Contains("ttbar") || JetType.Contains("DiTop")) {
	  
      // make the discriminants migration matrices:
      /*
      int nDBins = 50;
      double dmin = -300.;
      double dmax = 300.;
      */

        /*
      TString obj = "DiscriminantDmW";
      hold -> AddTH2D(obj + migratag, obj + migratag + ";" + ptclname + " Signed |D_{m_{W}}|^{1/6} [GeV];" + detname + " Signed |D_{m_{W}}|^{1/6} [GeV]", nDBins, dmin, dmax, nDBins, dmin, dmax);
      obj = "DiscriminantDmt";
      hold -> AddTH2D(obj + migratag, obj + migratag + ";" + ptclname + " Signed |D_{m_{t}}|^{1/6} [GeV];" + detname + " Signed |D_{m_{t}}|^{1/6} [GeV]", nDBins, 2.75*dmin, 2.75*dmax, nDBins, 2.75*dmin, 2.75*dmax);
      */

      // Pout etc.
      hold -> AddTH2D(JetType + "Pout" + migratag, JetType + "Pout" + migratag + ";" + ptclname + " |p_{out}|" + GeV + ";" + detname + " |p_{out}|" + GeV,
                      m_nPhysObjBins["DiTopPout"], m_PhysObjBins["DiTopPout"], m_nPhysObjBins["DiTopPout"], m_PhysObjBins["DiTopPout"] );
      // DeltaPhi:
      hold -> AddTH2D(JetType + "DeltaPhi" + migratag, JetType + "DeltaPhi" + migratag + ";" + ptclname + " #Delta#phi^{t#bar{t}}" + ";" + detname + " #Delta#phi^{t#bar{t}}",
                      m_nPhysObjBins["DiTopDeltaPhi"], m_PhysObjBins["DiTopDeltaPhi"], m_nPhysObjBins["DiTopDeltaPhi"], m_PhysObjBins["DiTopDeltaPhi"]);
      hold -> AddTH2D(JetType + "Delta" + migratag, JetType + "Delta" + migratag + ";" + ptclname + " #delta^{t#bar{t}}" + ";" + detname + " #delta^{t#bar{t}}" ,
                      m_nPhysObjBins["Delta"], m_PhysObjBins["Delta"], m_nPhysObjBins["Delta"], m_PhysObjBins["Delta"]);
      // Yboost:
      hold -> AddTH2D(JetType + "Yboost" + migratag, JetType + "Yboost" + migratag + ";" + ptclname + " y_{boost}^{t#bar{t}}" + ";" + detname + " y_{boost}^{t#bar{t}}",
                      m_nPhysObjBins["Yboost"], m_PhysObjBins["Yboost"], m_nPhysObjBins["Yboost"], m_PhysObjBins["Yboost"]);
      // Rttbar
      hold -> AddTH2D(JetType + "Rttbar" + migratag, JetType + "Rttbar" + migratag + ";" + ptclname + " R^{t2,t1}" + ";" + detname + " y_{boost}^{t#bar{t}}",
                      m_nPhysObjBins["Rttbar"], m_PhysObjBins["Rttbar"], m_nPhysObjBins["Rttbar"], m_PhysObjBins["Rttbar"]);
      // Chittbar:
      hold -> AddTH2D(JetType + "Chittbar" + migratag, JetType + "Chittbar" + migratag + ";" + ptclname + " #chi^{t#bar{t}}" + ";" + detname + " #chi^{t#bar{t}}",
                      m_nPhysObjBins["Chittbar"], m_PhysObjBins["Chittbar"], m_nPhysObjBins["Chittbar"], m_PhysObjBins["Chittbar"]);
      // Cos theta*
      hold -> AddTH2D(JetType + "CosThetaStar" + migratag, JetType + "CosThetaStar" + migratag + ";" + ptclname + " |cos#theta*|" + ";" + detname + " |cos#theta*|",
                        m_nPhysObjBins["CosThetaStar"], m_PhysObjBins["CosThetaStar"], m_nPhysObjBins["CosThetaStar"], m_PhysObjBins["CosThetaStar"]); // abs?!
	}

    TString topotag = "_" + hold -> GetName();
    bool IsTop = (!JetType.Contains("DiTop")) && (JetType.Contains("Top") );
    bool IsDiTop = JetType.Contains("DiTop");
    if (IsDiTop) {
        hold -> AddTH2D(JetType + "Pt" + migratag, JetType + "Pt" + migratag + ";" + ptclname + " p_{T}" + GeV + ";" + detname + " p_{T}" + GeV,
                        m_nPhysObjBins["DiTopPt"], m_PhysObjBins["DiTopPt"], m_nPhysObjBins["DiTopPt"], m_PhysObjBins["DiTopPt"]);

        hold -> AddTH2D(JetType + "Mass" + migratag, JetType + "Mass" + migratag + ";" + ptclname + " m" + GeV + ";" + detname + " m" + GeV,
                        m_nPhysObjBins["DiTopMass" + topotag], m_PhysObjBins["DiTopMass" + topotag], m_nPhysObjBins["DiTopMass" + topotag], m_PhysObjBins["DiTopMass" + topotag]);
    } else if (IsTop) {
        hold -> AddTH2D(JetType + "Pt" + migratag, JetType + "Pt" + migratag + ";" + ptclname + " p_{T}" + GeV + ";" + detname + " p_{T}" + GeV,
                        m_nPhysObjBins["TopPt" + topotag], m_PhysObjBins["TopPt" + topotag], m_nPhysObjBins["TopPt" + topotag], m_PhysObjBins["TopPt" + topotag]);
        hold -> AddTH2D(JetType + "Mass" + migratag, JetType + "Mass" + migratag + ";" + ptclname + " m" + GeV + ";" + detname + " m" + GeV,
                        nMassBins, Massmin, Massmax, nMassBins, Massmin, Massmax);
    } else {
        hold -> AddTH2D(JetType + "Pt" + migratag, JetType + "Pt" + migratag + ";" + ptclname + " p_{T}" + GeV + ";" + detname + " p_{T}" + GeV,
                        nPtBins, pTmin, pTmax, nPtBins, pTmin, pTmax);
        hold -> AddTH2D(JetType + "Mass" + migratag, JetType + "Mass" + migratag + ";" + ptclname + " m" + GeV + ";" + detname + " m" + GeV,
                        nMassBins, Massmin, Massmax, nMassBins, Massmin, Massmax);
    }

    if (IsTop || IsDiTop) {
        hold -> AddTH2D(JetType + "Rapidity" + migratag, JetType + "Rapidity" + migratag + ";" + ptclname + " y" + ";" + detname + " y",
                        m_nPhysObjBins["Rapidity"], m_PhysObjBins["Rapidity"], m_nPhysObjBins["Rapidity"], m_PhysObjBins["Rapidity"]);
    }
    else {

    hold -> AddTH2D(JetType + "Rapidity" + migratag, JetType + "Rapidity" + migratag + ";" + ptclname + " y" + ";" + detname + " y",
                    nRapidityBins, Rapiditymin, Rapiditymax, nRapidityBins, Rapiditymin, Rapiditymax);
}

    if (AddMore) {
        for (int ij = 1; ij <= AddMore; ++ij) {
            TString sij = Form("%i", ij);
            hold -> AddTH2D(JetType +  + sij + "Pt" + migratag, JetType +  + sij + "Pt" + migratag + ";" + ptclname +  + sij + " p_{T}" + GeV + ";" + detname +  + sij + " p_{T}" + GeV, nPtBins, pTmin, pTmax, nPtBins, pTmin, pTmax);
            hold -> AddTH2D(JetType +  + sij + "Rapidity" + migratag, JetType +  + sij + "Rapidity" + migratag + ";" + ptclname +  + sij + " y" + ";" + detname +  + sij + " y", nRapidityBins, Rapiditymin, Rapiditymax, nRapidityBins, Rapiditymin, Rapiditymax);
            hold -> AddTH2D(JetType +  + sij + "Mass" + migratag, JetType +  + sij + "Mass" + migratag + ";" + ptclname +  + sij + " m" + GeV + ";" + detname +  + sij + " m" + GeV, nMassBins, Massmin, Massmax, nMassBins, Massmin, Massmax);
        }
    }
    if (JetType.Contains("pseudotoplepton")) {
        // leptonic top mass vs hadronic top mass at the same level, book both m_levels;)
        TString objname = "pseudotopleptonVspseudotophadron";
        hold -> AddTH2D(objname + "Mass_" + m_level1, "pseudotopleptonVspseudotophadronMass_" + m_level1 + ";" + m_level1 + " pseudotophadron m" + GeV + ";" + m_level1 + " pseudotoplepton m" + GeV,
                        2*nMassBins, Massmin, Massmax, 2*nMassBins, Massmin, Massmax);
        hold -> AddTH2D(objname + "Mass_" + m_level2, "pseudotopleptonVspseudotophadronMass_" + m_level2 + ";" + m_level2 + " pseudotophadron m" + GeV + ";" + m_level2 + " pseudotoplepton m" + GeV,
                        2*nMassBins, Massmin, Massmax, 2*nMassBins, Massmin, Massmax);
        hold -> AddTH2D(objname + "Pt_" + m_level1, "pseudotopleptonVspseudotophadronPt_" + m_level1 + ";" + m_level1 + " pseudotophadron p_{T}" + GeV + ";" + m_level1 + " pseudotoplepton p_{T}" + GeV,
                        2*nPtBins, pTmin, pTmax, 2*nPtBins, pTmin, pTmax);
        hold -> AddTH2D(objname + "Pt_" + m_level2, "pseudotopleptonVspseudotophadronPt_" + m_level2 + ";" + m_level2 + " pseudotophadron p_{T}" + GeV + ";" + m_level2 + " pseudotoplepton p_{T}" + GeV,
                        2*nPtBins, pTmin, pTmax, 2*nPtBins, pTmin, pTmax);
        hold -> AddTH2D(objname + "Rapidity_" + m_level1, "pseudotopleptonVspseudotophadronRapidity_" + m_level1 + ";" + m_level1 + " pseudotophadron y" + GeV + ";" + m_level1 + " pseudotoplepton y" + GeV,
                        2*nRapidityBins, Rapiditymin, Rapiditymax, 2*nRapidityBins, Rapiditymin, Rapiditymax);
        hold -> AddTH2D(objname + "Rapidity_" + m_level2, "pseudotopleptonVspseudotophadronRapidity_" + m_level2 + ";" + m_level2 + " pseudotophadron y" + GeV + ";" + m_level2 + " pseudotoplepton y" + GeV,
                        2*nRapidityBins, Rapiditymin, Rapiditymax, 2*nRapidityBins, Rapiditymin, Rapiditymax);
    }

	if (FillSub) {
	  int nTau32Bins = 20;
	  double Tau32min = 0;
	  double Tau32max = 1;
	  int nTau21Bins = nTau32Bins;
	  double Tau21min = Tau32min;
      double Tau21max = Tau32max;
      hold -> AddTH2D(JetType + "Tau32" + migratag, JetType + "Tau32" + migratag + ";" + ptclname + " #tau_{32}" + ";" + detname + " #tau_{32}",
                      nTau32Bins, Tau32min, Tau32max, nTau32Bins, Tau32min, Tau32max);
      hold -> AddTH2D(JetType + "Tau21" + migratag, JetType + "Tau21" + migratag + ";" + ptclname + " #tau_{21}" + ";" + detname + " #tau_{21}",
                      nTau21Bins, Tau21min, Tau21max, nTau21Bins, Tau21min, Tau21max);
      if (AddMore) {
          for (int ij = 1; ij <= AddMore; ++ij) {
              TString sij = Form("%i", ij);
              hold -> AddTH2D(JetType + sij + "Tau32" + migratag, JetType + sij + "Tau32" + migratag + ";" + ptclname +  sij + " #tau_{32}" + ";" + detname + sij + " #tau_{32}", nTau32Bins, Tau32min, Tau32max, nTau32Bins, Tau32min, Tau32max);
              hold -> AddTH2D(JetType +  sij + "Tau21" + migratag, JetType + sij + "Tau21" + migratag + ";" + ptclname +  sij + " #tau_{21}" + ";" + detname + sij + " #tau_{21}", nTau21Bins, Tau21min, Tau21max, nTau21Bins, Tau21min, Tau21max);
        }
      }
	} // substructure
	 
     //gDirectory->cd("../");
	 
   } else
     cerr << "ERROR HistoMaker::MakeDijetHistos: making dijet histos!" << endl;

  return;
}


// __________________________________________________

void HistoMaker::MakeAlljetMigrations(TString cutname, TString m_level1, TString m_level2)
{

    auto iter = m_holders.find(cutname);
    if (iter != m_holders.end()) {

        auto hold = iter->second;
        if (!gDirectory->GetDirectory(hold->GetName()))
            gDirectory->mkdir(hold->GetName());
        gDirectory->cd(hold->GetName());
        TString migradir = "migrations";
        if (!gDirectory->GetDirectory(migradir)) {
            gDirectory->mkdir(migradir);
        }
        gDirectory->cd(migradir);

        // Nov2020
        this -> MakeGlobalMigrations(cutname, "Global", m_level1, m_level2);

        //this -> MakeJetMigrations(cutname, "Jet", m_level1, m_level2, false, 4);
        //this -> MakeJetMigrations(cutname, "LJet", m_level1, m_level2, false, 4);
        this -> MakeJetMigrations(cutname, "Top", m_level1, m_level2, true, 4);
        this -> MakeJetMigrations(cutname, "W", m_level1, m_level2, true, 4);
        this -> MakeJetMigrations(cutname, "DiTop", m_level1, m_level2, false, 2);
        this -> MakeJetMigrations(cutname, "FourTop", m_level1, m_level2, false, 2);

        gDirectory->cd("../");
        gDirectory->cd("../");

    } else
        cerr << "ERROR HistoMaker::MakeMigrations: making migration histos!" << endl;

    return;
}


// __________________________________________________

void HistoMaker::FillJetHistos(TString cutname, TString JetType, vector<TMyLorentzVector*> jets, double weight, kJetTypes jettype, bool FillSubs)
{
	

  int nMaxJetsToFill = ComputeMaxJetsToFill(JetType);

  HistoHolder* hold = 0;
  try {
    hold = m_holders[cutname];
  }  catch (exception& e) {
    cout << "ERROR in HistoMaker::FillJetHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
   }
  
    TString jname = m_level + JetType;
	
	int j=0;
	int nitems = 0;

    // fill leading and subleading jet mass and dR, 2.7.2021
    double mjj = -1;
    double dR = -1.;
    if (jets.size() > 1) {
        mjj = (*jets[0] + *jets[1]).M();
        dR = jets[0]->DeltaR(*jets[1]);
	if (hold -> GetTH1D("LeadDi" + jname + "Mass"))
	    hold -> FillTH1D("LeadDi" + jname + "Mass", mjj, weight);
	if (hold -> GetTH1D("LeadDi" + jname + "dR"))
	    hold -> FillTH1D("LeadDi" + jname + "dR", dR, weight);
        for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
            TString reptag = Form("_rep%i", irep);
            if (hold -> GetTH1D("LeadDi" + jname + "Mass" + reptag))
                hold -> FillTH1D("LeadDi" + jname + "Mass" + reptag, mjj, weight*m_replicaWeights[irep]);
            if (hold -> GetTH1D("LeadDi" + jname + "dR" + reptag))
                hold -> FillTH1D("LeadDi" + jname + "dR" + reptag, dR, weight*m_replicaWeights[irep]);

        }
    }

	for (auto jet : jets) {
	  
        if (!PassedJetType(jettype, jet))
            continue;


	  j=j+1;
	  nitems++;

	  hold -> FillTH1D(jname + "Pt", jet->Pt(), weight);
	  hold -> FillTH1D(jname + "Rapidity", jet->Rapidity(), weight);
	  hold -> FillTH1D(jname + "Mass", jet->M(), weight);

      // fill automatically replicas if exist
      for (unsigned int irep = 0; irep < m_nReplicas; ++irep) {
          TString reptag = Form("_rep%i", irep);
          if (hold -> GetTH1D(jname + "Pt" + reptag)) {
              hold -> FillTH1D(jname + "Pt" + reptag, jet->Pt(), weight*m_replicaWeights[irep]);
          }
          if (hold -> GetTH1D(jname + "Rapidity" + reptag)) {
              hold -> FillTH1D(jname + "Rapidity" + reptag, jet->Rapidity(), weight*m_replicaWeights[irep]);
          }
          if (hold -> GetTH1D(jname + "Mass" + reptag)) {
              hold -> FillTH1D(jname + "Mass" + reptag, jet->M(), weight*m_replicaWeights[irep]);
          }

          if (hold -> GetTH1D(jname + "Pt" + cMidPoint + reptag)) {
              hold -> FillTH1D(jname + "Pt" + cMidPoint + reptag, jet->Pt(), weight*m_replicaWeights[irep]);
          }
          if (hold -> GetTH1D(jname + "Rapidity" + cMidPoint + reptag)) {
              hold -> FillTH1D(jname + "Rapidity" + cMidPoint + reptag, jet->Rapidity(), weight*m_replicaWeights[irep]);
          }
          if (hold -> GetTH1D(jname + "Mass" + cMidPoint + reptag)) {
              hold -> FillTH1D(jname + "Mass" + cMidPoint + reptag, jet->M(), weight*m_replicaWeights[irep]);
          }
      }

	  
	  if (FillSubs) {
          hold -> FillTH1D(jname + "Tau21", jet->Tau21(), weight);
          hold -> FillTH1D(jname + "Tau32", jet->Tau32(), weight);
          hold -> FillTH2D(jname + "Tau21VsMass", jet->M(), jet->Tau21(), weight);
          hold -> FillTH2D(jname + "Tau32VsMass", jet->M(), jet->Tau32(), weight);
          hold -> FillTH2D(jname + "Tau21VsTau32", jet->Tau32(), jet->Tau21(), weight);
          // 2.7.2021
          hold -> FillTH2D(jname + "Tau21VsPt", jet->Pt(), jet->Tau21(), weight);
          hold -> FillTH2D(jname + "Tau32VsPt", jet->Pt(), jet->Tau32(), weight);
          hold -> FillTH2D(jname + "PtVsMass", jet->M(), jet->Pt(), weight);
          hold -> FillTH2D(jname + "PtVsMass", jet->M(), jet->Pt(), weight);
	  }
  
      if (j <= nMaxJetsToFill) {
	    hold -> FillTH1D(jname + TString(Form("%i", j)) + "Pt", jet->Pt(), weight);
	    hold -> FillTH1D(jname + TString(Form("%i", j)) + "Mass", jet->M(), weight);
	    hold -> FillTH1D(jname + TString(Form("%i", j)) + "Rapidity", jet->Rapidity(), weight);
	    
	    if (FillSubs) {
		  hold -> FillTH1D(jname + TString(Form("%i", j)) + "Tau21", jet->Tau21(), weight);	
		  hold -> FillTH1D(jname + TString(Form("%i", j)) + "Tau32", jet->Tau32(), weight);
	    }
	    
	  } 
	    
    } // jet loop
    
    // cout << "Will try to fill " << (jname + "N").Data() << endl;
	hold -> FillTH1D(jname + "N", nitems, weight);
	

    return;
}



// __________________________________________________

void HistoMaker::FillSingleObjectHistos(TString cutname, TString JetType, TMyLorentzVector& obj, double weight)
{
	
  HistoHolder* hold = 0;
  try {
    hold = m_holders[cutname];
  }  catch (exception& e) {
    cout << "ERROR in HistoMaker::FillJetHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
  }
	
  TString jname = m_level + JetType;
  hold -> FillTH1D(jname + "Pt", obj.Pt(), weight);
  hold -> FillTH1D(jname + "Rapidity", obj.Rapidity(), weight);
  hold -> FillTH1D(jname + "Mass", obj.M(), weight);

  return;
}


// __________________________________________________

void HistoMaker::FillJetMigrations(TString cutname, TString jettag, vector<TMyLorentzVector*> ptcljets,
                                   vector<TMyLorentzVector*> detjets, double weight,
                                   int AddMore, kJetTypes jettype, bool FillSubs)
{

	HistoHolder* hold = 0;
	try {
		hold = m_holders[cutname];
	}  catch (exception& e) {
    cout << "ERROR in HistoMaker::FillJetMigrations: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return; 
   }

    TString migratag = "_" + m_level1 + "_" + m_level2;
  
    int j=0;
	for (auto ptcljet : ptcljets) {
	  j=j+1;
      if (!PassedJetType(jettype, ptcljet))
          continue;
      if (j <= detjets.size()) {
          if (!PassedJetType(jettype, detjets[j-1]))
              continue;

          hold -> FillTH2D(jettag + "Pt" + migratag, ptcljet->Pt(), detjets[j-1]->Pt(), weight);
          hold -> FillTH2D(jettag + "Rapidity" + migratag, ptcljet->Rapidity(), detjets[j-1]->Rapidity(), weight);
          hold -> FillTH2D(jettag + "Mass" + migratag, ptcljet->M(), detjets[j-1]->M(), weight);
          if (FillSubs) {
              hold -> FillTH2D(jettag + "Tau21" + migratag, ptcljet->Tau21(), detjets[j-1]->Tau21(), weight);
              hold -> FillTH2D(jettag + "Tau32" + migratag, ptcljet->Tau32(), detjets[j-1]->Tau32(), weight);
          }
          if (AddMore && j <= AddMore) {
              hold -> FillTH2D(jettag + TString(Form("%i", j)) + "Pt" + migratag, ptcljet->Pt(), detjets[j-1]->Pt(), weight);
              hold -> FillTH2D(jettag + TString(Form("%i", j)) + "Rapidity" + migratag, ptcljet->Rapidity(), detjets[j-1]->Rapidity(), weight);
              hold -> FillTH2D(jettag + TString(Form("%i", j)) + "Mass" + migratag, ptcljet->M(), detjets[j-1]->M(), weight);
              if (FillSubs) {
                  hold -> FillTH2D(jettag + TString(Form("%i", j)) + "Tau21" + migratag, ptcljet->Tau21(), detjets[j-1]->Tau21(), weight);
                  hold -> FillTH2D(jettag + TString(Form("%i", j)) + "Tau32" + migratag, ptcljet->Tau32(), detjets[j-1]->Tau32(), weight);
              }
          }
	  
	  } // sizes consistent
	
    } // jet loop

  return;
}

// __________________________________________________

void HistoMaker::FillAlljetMigrations(TString cutname,
                                      vector<TMyLorentzVector*> ptclLjets,
                                      vector<TMyLorentzVector*> ptcltops,
                                      vector<TMyLorentzVector*> ptclDiTops,
                                      vector<TMyLorentzVector*> ptclFourTops,
                                      vector<TMyLorentzVector*> detLjets,
                                      vector<TMyLorentzVector*> dettops,
                                      vector<TMyLorentzVector*> detDiTops,
                                      vector<TMyLorentzVector*> detFourTops, double weight)
{

    // fill eventually leading and subleading tops?!
    // would have to keep track of the possible swap in pTs between BS...

    HistoHolder* hold = 0;
    try {
        hold = m_holders[cutname];
    }  catch (exception& e) {
    cout << "ERROR in HistoMaker::FillJetMigrations: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
    return;
   }

    this -> FillJetMigrations(cutname, "W", ptclLjets, detLjets, weight, 2, kWTagJets, true);
    this -> FillJetMigrations(cutname, "Top", ptcltops, dettops, weight, 4, kallJets, true);
    this -> FillJetMigrations(cutname, "DiTop", ptclDiTops, detDiTops, weight, 2, kallJets, false);
    this -> FillJetMigrations(cutname, "FourTop", ptclFourTops, detFourTops, weight, 0, kallJets, false);


  return;

}


// __________________________________________________
// __________________________________________________

void HistoMaker::FillJESHistos(TString cutname, vector<TMyLorentzVector*> ptcljets, vector<TMyLorentzVector*> detjets, bool IsClosure, bool isLJets, double weight)
{
    HistoHolder* hold = 0;
    try {
        hold = m_holders[cutname];
    }  catch (exception& e) {
        cout << "ERROR in HistoMaker::FillJESHisto: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
        return;
    }

    TString jtag = isLJets ? "LJets" : "Jets";

    for (auto djet : detjets) {
        int iptcl = -1;
        double minDR = 9999;
        int iptclMin = -1;
        for (auto pjet : ptcljets) {
            iptcl++;
            double dr = djet -> DeltaR(*pjet);
            if (dr < 0.2 && dr < minDR) {
                iptclMin = iptcl;
                minDR = dr;
            }
        } // particle jets
        if (iptclMin >= 0) {
            double pt = ptcljets[iptclMin] -> Pt();
            double E = ptcljets[iptclMin] -> E();
            double eta = ptcljets[iptclMin] -> Eta();
            if (pt > 0. && E > 0.) {
                double Rpt = djet->Pt() / pt;
                double RE = djet->E() / E;
                double DEta = djet->Eta() - eta;
                double DPhi = fabs(djet->DeltaPhi(*ptcljets[iptclMin]));

                TString closuretag = "";
                if (IsClosure)
                    closuretag = "Closure";
                hold -> FillTH2D("JESRPt" + closuretag + jtag + "Eta", djet -> Eta(), Rpt, weight);
                hold -> FillTH2D("JESRPt" + closuretag + jtag + "Pt", djet -> Pt(), Rpt, weight);
                hold -> FillTH2D("JESRPt" + closuretag + jtag + "E", djet -> E(), Rpt, weight);
                hold -> FillTH3D("JESRPt" + closuretag + jtag + "EtaPt", djet -> Eta(), djet -> Pt(), Rpt, weight);
                hold -> FillTH3D("JESRPt" + closuretag + jtag + "EtaE", djet -> Eta(), djet -> E(), Rpt, weight);

                hold -> FillTH2D("JESRE" + closuretag + jtag + "Eta", djet -> Eta(), RE, weight);
                hold -> FillTH2D("JESRE" + closuretag + jtag + "Pt", djet -> Pt(), RE, weight);
                hold -> FillTH2D("JESRE" + closuretag + jtag + "E", djet -> E(), RE, weight);
                hold -> FillTH3D("JESRE" + closuretag + jtag + "EtaPt", djet -> Eta(), djet -> Pt(), RE, weight);
                hold -> FillTH3D("JESRE" + closuretag + jtag + "EtaE", djet -> Eta(), djet -> E(), RE, weight);

                /*
                hold -> FillTH2D("JESDEta" + closuretag + jtag + "Eta", djet -> Eta(), DEta, weight);
                hold -> FillTH2D("JESDEta" + closuretag + jtag + "Pt", djet -> Pt(), DEta, weight);
                hold -> FillTH2D("JESDEta" + closuretag + jtag + "E", djet -> E(), DEta, weight);

                hold -> FillTH2D("JESDPhi" + closuretag + jtag + "Eta", djet -> Eta(), DPhi, weight);
                hold -> FillTH2D("JESDPhi" + closuretag + jtag + "Pt", djet -> Pt(), DPhi, weight);
                hold -> FillTH2D("JESDPhi" + closuretag + jtag + "E", djet -> E(), DPhi, weight);
                */

            } // nonzero pT,E

        } // matched

    } // detjets
}



// __________________________________________________


void HistoMaker::FillMigrations(TString cutname, vector<TMyLorentzVector*> ptclLjets, vector<TMyLorentzVector*> ptcljets, vector<TMyLorentzVector*> ptclDijets, 
									             vector<TMyLorentzVector*> detLjets, vector<TMyLorentzVector*> detjets, vector<TMyLorentzVector*> detDijets, double weight)
{
	
  this -> FillJetMigrations(cutname, "Jet", ptcljets, detjets, weight, 4);
  this -> FillJetMigrations(cutname, "LJet", ptclLjets, detLjets, weight, 4);
  this -> FillJetMigrations(cutname, "Dijet", ptclDijets, detDijets, weight, 2); // TO CHECK!
											 
  return;											 
}



// __________________________________________________

void HistoMaker::FillSingleTagDRHistos( mindrdetadphi drt, mindrdetadphi drw, double pt, double weight)
{
    TString cutname = "TagHistos";

    HistoHolder* hold = 0;
    try {
        hold = m_holders[cutname];
    }  catch (exception& e) {
        cout << "ERROR in HistoMaker::FillSingleTagDRHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
        return;
    }

    if (drt.mindr < 100.) {
        hold -> FillTH1D("MinDRLjetTopParton" + m_level, drt.mindr, weight);
        hold -> FillTH1D("MinDRLjetTopPartonZoom" + m_level, drt.mindr, weight);
        hold -> FillTH2D("LjetPtVsMinDRLjetTopParton" + m_level, drt.mindr, pt, weight);
        hold -> FillTH2D("LjetPtVsMinDRLjetTopPartonZoom" + m_level, drt.mindr, pt, weight);

        hold -> FillTH1D("MinDEtaLjetTopParton" + m_level, drt.mindeta, weight);
        hold -> FillTH1D("MinDEtaLjetTopPartonZoom" + m_level, drt.mindeta, weight);
        hold -> FillTH2D("LjetPtVsMinDEtaLjetTopParton" + m_level, drt.mindeta, pt, weight);
        hold -> FillTH2D("LjetPtVsMinDEtaLjetTopPartonZoom" + m_level, drt.mindeta, pt, weight);

        hold -> FillTH1D("MinDPhiLjetTopParton" + m_level, drt.mindphi, weight);
        hold -> FillTH1D("MinDPhiLjetTopPartonZoom" + m_level, drt.mindphi, weight);
        hold -> FillTH2D("LjetPtVsMinDPhiLjetTopParton" + m_level, drt.mindphi, pt, weight);
        hold -> FillTH2D("LjetPtVsMinDPhiLjetTopPartonZoom" + m_level, drt.mindphi, pt, weight);
    }
    if (drw.mindr < 100.) {
        hold -> FillTH1D("MinDRLjetWParton" + m_level, drw.mindr, weight);
        hold -> FillTH1D("MinDRLjetWPartonZoom" + m_level, drw.mindr, weight);
        hold -> FillTH2D("LjetPtVsMinDRLjetWParton" + m_level, drw.mindr, pt, weight);
        hold -> FillTH2D("LjetPtVsMinDRLjetWPartonZoom" + m_level, drw.mindr, pt, weight);

        hold -> FillTH1D("MinDEtaLjetWParton" + m_level, drw.mindeta, weight);
        hold -> FillTH1D("MinDEtaLjetWPartonZoom" + m_level, drw.mindeta, weight);
        hold -> FillTH2D("LjetPtVsMinDEtaLjetWParton" + m_level, drw.mindeta, pt, weight);
        hold -> FillTH2D("LjetPtVsMinDEtaLjetWPartonZoom" + m_level, drw.mindeta, pt, weight);

        hold -> FillTH1D("MinDPhiLjetWParton" + m_level, drw.mindphi, weight);
        hold -> FillTH1D("MinDPhiLjetWPartonZoom" + m_level, drw.mindphi, weight);
        hold -> FillTH2D("LjetPtVsMinDPhiLjetWParton" + m_level, drw.mindphi, pt, weight);
        hold -> FillTH2D("LjetPtVsMinDPhiLjetWPartonZoom" + m_level, drw.mindphi, pt, weight);
    }
}

// __________________________________________________

void HistoMaker::FillSingleTagKinemHistos(TString objname, TMyLorentzVector *jet, double weight)
{

    TString cutname = "TagHistos";

    HistoHolder* hold = 0;
    try {
        hold = m_holders[cutname];
    }  catch (exception& e) {
        cout << "ERROR in HistoMaker::FillSingleTagKinemHistos: did not find HistoHolder " << cutname.Data() << " Error: " << e.what() << endl;
        return;
    }
    hold -> FillTH1D(objname + "Pt" + m_level, jet -> Pt(), weight);
    hold -> FillTH1D(objname + "Mass" + m_level, jet -> M(), weight);
    hold -> FillTH1D(objname + "Tau32" + m_level, jet -> Tau32(), weight);
    hold -> FillTH1D(objname + "Tau21" + m_level, jet -> Tau21(), weight);
    hold -> FillTH2D(objname + "Tau32VsMass" + m_level, jet -> M(), jet -> Tau32(), weight);
    hold -> FillTH2D(objname + "Tau21VsMass" + m_level, jet -> M(), jet -> Tau21(), weight);
    hold -> FillTH2D(objname + "Tau32VsPt" + m_level, jet -> Pt(), jet -> Tau32(), weight);
    hold -> FillTH2D(objname + "Tau21VsPt" + m_level, jet -> Pt(), jet -> Tau21(), weight);
    hold -> FillTH2D(objname + "PtVsMass" + m_level, jet -> M(), jet -> Pt(), weight);

}

// __________________________________________________

void HistoMaker::MakeSingleDRHistos(TString cutname)
{
    auto iter = m_holders.find(cutname);
    int npt = 100;
    double minpt = 0.;
    double maxpt = 500;
    if (iter != m_holders.end()) {
        auto hold = iter->second;
        if (!gDirectory->GetDirectory(hold->GetName())) {
            gDirectory->mkdir(hold->GetName());
        }
        gDirectory->cd(hold->GetName());

        TString vars[] = {"MinDR", "MinDEta", "MinDPhi"};
        TString labs[] = {"#DeltaR_{min}", "#Delta#eta_{min}", "#Delta#phi_{min}"};
        double ranges[][2] = {  {0, 5}, {-4, 4}, {0, TMath::Pi()} };
        int n = sizeof(vars) / sizeof(TString);
        for (int i = 0; i < n; ++i) {

            hold -> AddTH1D(vars[i] + "LjetTopParton" + m_level, vars[i] + "LjetTopParton;" + labs[i] + "^{t,Ljet}, " + m_level + " level", 250, ranges[i][0], ranges[i][1]);
            hold -> AddTH1D(vars[i] + "LjetWParton" + m_level, vars[i] + "LjetWParton;" + labs[i] + "^{W,Ljet}" + m_level + " level", 250, ranges[i][0], ranges[i][1]);
            hold -> AddTH1D(vars[i] + "LjetTopPartonZoom" + m_level, vars[i] + "LjetTopParton;" + labs[i] + "^{t,Ljet}, " + m_level + " level", 100, ranges[i][0]/5, ranges[i][1]/5);
            hold -> AddTH1D(vars[i] + "LjetWPartonZoom" + m_level, vars[i] + "LjetWParton;" + labs[i] + "^{W,Ljet}, " + m_level + " level", 100, ranges[i][0]/5, ranges[i][1]/5);

            hold -> AddTH2D("LjetPtVs" + vars[i] + "LjetTopParton" + m_level, "LjetPt" + vars[i] + "LjetTopParton;" + labs[i] + "^{t,Ljet};" +  m_level + " p_{T}^{Ljet}", 250, ranges[i][0], ranges[i][1], npt, minpt, maxpt );
            hold -> AddTH2D("LjetPtVs" + vars[i] + "LjetWParton" + m_level, "LjetPtVs" + vars[i] + "LjetWParton;" + labs[i] + "^{W,Ljet};" +  m_level + " p_{T}^{Ljet}", 250, ranges[i][0], ranges[i][1], npt, minpt, maxpt);
            hold -> AddTH2D("LjetPtVs" + vars[i] + "LjetTopPartonZoom" + m_level, "LjetPtVs" + vars[i] + "LjetTopPartonZoom;" + labs[i] + "^{t,Ljet};" +  m_level + " p_{T}^{Ljet}", 100, ranges[i][0]/5, ranges[i][1]/5, npt, minpt, maxpt);
            hold -> AddTH2D("LjetPtVs" + vars[i] + "LjetWPartonZoom" + m_level, "LjetPtVs" + vars[i] + "LjetWPartonZoom;" + labs[i] + "^{W,Ljet};" +  m_level + " p_{T}^{Ljet}", 100, ranges[i][0]/5, ranges[i][1]/5, npt, minpt, maxpt);
        }
        gDirectory->cd("../");
    }
}

// __________________________________________________


void HistoMaker::MakeSingleTagKinemHistos(TString cutname, TString objname)
{

    int nMassBins = 79;
    double Massmin = 0;
    double Massmax = 300;

    TString objtag = objname;
    objtag.ReplaceAll("ttruthMatched", "t-match");
    objtag.ReplaceAll("nontMatched", "t-non-match");
    objtag.ReplaceAll("WtruthMatched", "W-match");
    objtag.ReplaceAll("nonWMatched", "W-match");
    objtag.ReplaceAll("notTopTag", "non-top-tag ");
    objtag.ReplaceAll("TopTag", "top-tag ");
    objtag.ReplaceAll("notWTag", "non-W-tag ");
    objtag.ReplaceAll("WTag", "W-tag ");
    objtag.ReplaceAll("Ljets", " ");
    auto iter = m_holders.find(cutname);
    if (iter != m_holders.end()) {
        auto hold = iter->second;
        if (!gDirectory->GetDirectory(hold->GetName())) {
            gDirectory->mkdir(hold->GetName());
        }
        gDirectory->cd(hold->GetName());
        //cout << "Making tagging histos in directory " << gDirectory->GetName() << endl;
        hold -> AddTH1D(objname + "Pt" + m_level, objname + "Pt" + m_level + ";" + objtag + " p_{T}^{J}" + m_level, m_nPhysObjBins["LJetPt"], (double*)m_PhysObjBins["LJetPt"]);
        hold -> AddTH1D(objname + "Mass" + m_level, objname + "Mass"  + m_level + ";" + objtag + " m^{J}" + m_level, m_nPhysObjBins["LJetMass"], (double*)m_PhysObjBins["LJetMass"]);
        hold -> AddTH1D(objname + "Tau32" + m_level, objname + "Tau32" + m_level + ";" + objtag + " #tau_{32}^{J}", 100, m_Taumin, m_Taumax);
        hold -> AddTH1D(objname + "Tau21" + m_level, objname + "Tau21" +  + m_level + ";" + objtag + " #tau_{21}^{J}", 100, m_Taumin, m_Taumax);

        hold -> AddTH2D(objname + "Tau32VsPt" + m_level, objname + "Tau32VsMass;" + objtag + m_level + " p_{T}^{J};" + objtag + m_level + " #tau_{32}^{J}",
                        m_nPhysObjBins["LJetPt"], (double*)m_PhysObjBins["LJetPt"], m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(objname + "Tau21VsPt" + m_level, objname + "Tau21VsMass;" + objtag +  + m_level + " p_{T}^{J};" + objtag + m_level + " #tau_{21}^{J}",
                        m_nPhysObjBins["LJetPt"], (double*)m_PhysObjBins["LJetPt"], m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(objname + "Tau32VsMass" + m_level, objname + "Tau32VsMass;" + objtag + m_level + " m^{J};" + objtag + m_level + " #tau_{32}^{J}",
                        nMassBins, Massmin, Massmax, m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(objname + "Tau21VsMass" + m_level, objname + "Tau21VsMass;" + objtag +  + m_level + " m^{J};" + objtag + m_level + " #tau_{21}^{J}",
                        nMassBins, Massmin, Massmax, m_nTauBins, m_Taumin, m_Taumax);
        hold -> AddTH2D(objname + "PtVsMass" + m_level, objname + "PtVsMass;" + objtag + m_level + " m^{J};" + objtag + m_level + " p_{T}^{J}",
                        nMassBins, Massmin, Massmax, m_nPhysObjBins["LJetPt"], (double*)m_PhysObjBins["LJetPt"]);
/*
        hold -> AddTH2D(objname + "Tau21VsTau32" + m_level, objname + "Tau21VsTau32;" + objtag +  + m_level + " #tau_{32};" + objtag + m_level + " #tau_{21}^{J}",
                        m_nTauBins, m_Taumin, m_Taumax, m_nTauBins, m_Taumin, m_Taumax);
*/

        // add tau32 vs tau21?
        // add th3d tauXY, pT, M?

        gDirectory->cd("../");
    }
}
// __________________________________________________
void HistoMaker::MakeTaggingHistos()
{

    TString objs[] = {"ttruthMatchedLjets", "WtruthMatchedLjets",
                      "TopTagttruthMatchedLjets", "WTagWtruthMatchedLjets",
                      "TopTagnontMatchedLjets", "WTagnonWMatchedLjets",
                      "notTopTagttruthMatchedLjets", "notWTagWtruthMatchedLjets",
                      "notTopTagnontMatchedLjets", "notWTagnonWMatchedLjets",
                      "TopTagLjets", "WTagLjets"};
    // fragile!
    TString tagdirname = "TagHistos";
    this -> AddCutLevel(tagdirname);
    for (auto obj : objs) {
        this -> MakeSingleTagKinemHistos(tagdirname, obj);
    }
    this -> MakeSingleDRHistos(tagdirname);
}

// __________________________________________________
// __________________________________________________

						

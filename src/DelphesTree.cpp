#include "DelphesTree.h"


DelphesTree::DelphesTree(TString path, TString pattern, TString weightStr) : fChain(0)
{

    TChain *tree = new TChain("Delphes");
    cout << "Adding pattern " << (path + "/*" + pattern + "*.root").Data() << endl;
    tree -> Add(path + "/*" + pattern + "*.root");
    Init(tree);

    m_ljets = kUndef;
    m_PseudotopType = kNonePseudo;
    m_isSignal = pattern.Contains("pp_2tj_allhad_NLO");
    int ppindex = pattern.Index("pp_");
    if (ppindex < 0)
        ppindex = pattern.Index("zp_");
    if (ppindex < 0)
        ppindex = 0;
    int Eindex = pattern.Index("14TeV");
    m_sampleTag = pattern(ppindex,Eindex - ppindex - 1);
    cout << "Decifered sample tag " << m_sampleTag.Data() << endl;

    m_ApplyXsectWeights = false;
    weightStr.ToUpper();
    if (weightStr.Contains("WEIGHT")) {
        cout << "Asked to apply the xsection weights!" << endl;
        m_ApplyXsectWeights = true;
    }

}

DelphesTree::~DelphesTree()
{
    if (!fChain) return;
    // delete fChain->GetCurrentFile();
}

void DelphesTree::Init(TChain *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    m_rand = new TRandom3();

    m_JES_calibFile = new TFile("data/JESfits.root", "read");
    m_JES_LJets_response = (TF2*) m_JES_calibFile -> Get("fit2d_JESRPtLJetsEtaE");
    m_JES_SJets_response = (TF2*) m_JES_calibFile -> Get("fit2d_JESRPtJetsEtaE");

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    m_treeReader = new ExRootTreeReader(fChain);

    m_branchEvent = m_treeReader->UseBranch("Event");
    //m_branchParticle = m_treeReader->UseBranch("Particle");
    m_branchElectron = m_treeReader->UseBranch("Electron");
    m_branchPhoton = m_treeReader->UseBranch("Photon");
    m_branchMuon = m_treeReader->UseBranch("Muon");
    m_branchJetJES = m_treeReader->UseBranch("JetJES");
    m_branchLJet = m_treeReader->UseBranch("LJet");
    m_branchGenJet = m_treeReader->UseBranch("GenJet");
    m_branchGenLJet = m_treeReader->UseBranch("GenLJet");

    m_branchMet = m_treeReader->UseBranch("MissingET");
    m_branchGenMet = m_treeReader->UseBranch("GenMissingET");

    m_branchGenElectron = m_treeReader->UseBranch("GenElectron");
    m_branchGenPhoton = m_treeReader->UseBranch("GenPhoton");
    m_branchGenMuon = m_treeReader->UseBranch("GenMuon");

    m_branchGenBhadrons = m_treeReader->UseBranch("GenBhadrons");
    m_branchGenTop = m_treeReader->UseBranch("GenTop");
    m_branchGenW = m_treeReader->UseBranch("GenW");


}


#include "DelphesTree.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[])

{
	
   TString path = "";
  TString pattern = "";
  if (argc < 3) {
      cout << "Usage " << argv[0] << " path pattern LjetType [weight] [channel] [pseudotopType]" << endl;
      cout << "     LjetTypeStr can be Default, SoftDropped, Pruned, Trimmed" << endl;
      //cout << "      where channel is ejets or mujets" << endl;
      //cout << "      where pseudotop type is Standard, CloseMt, SameMt, TwoStep, TwoStepII, BestBsAndNu" << endl;
	  return -1;
  }
  path = argv[1];
  pattern = argv[2];
  TString LjetTypeStr = argv[3];
  TString weightStr = "";
  if (argc > 4)
      weightStr = argv[4];


  TString channel = "";
  if (argc > 5)
      channel = argv[5];
  TString pseudo = "";
  if (argc > 6)
      pseudo = argv[6];

  // !!!
  // 2.9.2020
  int nReplicas = 0; // 0 // 100; // 100

  cout << "Running over " << path.Data() << "/*" << pattern.Data() << "*.root" << endl;

  //   local:
  // gSystem->Load("/home/qitek/install/delphes/libDelphesNoFastJet.so");
  // IGA:
  gSystem->Load("/home/qitek/work/delphes/libDelphes.so");


  DelphesTree *delphesTree = new DelphesTree(path, pattern, weightStr);
  TString outtag = pattern;
  outtag.ReplaceAll("?", "X");
  outtag.ReplaceAll("*", "STAR");
  // Loop! ;-)
  delphesTree -> LoopBoosted(outtag, LjetTypeStr, nReplicas);
  delete delphesTree;


  /*
  TChain *chain = new TChain("Delphes");
  chain->Add(path + "*" + pattern + "*.root");
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  TClonesArray *branchLJet = treeReader->UseBranch("LJet");
  Jet *LJet;

  Long64_t allEntries = treeReader->GetEntries();
  // Loop over all events
  for(Long64_t entry = 0; entry < allEntries; ++entry)
  {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      // Loop over all electrons in event
      for(int i = 0; i < branchLJet->GetEntriesFast(); ++i)
      {
          LJet = (Jet*) branchLJet->At(i);
          cout << " LJetMass: "           << LJet->Mass << " "
               << " LJet_SoftDropped: "   << LJet->SoftDroppedJet.M() << " "
               << " LJet_SoftDroppedP4: " << LJet->SoftDroppedP4[0].M()
               << " LJet_PrunedP4: "      << LJet->PrunedP4[0].M()
               << " LJet_TrimmedP4: "     << LJet->TrimmedP4[0].M()
               << endl;

      } // jets

  } // entries
*/


  return 0;
}

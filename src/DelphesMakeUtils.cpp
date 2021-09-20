#include "DelphesTree.h"
// __________________________________________________________________


bool DelphesTree::IsBtaggedPtclLevel(double jpt, double jeta, double jphi, double jm)
{
    double DR = -1;
    TLorentzVector jet;
    jet.SetPtEtaPhiM(jpt, jeta, jphi, jm);
    TLorentzVector bhad;
    for (int ib = 0; ib < m_branchGenBhadrons->GetEntries(); ++ib)  {
        GenParticle *GenBhadrons = (GenParticle*) m_branchGenBhadrons -> At(ib);
        if (GenBhadrons->PT < 5.) // c.f. ATLAS paper!
            continue;
        bhad.SetPtEtaPhiM(GenBhadrons->PT, GenBhadrons->Eta, GenBhadrons->Phi, GenBhadrons->Mass);
        DR = jet.DeltaR(bhad);
        if (DR < 0.4)
            return true;
    }
    return false;
}

// __________________________________________________________________

double JEScorr(double eta, double pt)
{
    // Delphes jet response formula???
    return sqrt( pow(3.0 - 0.2*(abs(eta)), 2) / pt + 1.0 );
}



// __________________________________________________________________

bool IsTopTagged(double mass, double tau21, double tau32)
{
    //return ( mass - kmt) < kTopMassWindow) && tau32 < kTau32Cut);
    return ( mass < kmtcut2 &&
             mass > kmtcut1 &&
             tau32 < kTau32tCut2 &&
             tau32 > kTau32tCut1 &&
             tau21 < kTau21tCut2 &&
             tau21 > kTau21tCut1
             );
}

// __________________________________________________________________

bool IsWTagged(double mass, double tau21, double tau32)
{
//    return ( (fabs(mass - kmW) < kWMassWindow) && tau21 < kTau21Cut);
    return ( mass < kmWcut2 &&
             mass > kmWcut1 &&
             tau32 > kTau32WCut1 &&
             tau32 < kTau32WCut2 &&
             tau21 < kTau21WCut2 &&
             tau21 > kTau21WCut1
             );
}

// __________________________________________________________________

vector<TMyLorentzVector*> DelphesTree::MakeGenTops(double ptcut, double etacut, int status)
{

    // make a vector of nlast truth tops in the truth record chain
    vector<TMyLorentzVector*> gentops;
    /*
   cout << "=== in DelphesTree::MakeGenTops ===" << endl;
   cout << " Parton tops: " << GenTop->size << endl;
   for (unsigned int i = 0; i < GenTop->size; ++i) {
     cout << "GenTop " << i << " status " << GenTop->Status[i] << " pdgid: " << GenTop->PID[i] << " pT: " << GenTop->PT[i] << endl;
    // before FSR:
     //     if (GenTop->Status[i] < 29) {
     //       TMyLorentzVector *gentop = new TMyLorentzVector(GenTop->PT[i], GenTop->Eta[i], GenTop->Phi[i], GenTop->Mass[i], 0., 0.);
     //       gentops.push_back(gentop);
     //}
   }
   */
    // take last tops in chain:
   //// cout << "GenTops: " << m_branchGenTop->GetEntries() << endl;
    //int i = GenTop->size - 1;
    //while (i < m_branchGenTop->GetEntries() && i >= m_branchGenTop->GetEntries() - nlast && i >= 0) {
    for (unsigned int i = 0; i < m_branchGenTop->GetEntries(); ++i) {
        GenParticle* GenTop = (GenParticle*) m_branchGenTop->At(i);
        ////cout << "  i=" << i << " charge=" << GenTop->Charge << " PID=" << GenTop->PID << " status=" << GenTop->Status << " pt=" << GenTop->PT << " eta=" << fabs(GenTop->Eta) << " phi="<< fabs(GenTop->Phi) <<  endl;
            if (GenTop->PT > ptcut && fabs(GenTop->Eta) < etacut && GenTop->Status == status )  {
                //cout << "  i=" << i << " charge=" << GenTop->Charge << " PID=" << GenTop->PID << " status=" << GenTop->Status << " pt=" << GenTop->PT << " eta=" << fabs(GenTop->Eta) << " phi="<< fabs(GenTop->Phi) <<  endl;
                TMyLorentzVector *gentop = new TMyLorentzVector(GenTop->PT, GenTop->Eta, GenTop->Phi, GenTop->Mass, 0., 0.);
                gentops.push_back(gentop);
            }
            //i--;
    }
    return gentops;
}


// __________________________________________________________________

vector<TMyLorentzVector*> DelphesTree::MakeGenWs(double ptcut, double etacut, int status)
{
    vector<TMyLorentzVector*> genWs;
    // take last Ws in chain:
    //// cout << "GenWs: " << m_branchGenW->GetEntries() << endl;
    //int i = GenW->size - 1;
    //while (i < m_branchGenW->GetEntries() && i >= m_branchGenW->GetEntries() - nlast && i >= 0) {
    for (unsigned int i = 0; i < m_branchGenW->GetEntries(); ++i) {
        GenParticle* GenW = (GenParticle*) m_branchGenW->At(i);
       //// cout << "  i=" << i << " charge=" << GenW->Charge << " PID=" << GenW->PID << " status=" << GenW->Status << " pt=" << GenW->PT << " eta=" << fabs(GenW->Eta) << " phi="<< fabs(GenW->Phi) <<  endl;
        if (GenW->PT > ptcut && fabs(GenW->Eta) < etacut  && GenW->Status == status)  {
            //cout << "  i=" << i << " charge=" << GenW->Charge << " PID=" << GenW->PID << " status=" << GenW->Status << " pt=" << GenW->PT << " eta=" << fabs(GenW->Eta) << " phi="<< fabs(GenW->Phi) <<  endl;
            TMyLorentzVector *genW = new TMyLorentzVector(GenW->PT, GenW->Eta, GenW->Phi, GenW->Mass, 0., 0.);
            genWs.push_back(genW);
        }
        //i--;
    }
    return genWs;
}

// __________________________________________________________________

vector<TMyLorentzVector*> DelphesTree::MakeGenDiTops(vector<TMyLorentzVector*> gentops)
{
    // warning, so far only the first pair is returned!!
    vector<TMyLorentzVector*> ditops;
    if (gentops.size() > 1) {
        TMyLorentzVector *mygendijet = new TMyLorentzVector(*gentops[0]);
        (*mygendijet)  += *gentops[1];
        ditops.push_back(mygendijet);
    }
    return ditops;
}


// __________________________________________________________________

vector<TMyLorentzVector*> DelphesTree::MakeJets(bool DetLevel, bool LJets,
                                                double ptcut, double etacut, bool CorrectJES,
                                                bool runFakeBtag,
                                                kLJetChoice LJetChoice)
{

    // NOTE: IN CASE I GO BACK TO LJETS CHANNEL EVER AGAIN,
    // THE JETS MUST BE CLEANED OFF FOR JETS OVERLAPPING WITH LEPTONS!
    // see the vector<TMyLorentzVector*>  DelphesTree::RemoveParticleJetsOverlapWithLeadLepton(kLjets ljets, vector<TMyLorentzVector*>& jets, double ptcut, bool debug)
    // function in alljets_analysis, which is a ljets pseudotop NIM analysis, BTW;-)


    int nn = 0;
    if (LJets)
        nn = DetLevel ? m_branchLJet->GetEntries() : m_branchGenLJet->GetEntries();
    else
        nn = DetLevel ? m_branchJetJES->GetEntries() : m_branchGenJet->GetEntries();

    vector<TMyLorentzVector*> topjets;
    vector<TMyLorentzVector*> wjets;
    vector<TMyLorentzVector*> ljets;
    vector<TMyLorentzVector*> smalljets;

    // loop over jets

    for (int ijet = 0; ijet < nn; ++ijet) {

        double jesc = 1.;
        if (CorrectJES && DetLevel) {
            if (LJets) {
                Jet *LJet = (Jet*)m_branchLJet -> At(ijet);
                jesc = 1./ m_JES_LJets_response -> Eval(LJet->PT * cosh(LJet->Eta), LJet->Eta);
            } else {
                // hmmm, was just Jet here in the old MakeClass fwk?! jk 14.10.2020
                // ...this means the JES correction was not evaluated at proper energy!
                Jet *jet = (Jet*)m_branchJetJES -> At(ijet);
                jesc = 1./ m_JES_SJets_response -> Eval(jet->PT * cosh(jet->Eta), jet->Eta);
            }
        }

        /*
        if (LJets && nn > 0) {
          cout << "TEST: ";
          if (DetLevel)
        cout << "detector";
          else
        cout << "particle";
          cout << " Jet_Constituents " << Jet_Constituents << " and size: " << Jet_Constituents->GetSize() << endl;
          }
        */
            double Pt, Eta, Phi, M, tau32, tau21;
            bool btag = false;
            bool wtag = false;
            bool toptag = false;
            if (LJets) {
                // large jets
                Jet *LJet = (Jet*)m_branchLJet -> At(ijet);
                Jet *GenLJet = (Jet*)m_branchGenLJet -> At(ijet);

                if (LJetChoice == kDefault || LJetChoice == kResolved) {
                    Pt  = (DetLevel ? jesc*LJet->PT   : GenLJet->PT);
                    Eta = (DetLevel ? LJet->Eta       : GenLJet->Eta);
                    Phi = (DetLevel ? LJet->Phi       : GenLJet->Phi);
                    M   = (DetLevel ? jesc*LJet->Mass : GenLJet->Mass);
                } else if (LJetChoice == kSoftDroppedP4) {
                    Pt  = (DetLevel ? jesc*LJet->SoftDroppedP4[0].Pt()   : GenLJet->SoftDroppedP4[0].Pt());
                    Eta = (DetLevel ? LJet->SoftDroppedP4[0].Eta()       : GenLJet->SoftDroppedP4[0].Eta());
                    Phi = (DetLevel ? LJet->SoftDroppedP4[0].Phi()       : GenLJet->SoftDroppedP4[0].Phi());
                    M   = (DetLevel ? jesc*LJet->SoftDroppedP4[0].M()    : GenLJet->SoftDroppedP4[0].M());
                } else if (LJetChoice == kPrunedP4) {
                    Pt  = (DetLevel ? jesc*LJet->PrunedP4[0].Pt()   : GenLJet->PrunedP4[0].Pt());
                    Eta = (DetLevel ? LJet->PrunedP4[0].Eta()       : GenLJet->PrunedP4[0].Eta());
                    Phi = (DetLevel ? LJet->PrunedP4[0].Phi()       : GenLJet->PrunedP4[0].Phi());
                    M   = (DetLevel ? jesc*LJet->PrunedP4[0].M()    : GenLJet->PrunedP4[0].M());
                } else if (LJetChoice == kTrimmedP4) {
                    Pt  = (DetLevel ? jesc*LJet->TrimmedP4[0].Pt()   : GenLJet->TrimmedP4[0].Pt());
                    Eta = (DetLevel ? LJet->TrimmedP4[0].Eta()       : GenLJet->TrimmedP4[0].Eta());
                    Phi = (DetLevel ? LJet->TrimmedP4[0].Phi()       : GenLJet->TrimmedP4[0].Phi());
                    M   = (DetLevel ? jesc*LJet->TrimmedP4[0].M()    : GenLJet->TrimmedP4[0].M());
                }
                tau32 = (DetLevel ? LJet->Tau[2] / LJet->Tau[1] : GenLJet->Tau[2] / GenLJet->Tau[1]);
                tau21 = (DetLevel ? LJet->Tau[1] / LJet->Tau[0] : GenLJet->Tau[1] / GenLJet->Tau[0]);
                // we don't need to btag Ljets; JK 30.3.2020, coronavirus times
                if (DetLevel) {
                    btag = LJet->BTag;// ZERO?? use matching to a small btagged jet?
                } else {
                    btag = IsBtaggedPtclLevel(GenLJet->PT, GenLJet->Eta, GenLJet->Phi, GenLJet->Mass);
                }
                wtag = IsWTagged(M, tau21, tau32);
                toptag = IsTopTagged(M, tau21, tau32);

            } else {
                // small Jets
                Jet *JetJES = (Jet*)m_branchJetJES -> At(ijet);
                Jet *GenJet = (Jet*)m_branchGenJet -> At(ijet);
                Pt  = (DetLevel ? jesc*JetJES->PT   : GenJet->PT);
                Eta = (DetLevel ? JetJES->Eta       : GenJet->Eta);
                Phi = (DetLevel ? JetJES->Phi       : GenJet->Phi);
                M   = (DetLevel ? jesc*JetJES->Mass : GenJet->Mass);
                tau32 = (DetLevel ? JetJES->Tau[2] / JetJES->Tau[1] : GenJet->Tau[2] / GenJet->Tau[1]);
                tau21 = (DetLevel ? JetJES->Tau[1] / JetJES->Tau[0] : GenJet->Tau[1] / GenJet->Tau[0]);
                if (DetLevel) {
                    btag = JetJES->BTag;
                } else {
                    btag = IsBtaggedPtclLevel(GenJet->PT, GenJet->Eta, GenJet->Phi, GenJet->Mass);
                }
                if (!btag && runFakeBtag) {
                    // add some b-tagged jets, emulate fake btag on light; both ptcl and detector levels
                    //
                    double randy = m_rand->Uniform(0., 1.);
                    if (randy < kbtagFakeEff)
                        btag = true;
                }
            }

            //cout << "  Pt=" << Pt << " Eta=" << Eta << " Phi=" << Phi << " M=" << M << " tau32=" << tau32 << " tau21=" << tau21 << endl;

            if ( (Pt > ptcut) && (fabs(Eta) < etacut) ) {

                if (toptag) {
                    TMyLorentzVector *jet = new TMyLorentzVector(Pt, Eta, Phi, M, tau32, tau21);
                    topjets.push_back(jet);
                }

                if (wtag) {
                    TMyLorentzVector *jet = new TMyLorentzVector(Pt, Eta, Phi, M, tau32, tau21);
                    wjets.push_back(jet);
                }

                if (!LJets) {
                    TMyLorentzVector *jet = new TMyLorentzVector(Pt, Eta, Phi, M, 0., 0., btag);
                    /*
                    if (DetLevel)
                        cout << " ........small det jet " << ijet << " :: Pt=" << Pt << " Eta=" << Eta << " Phi=" << Phi << " M=" << M << " btag=" <<  btag  << endl;
                    else
                        cout << " ........small ptcl jet " << ijet << " :: Pt=" << Pt << " Eta=" << Eta << " Phi=" << Phi << " M=" << M << " btag=" <<  btag  << endl;
                        */
                    smalljets.push_back(jet);
                } else {
                    /*
                    if (DetLevel)
                        cout << " ........large det jet " << ijet << " :: Pt=" << Pt << " Eta=" << Eta << " Phi=" << Phi << " M=" << M << " btag=" <<  btag  << endl;
                    else
                        cout << " ........large ptcl jet " << ijet << " :: Pt=" << Pt << " Eta=" << Eta << " Phi=" << Phi << " M=" << M << " btag=" <<  btag  << endl;
                        */
                    TMyLorentzVector *jet = new TMyLorentzVector(Pt, Eta, Phi, M, tau32, tau21, btag, wtag, toptag);
                    ljets.push_back(jet);
                }
                //cout << "  jet Pt,Eta,Phi,M: " << jet->Pt() << " " <<  jet->Eta() << " " <<  jet->Phi() << " " <<  jet->M() << " " <<  endl;
            }

            //cout << "  Created " << jets.size() << " jets." << endl;

        }



    // return for the moment all large jets!
    if (LJets) {
        //if (DetLevel) printJetArray(ljets, "unsorted LJets");
        if (CorrectJES)
            std::sort(ljets.begin(), ljets.end(), cmpPt);
        //if (DetLevel) printJetArray(ljets, "sorted LJets");
        return ljets;
    } else {
        //if (DetLevel) printJetArray(smalljets, "unsorted jets");
        if (CorrectJES)
            std::sort(smalljets.begin(), smalljets.end(), cmpPt);
        //if (DetLevel) printJetArray(smalljets, "sorted jets");
        return smalljets;
    }

   
}




// __________________________________________________________________

vector<TMyLorentzVector*> DelphesTree::MakeTops(vector<TMyLorentzVector*> &jets, vector<TMyLorentzVector*> &ljets, double mindrbjet, int debug)
{

    int ntops = Counttoptags(ljets);
    int nWs = CountWtags(ljets);
    int nbs = CountBtags(jets);

    if (debug > 0)
    cout << "In MakeTop: "
         << " ntops=" << ntops
         << " nWs=" << nWs
         << " nbs=" << nbs
         << endl;


    vector<TMyLorentzVector*> tops;
    TMyLorentzVector* top;

    for (auto ljet : ljets) {
        if (ljet -> Toptag()) {
            top = new TMyLorentzVector(*ljet);
            top->SetToptag(true);
            tops.push_back(top);
        }
    } // ljets

    if (ntops < 2 && nbs > 0 && nWs > 0) {
        if (debug > 0) cout << "OK, trying semiboosted reconstruction!" << endl;

        for (auto ljet : ljets) {
            if (ljet -> Wtag()) {
                if (debug > 0) cout << "  found Wtagged jet, will loop over " << jets.size() << " small jets..."  << endl;
                int bijet = -1;
                int ijet = -1;
                double mindr = 999;
                for (auto jet : jets) {
                    ijet++;
                    double dr = ljet -> DeltaR(*jet);
                    if (debug > 0) cout << "   btag=" << jet -> Btag() << " dr=" << dr   << endl;
                    if (jet -> Btag() && !jet->GetUsed() && dr > mindrbjet) {
                        // now make sure this does not overlap with another top/W tagged ljet!
                        // later also check overlap with lepton?? ;-)
                        bool overlaps = false;
                        for (auto ljet2 : ljets) {
                            if ( (ljet2 -> Toptag() || ljet2 -> Wtag()) && ljet2->DeltaR(*jet) < mindrbjet) {
                                overlaps = true;
                                break;
                            }
                        }
                        if (!overlaps && dr < mindr) {
                            mindr = dr;
                            bijet = ijet;
                        } // overlap and better dr
                        else {
                            if (debug > 0) cout << "  Sorry, overlap with another top tagged jet!" << endl;
                        }
                    } // btag
                } // jets

                // add the best closet-in-DR btagged small jet tot he Wtag jet to form the top candidate;)
                if (bijet >= 0) {
                    if (debug > 0) cout << "    adding a bjet to the Wtag jet!" << endl;
                    jets[bijet] -> SetUsed(true);
                    top = new TMyLorentzVector(*ljet);
                    top -> Add(*jets[bijet]);
                    top -> SetToptag(true);
                    top -> SetWtag(false);
                    tops.push_back(top);
                }

            } // Wtag
        } // ljets
    }

    return tops;
}


// __________________________________________________________________

vector<TMyLorentzVector*> DelphesTree::MakeDiTops(vector<TMyLorentzVector*> &tops, int debug)
{
    vector<TMyLorentzVector*> ditops;
    TMyLorentzVector* ditop;
    if (tops.size() > 1) {
        bool lookingforfirst = true;
        for (auto top : tops) {
            //if (top -> Toptag()) {
                if (lookingforfirst) {
                  ditop = new TMyLorentzVector(*top);
                } else {
                    ditop -> Add(*top);
                    ditop -> SetToptag(true);
                    ditops.push_back(ditop);
                }
                lookingforfirst = !lookingforfirst;
            //}


        }
    }
    if (debug > 0) cout << "MakeDiTops: got " << tops.size() << " tops, returning " << ditops.size() << "ditops!" << endl;
    return ditops;
}

// __________________________________________________________________

vector<TMyLorentzVector*> DelphesTree::MakeFourTops(vector<TMyLorentzVector*> tops, int debug)
{
    vector<TMyLorentzVector*> fourtops;
    TMyLorentzVector* fourtop = 0;
    if (tops.size() > 3) {
        int nobjs = 0;
        for (auto top : tops) {
            if (top -> Toptag()) {
                if (nobjs == 0) {
                  fourtop = new TMyLorentzVector(*top);
                } else {
                    fourtop -> Add(*top);
                }
                nobjs++;
                if (nobjs == 4) {
                    fourtops.push_back(fourtop);
                    nobjs = 0;
                }
            } // tagger
        } // tops
    } // size
    return fourtops;
}

// __________________________________________________________________

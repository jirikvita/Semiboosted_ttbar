#include "KinemUtils.h"

// __________________________________________________________________

bool cmpPt(const TMyLorentzVector *a, const TMyLorentzVector *b)
{
    return a->Pt() > b->Pt();
}

// __________________________________________________________________

void printJetArray(vector<TMyLorentzVector*> jets, TString tag)
{
    int i = 0;
    cout << "Printing jet array " << tag.Data() << endl;
    for (auto jet : jets) {
        cout << "Pt=" << jet->Pt() << endl;
        i++;
    }
}


// __________________________________________________________________

bool MatchedWithinDR(TMyLorentzVector *vec, vector<TMyLorentzVector*> lvecs, mindrdetadphi &vals, double DRcut)
{
   vals.mindr = 999.;
    double dr = -1.;
    bool found = false;
    for (auto lvec : lvecs) {
        dr = vec->DeltaR(*lvec);
        if (dr < vals.mindr) {
            vals.mindr = dr;
            vals.mindeta = vec->Eta() - lvec->Eta();
            vals.mindphi = fabs(vec->DeltaPhi(*lvec));
        } if (dr < DRcut) {
            found = true;
        }
    }
    return found;
}

// __________________________________________________________________
// __________________________________________________________________
// __________________________________________________________________


//
//
// code below adopted from the RIVET routine
// 
// 


// __________________________________________________________________

double computeneutrinoz(const TMyLorentzVector lepton, double metx, double mety, bool invertSolutions, double& discriminant) {
  //computing z component of neutrino momentum given lepton and met
  double pzneutrino;
  double met = sqrt ( sqr(metx) + sqr(mety) );
  double m_W = 80.399; // in GeV, given in the paper
  double k = (( sqr( m_W ) - sqr( lepton.M() ) ) / 2 ) + (lepton.Px() * metx + lepton.Py() * mety);
  double a = sqr ( lepton.E() )- sqr ( lepton.Pz() );
  double b = -2*k*lepton.Pz();
  double c = sqr( lepton.E() ) * sqr( met ) - sqr( k );
  discriminant = sqr(b) - 4 * a * c;
  double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
  if (discriminant < 0)  pzneutrino = - b / (2 * a); //if the discriminant is negative
  else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
    double absquad[2];
    for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
    if (absquad[0] < absquad[1]) {
      if (!invertSolutions)
	pzneutrino = quad[0];
      else 
	pzneutrino = quad[1];
    } else {
      if (!invertSolutions)
	pzneutrino = quad[1];
      else 
	pzneutrino = quad[0];
    }
  } // discr
           
  return pzneutrino;
}

// __________________________________________________________________

double computeneutrinoz_from_mt(double mthad, const TMyLorentzVector lepton, const TMyLorentzVector blep, double metx, double mety, bool invertSolutions, double& discriminant) {
  //computing z component of neutrino momentum from mthad = mtlep constraint
  double pzneutrino;
  
  TVector3 vlepton = lepton.Vect(); 
  TVector3 vblep = blep.Vect(); 

  double DeltaSq = vlepton.Dot(vblep);
  double SigmaSq = sqr(blep.M()) + sqr(lepton.M()) - 2 * DeltaSq + 2 * lepton.E()*blep.E() - 2 * ( (lepton.Px() + blep.Px() ) * metx + ( lepton.Py() + blep.Py() ) * mety );
  double EsumSq = sqr(lepton.E() + blep.E());
  double pzSum = lepton.Pz() + blep.Pz();
  double massSum = sqr(mthad) - SigmaSq;
  double METSq = sqr(metx) + sqr(mety);
      
  double a = 4 * ( sqr(pzSum) - EsumSq );
  double b = 4 * massSum*pzSum;
  double c = sqr(massSum) - 4*METSq*EsumSq;
      
  discriminant = sqr(b) - 4 * a * c;
  double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
  if (discriminant < 0)  pzneutrino = - b / (2 * a); //if the discriminant is negative
  else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
    double absquad[2];
    for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
    if (absquad[0] < absquad[1]) {
      if (!invertSolutions)
	pzneutrino = quad[0];
      else 
	pzneutrino = quad[1];
    } else {
      if (!invertSolutions)
	pzneutrino = quad[1];
      else 
	pzneutrino = quad[0];
    }
  } // discr
           
  return pzneutrino;
}


// __________________________________________________________________

double mT(TMyLorentzVector l, TMyLorentzVector nu) {
  return sqrt( 2 * l.Pt() * nu.Pt() * (1 - cos(l.DeltaPhi(nu) ) ) );
}

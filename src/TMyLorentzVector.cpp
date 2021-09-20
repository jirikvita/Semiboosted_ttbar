#include "TMyLorentzVector.h"

// __________________________________________________________________


int CountBtags(vector<TMyLorentzVector*> &jets)
{
    int n = 0;
    for (auto jet : jets) {
        if (!jet) continue;
        if (jet->Btag())
            n++;
    }
    return n;
}

int CountWtags(vector<TMyLorentzVector*> &jets)
{
    int n = 0;
    for (auto jet : jets) {
        if (!jet) continue;
        if (jet->Wtag())
            n++;
    }
    return n;
}

int Counttoptags(vector<TMyLorentzVector*> &jets)
{
    int n = 0;
    for (auto jet : jets) {
        if (!jet) continue;
        if (jet->Toptag())
            n++;
    }
    return n;
}

// __________________________________________________________________

TMyLorentzVector::TMyLorentzVector() : m_btag(false), m_wtag(false), m_toptag(false), m_used(false)
{
    this -> Set(0., 0., 0., 0., 0., 0.);

}
// __________________________________________________________________

TMyLorentzVector::TMyLorentzVector(TLorentzVector vec) : TLorentzVector(vec)
{
  m_tau32 = 0;  
  m_tau21 = 0;
  m_btag = false;
  m_wtag = false;
  m_toptag = false;
  m_used = false;

}
// __________________________________________________________________

TMyLorentzVector::TMyLorentzVector(TLorentzVector vec, double tau32, double tau21, bool btag, bool wtag, bool toptag): TLorentzVector(vec)
{
  m_tau32 = tau32;  
  m_tau21 = tau32;
  m_btag = btag;
  m_wtag = wtag;
  m_toptag = toptag;
  m_used = false;

}

// __________________________________________________________________

void TMyLorentzVector::Set(TLorentzVector vec, double tau32, double tau21, bool btag, bool wtag, bool toptag)
{
  this -> SetPtEtaPhiM(vec.Pt(), vec.Eta(), vec.Phi(), vec.M());
  m_tau32 = tau32;  
  m_tau21 = tau32;
  m_btag = btag;
  m_wtag = wtag;
  m_toptag = toptag;
  m_used = false;
}

// __________________________________________________________________

void TMyLorentzVector::Set(TLorentzVector vec, bool btag, bool wtag, bool toptag)
{
  this -> SetPtEtaPhiM(vec.Pt(), vec.Eta(), vec.Phi(), vec.M());
  m_tau32 = 0.;  
  m_tau21 = 0.;
  m_btag = btag;
  m_wtag = wtag;
  m_toptag = toptag;
  m_used = false;
}
// __________________________________________________________________

TMyLorentzVector::TMyLorentzVector(double Pt, double Eta, double Phi, double M, double tau32, double tau21, bool btag, bool wtag, bool toptag)
{
   this -> Set(Pt, Eta, Phi, M, tau32, tau21);
   m_btag = btag;
   m_wtag = wtag;
   m_toptag = toptag;
   m_used = false;
}
  
 // __________________________________________________________________
  
  void TMyLorentzVector::Set(double Pt, double Eta, double Phi, double M, double tau32, double tau21, bool btag, bool wtag, bool toptag)
  {
    this -> SetPtEtaPhiM(Pt, Eta, Phi, M);
    m_tau32 = tau32;  
    m_tau21 = tau21;
    m_btag = btag;
    m_wtag = wtag;
    m_toptag = toptag;
    m_used = false;
  }
// __________________________________________________________________

void TMyLorentzVector::Set(double Pt, double Eta, double Phi, double M, bool btag, bool wtag, bool toptag)
  {
    this -> SetPtEtaPhiM(Pt, Eta, Phi, M);
    m_tau32 = 0.;  
    m_tau21 = 0.;
    m_btag = btag;
    m_wtag = wtag;
    m_toptag = toptag;
    m_used = false;
}

// __________________________________________________________________
  
TMyLorentzVector::~TMyLorentzVector() {}

// __________________________________________________________________
  bool TMyLorentzVector::Btag() {return m_btag;}
  bool TMyLorentzVector::Wtag() {return m_wtag;}
  bool TMyLorentzVector::Toptag() {return m_toptag;}
  double TMyLorentzVector::Tau32() {return m_tau32;}
  double TMyLorentzVector::Tau21() {return m_tau21;}
  bool TMyLorentzVector::GetUsed() {return m_used;}
  // __________________________________________________________________

  void TMyLorentzVector::SetBtag(bool btag) {m_btag = btag;}
  void TMyLorentzVector::SetWtag(bool wtag) {m_wtag = wtag;}
  void TMyLorentzVector::SetToptag(bool toptag) {m_toptag = toptag;}
  void TMyLorentzVector::SetTau32(double tau32) {m_tau32 = tau32;}
  void TMyLorentzVector::SetTau21(double tau21) {m_tau21 = tau21;}
  void TMyLorentzVector::SetUsed(bool used) { m_used = used;}
// __________________________________________________________________

TLorentzVector TMyLorentzVector::GetTLorentzVector()
{
     TLorentzVector vec;
     vec.SetPtEtaPhiM(this->Pt(), this->Eta(), this->Phi(), this->M());  
     return vec;
}

// __________________________________________________________________

void TMyLorentzVector::Add(TLorentzVector vec)
{
  
  TLorentzVector tmp = this->GetTLorentzVector();
  tmp += vec;
  this->Set(tmp);
  return;
  
}


// __________________________________________________________________

TMyLorentzVector TMyLorentzVector::operator+ (TMyLorentzVector right)
{

  right.Add(this->GetTLorentzVector());
  return right;
}



// __________________________________________________________________
// __________________________________________________________________
// __________________________________________________________________



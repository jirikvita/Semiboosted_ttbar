#include "TLorentzVector.h"

#ifndef TMyLorentzVector_HPP
#define TMyLorentzVector_HPP

#include <vector>

using std::vector;

class TMyLorentzVector : public TLorentzVector
{
  private:
  
    double m_tau32;
    double m_tau21;
    bool m_btag;
    bool m_toptag;
    bool m_wtag;
    bool m_used; // e.g. for building composed objecs
    
  public:
  
  TMyLorentzVector();
  TMyLorentzVector(TLorentzVector vec);
  TMyLorentzVector(TLorentzVector vec, double tau32, double tau21, bool btag = false, bool wtag = false, bool toptag = false);
  TMyLorentzVector(double Pt, double Eta, double Phi, double M, double tau32, double tau21, bool btag = false, bool wtag = false, bool toptag = false);

  void Set(double Pt, double Eta, double Phi, double M, double tau32, double tau21, bool btag = false, bool wtag = false, bool toptag = false);
  void Set(double Pt, double Eta, double Phi, double M, bool btag = false, bool wtag = false, bool toptag = false);
  void Set(TLorentzVector vec, double tau32, double tau21, bool btag = false, bool wtag = false, bool toptag = false);
  void Set(TLorentzVector vec, bool btag = false, bool wtag = false, bool toptag = false);

  ~TMyLorentzVector();

  bool Btag();
  bool Wtag();
  bool Toptag();
  double Tau32();
  double Tau21();

  void SetBtag(bool btag);
  void SetToptag(bool toptag);
  void SetWtag(bool wtag);
  void SetTau32(double tau32);
  void SetTau21(double tau21);

  void SetUsed(bool used);
  bool GetUsed();


  inline TLorentzVector GetTLorentzVector();
  inline void Add(TLorentzVector vec);
  TMyLorentzVector operator+ (TMyLorentzVector right);
  
};



int CountBtags(vector<TMyLorentzVector*> &jets);
int CountWtags(vector<TMyLorentzVector*> &jets);
int Counttoptags(vector<TMyLorentzVector*> &jets);



#endif

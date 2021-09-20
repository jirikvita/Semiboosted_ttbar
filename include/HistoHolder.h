/* File HistoHolder.hpp
 *
 * Created       : Mon Jan 16 09:01:36 CST 2006
 * Author        : Dag Gillberg, Jiri Kvita
 * Purpose       : 
 * Last modified : 
 * Comments      : 
 */

#ifndef HistoHolder_HPP_
#define HistoHolder_HPP_

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TString.h"

#include <map>
#include <iostream>


using namespace std;

const TString cMidPoint = "_denser";

class HistoHolder // : public TObject 
{
  
 private:
  TString m_name;
  map<TString, TH1D*>     _HistoMapTH1D;
  map<TString, TProfile*> _HistoMapTProfile;
  map<TString, TH2D*>     _HistoMapTH2D;
  map<TString, TH3D*>     _HistoMapTH3D;
  
  bool m_jetsOnly;
    
  
 public:
  
  // Constructors, destructor: 
  HistoHolder(TString name, bool jetsOnly = false);
  TString GetName();
  
  inline TH1D* GetTH1D( TString HName ) 
  { 
    return _HistoMapTH1D[HName]; 
  }
  inline TH2D* GetTH2D( TString HName ) 
  { 
    return _HistoMapTH2D[HName]; 
  }
  inline TH3D* GetTH3D( TString HName )
  {
    return _HistoMapTH3D[HName];
  }
  inline TProfile* GetTProfile( TString HName ) 
  { 
    return _HistoMapTProfile[HName]; 
  }

  bool JetsOnly();
  bool SetJetsOnly(bool jetsOnly);
  
  void AddTH1D( TString HName, TString HTitle, Int_t Nbins, Double_t xlow, Double_t xhigh);
  void AddTH1D ( TString HName, TString HTitle, Int_t Nbins, Double_t *xbins, TString addMidpoints = "false");

  void AddTH2D ( TString HName, TString HTitle, 
		 Int_t Nbinsx, Double_t xlow, Double_t xhigh, 
		 Int_t Nbinsy, Double_t ylow, Double_t yhigh);
  void AddTH2D ( TString HName, TString HTitle,
         Int_t Nbinsx, Double_t *xbins,
         Int_t Nbinsy, Double_t *ybins);

  void AddTH2D ( TString HName, TString HTitle,
                 Int_t Nbinsx, Double_t xlow, Double_t xhigh,
                 Int_t Nbinsy, Double_t *ybins);
  void AddTH2D ( TString HName, TString HTitle,
                 Int_t Nbinsx, Double_t *xbins,
                 Int_t Nbinsy, Double_t ylow, Double_t yhigh);

  void AddTH3D ( TString HName, TString HTitle,
                 Int_t Nbinsx, Double_t *xbins,
                 Int_t Nbinsy, Double_t *ybins,
                 Int_t Nbinsz, Double_t* zbins
                 );
  void AddTH3D ( TString HName, TString HTitle,
                 Int_t Nbinsx, Double_t xlow, Double_t xhigh,
                 Int_t Nbinsy, Double_t ylow, Double_t yhigh,
                 Int_t Nbinsz, Double_t zlow, Double_t zhigh
                 );

  void AddTProfile ( TString HName, TString HTitle, 
		     Int_t Nbinsx, Double_t xlow, Double_t xhigh, 
		     Double_t ylow, Double_t yhigh);
  void AddTProfile ( TString HName, TString HTitle, 
		     Int_t Nbinsx, Double_t *xbins, 
		     Double_t ylow, Double_t yhigh);
  void AddTProfile ( TString HName, TString HTitle, 
		     Int_t Nbinsx, Double_t xlow, Double_t xhigh);
  void AddTProfile ( TString HName, TString HTitle, 
		     Int_t Nbinsx, Double_t *xbins);
  
  void FillTH1D( TString HName, Double_t val, Double_t weight = 1.);
  void FillTH2D( TString HName, Double_t valx, Double_t valy, Double_t weight = 1.);
  void FillTH3D( TString HName, Double_t valx, Double_t valy, Double_t valz, Double_t weight = 1.);
  void FillTProfile( TString HName, Double_t valx, Double_t valy, Double_t weight = 1.);
  
  void ThrowError(TString HName);

  // To be called from "inputFileOpened"
  //  ClassDef(HistoHolder,1);
  
};



#endif

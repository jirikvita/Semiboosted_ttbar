/* File HistoHolder.hpp
 *
 * Created       : Mon Jan 16 09:01:36 CST 2006
 * Author        : Dag Gillberg, Jiri Kvita
 * Purpose       : 
 * Last modified : 
 * Comments      : 
 */

#include <iostream>

#include "HistoHolder.h"

using namespace std;


// Constructor
HistoHolder::HistoHolder(TString name, bool jetsOnly) : m_name(name), m_jetsOnly(jetsOnly)
{
}


TString HistoHolder::GetName()
{
	return m_name;
}

// -- adding -- //

// with option:




bool HistoHolder::JetsOnly()
{
  return m_jetsOnly;
}

bool HistoHolder::SetJetsOnly(bool jetsOnly)
{
  m_jetsOnly = jetsOnly;
}




void HistoHolder::AddTH3D ( TString HName, TString HTitle,
                            Int_t Nbinsx, Double_t *xbins,
                            Int_t Nbinsy, Double_t *ybins,
                            Int_t Nbinsz, Double_t* zbins
                            )
{
    _HistoMapTH3D[HName] = new TH3D(HName, HTitle, Nbinsx, xbins, Nbinsy, ybins, Nbinsz, zbins);
    _HistoMapTH3D[HName] -> Sumw2();
}

void HistoHolder::AddTH3D ( TString HName, TString HTitle,
                            Int_t Nbinsx, Double_t xlow, Double_t xhigh,
                            Int_t Nbinsy, Double_t ylow, Double_t yhigh,
                            Int_t Nbinsz, Double_t zlow, Double_t zhigh
                            )
{
    _HistoMapTH3D[HName] = new TH3D(HName, HTitle, Nbinsx, xlow, xhigh, Nbinsy, ylow, yhigh, Nbinsz, zlow, zhigh);
    _HistoMapTH3D[HName] -> Sumw2();
}

void HistoHolder::AddTH2D ( TString HName, TString HTitle,
               Int_t Nbinsx, Double_t xlow, Double_t xhigh,
               Int_t Nbinsy, Double_t *ybins)
{
  //cout << "    booking TH2 A " << HName.Data() << " nbinx=" << Nbinsx << "  nbinsy=" << Nbinsy << endl;
  _HistoMapTH2D[HName] = new TH2D(HName, HTitle, Nbinsx, xlow, xhigh, Nbinsy, ybins);
    _HistoMapTH2D[HName] -> Sumw2();
}

void HistoHolder::AddTH2D ( TString HName, TString HTitle,
               Int_t Nbinsx, Double_t *xbins,
               Int_t Nbinsy, Double_t ylow, Double_t yhigh)
{
  //    cout << "    booking TH2 B " << HName.Data() << " nbinx=" << Nbinsx << "  nbinsy=" << Nbinsy << endl;
    _HistoMapTH2D[HName] = new TH2D(HName, HTitle, Nbinsx, xbins, Nbinsy, ylow, yhigh);
    _HistoMapTH2D[HName] -> Sumw2();
}


// standard: (follows THXD constructors)

void HistoHolder::AddTH1D( TString HName, TString HTitle, Int_t Nbins, Double_t xlow, Double_t xhigh)
{
   // cout << "    booking " << HName.Data() << " " << Nbins << " " << xlow << " " << xhigh << endl;
    _HistoMapTH1D[HName] =  new TH1D(HName, HTitle, Nbins, xlow, xhigh);
    // cout << _HistoMapTH1D[HName] << endl;
    _HistoMapTH1D[HName] -> Sumw2();
}

void HistoHolder::AddTH1D( TString HName, TString HTitle, Int_t Nbins, Double_t *xbins, TString addMidpoints)
{
  // cout << "    booking TH1 " << HName.Data() << " nbins=" << Nbins << endl;
  _HistoMapTH1D[HName] =  new TH1D(HName, HTitle, Nbins, xbins);
  _HistoMapTH1D[HName] -> Sumw2();
  if (addMidpoints != "") {
      // cout << HName.Data() << " denser bins: ";
      int Nnewbins = 2*Nbins;
      double *newbins = new double[Nnewbins+1];
      for (int i = 0; i < Nnewbins+1; ++i) {
          if (i % 2 == 0)
              newbins[i] = xbins[i/2];
          else
              newbins[i] = 0.5 * (xbins[(i-1)/2] + xbins[(i-1)/2+1]);
          // cout << newbins[i] << " ";
      }
      // cout << endl;
      // automatically create also a denser (midpoint) version of the non-uniformly binned histo:
      TString newname = HName + cMidPoint;
      _HistoMapTH1D[newname] =  new TH1D(newname, HTitle, Nnewbins, newbins);
      _HistoMapTH1D[newname] -> Sumw2();
      delete [] newbins;
  }
}

void HistoHolder::AddTH2D ( TString HName, TString HTitle,
                Int_t Nbinsx, Double_t xlow, Double_t xhigh,
                Int_t Nbinsy, Double_t ylow, Double_t yhigh)
{
  _HistoMapTH2D[HName] =  new TH2D(HName, HTitle, Nbinsx, xlow, xhigh, Nbinsy, ylow, yhigh);
  _HistoMapTH2D[HName] -> Sumw2();
}

void HistoHolder::AddTH2D ( TString HName, TString HTitle,
                Int_t Nbinsx, Double_t *xbins,
                Int_t Nbinsy, Double_t *ybins)
{
  //   cout << "    booking TH2 C " << HName.Data() << " nbinx=" << Nbinsx << "  nbinsy=" << Nbinsy << endl;
  _HistoMapTH2D[HName] =  new TH2D(HName, HTitle, Nbinsx, xbins, Nbinsy, ybins);
  _HistoMapTH2D[HName] -> Sumw2();
}

void HistoHolder::AddTProfile ( TString HName, TString HTitle, 
				Int_t Nbinsx, Double_t xlow, Double_t xhigh, 
				Double_t ylow, Double_t yhigh)
{
  _HistoMapTProfile[HName] =  new TProfile(HName, HTitle, Nbinsx, xlow, xhigh, ylow, yhigh);
}

void HistoHolder::AddTProfile ( TString HName, TString HTitle, 
				Int_t Nbinsx, Double_t *xbins, 
				Double_t ylow, Double_t yhigh)
{
  _HistoMapTProfile[HName] =  new TProfile(HName, HTitle, Nbinsx, xbins, ylow, yhigh);
}

void HistoHolder::AddTProfile ( TString HName, TString HTitle, 
				Int_t Nbinsx, Double_t xlow, Double_t xhigh)
{
  _HistoMapTProfile[HName] =  new TProfile(HName, HTitle, Nbinsx, xlow, xhigh);
}

void HistoHolder::AddTProfile ( TString HName, TString HTitle, 
				Int_t Nbinsx, Double_t *xbins)
{
  _HistoMapTProfile[HName] =  new TProfile(HName, HTitle, Nbinsx, xbins);
}


// -- filling -- //


void HistoHolder::FillTH1D( TString HName, Double_t val, Double_t weight )
{
  if ( TH1D *h = GetTH1D(HName) )
    h -> Fill(val, weight);
  else
    ThrowError(HName);

  if ( TH1D *h = GetTH1D(HName + cMidPoint) )
    h -> Fill(val, weight);


}



void HistoHolder::FillTH2D( TString HName, Double_t valx, Double_t valy, Double_t weight )
{
  if ( TH2D *h = GetTH2D(HName) )
    h -> Fill(valx, valy, weight);
  else
    ThrowError(HName);
}

void HistoHolder::FillTH3D( TString HName, Double_t valx, Double_t valy, Double_t valz, Double_t weight )
{
  if ( TH3D *h = GetTH3D(HName) )
    h -> Fill(valx, valy, valz, weight);
  else
    ThrowError(HName);
}


void HistoHolder::FillTProfile( TString HName, Double_t valx, Double_t valy, Double_t weight )
{
  if ( TProfile *h = GetTProfile(HName) )
    h -> Fill(valx, valy, weight);
  else
    ThrowError(HName);
}


void HistoHolder::ThrowError(TString HName)
{
  cout << "Error trying to fill histo \"" << HName << "\" while no such created!" << endl;
}


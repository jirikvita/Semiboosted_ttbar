/* File HistoTools.hpp
 *
 * Created       : Fri Sep 14 16:05:31 CDT 2007
 * Author        : kvita
 * Purpose       : 
 * Last modified : 
 * Comments      : 
 */

#ifndef HistoTools_HPP_
#define HistoTools_HPP_

#include "TH1D.h"
#include "TH2D.h"
# include "TProfile.h"

#include "TLine.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraphErrors.h"

// #include "T.h"

#include "TLegend.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TPaveText.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

namespace HistoTools {


  // ___________________________________
  
  TH1D* DivideHistos(TH1D* h1, TF1 *f2, std::string tag = "");
  TH1D* DivideHistos(TH1D* h1, TH1D *h2, std::string tag = "");
  TH1D* DivideHistos(TF1 *f2, TH1D* h1, std::string tag = "");
  TH1D* DivideHistosErrors(TH1D* h1, TH1D *h2, std::string tag = "");
  // ___________________________________


  TString GetSubstring(TString name1, TString name2, TString sub1, TString sub2);
  void PrintCanvasAs(TCanvas *can, std::string suffix = ".eps", bool SaveIndividual = false, std::string name = "");

  bool IsUniformlyBinned(TH1* h1);
  void InvertHisto(TH1D* h1);

  int DivideByBinWidth(TH1* h, int debug = 0);
  int MultiplyByBinWidth(TH1 *h1, int debug = 0);
  int DivideByBinCenter(TH1 *h1, int debug = 0);
  int MultiplyByBinCenter(TH1 *h1, int debug);
  int PrintDivideIntegrals(TH1* h1, TH1* h2);
  int PrintIntegrals(TH1* h1);
  int NormaliseToUnity(TH1* h1, bool UseWidth = false, bool UseOverflows = false);

  void SetHistoStyle(TH1D *h, int col, int width, int mark, float size, int rebin, bool _DivideByBinWidth);
  void SetHistoStyle(TH1D *h, int lcol, int lwidth, int lstyle, int mark, int mcol, float msize);
  void CopyHistoStyle(TH1D *h1, TH1D *h2);
  void CopyHistoStyle(TH1D *h1, TGraph *h2);
  void CopyHistoStyle(TGraph *h1, TH1D *h2);
  int PrintBinContent(TH1* h1);
  int PrintBinContentTeXSimple(TH1* h1, TString outfilename, bool vertical = false, int precision = 3);
  int PrintBinContentTeXSimple2D(TH2* h2, TString outfilename, int precision = 3);
  int PrintBinContentTeX(int nHistos, TH1** h1, TString *Tags, bool *PrintErrorInsteadOfBinContent, TString outfilename, TString caption, TString label, 
			 int PrintError, bool PrintSum, bool UseWidth, bool UseOverflows, bool SkipEmptyBins, 	
			 int precision = 3, int error_precision = 2, double SubtractVal = 0.,	
			 TH1D *ScaleHisto = 0, TH1D *CentralHisto = 0, bool Symmetrise = true);
  void SetErrorBarsRelative(TH1D *from, TH1D* to);
  void AddErrorBarsRelative(TH1D *from, double CentralConstant, TH1D* to, TH1D *scale = 0);
  void ScaleErrorBars(TH1D* to, double errorSF);

  TH1D *MakeFitResidualsHisto(TH1D* h1, TF1 *fit, std::string tag = "_res", 
			      bool _relative = true, float xmin = -10, float xmax = -100 );
  TH1D *MakeFitResidualsProjection(TH1D* h1, TF1 *fit, std::string tag = "_res", bool _relative = true,
				   float x1 = -3., float x2 = 3., int Nbins = 100, bool weightBySigma = false,
				   float xmin = -10, float xmax = -100);

  TH2D* MergeEnsembleHistos(int Nesb, TH1D **hist, int nbinsy, std::string name = "", float maxmaxy = 999, int ibin0 = 0);

  bool CheckFileOpen(TFile *f);
 
  
  Double_t GetNeff(TH1D* h1);
  Double_t GetMCError(TH1D* h1, Float_t *N_eff);

  Double_t HistMeanError (TH1* h, Int_t axis=1);
  Double_t HistRMSError (TH1* h, Int_t axis=1);
  Float_t GetYError(TH1D* dist, Float_t *N_eff);

  void PrintTF1Pars(TF1 *f);

  // ___________________________________

  bool FitGauss(TH1D* histo, int rebin, float width, float* mean, float* sigma, float* error, const char* fitoption,
		float _MinimumEntries, bool getfit, TF1 **fit);
  
  int PlotProfile(TH2 *h2, 
		  TGraphErrors **profile_graph, TGraphErrors **resolution_graph, 
		  TPad **profile_can, TPad **resolution_can,
		  bool doGauss, bool plotSlices, 
		  TString title = "graph", TString SlicesTag = "",
		  bool zoom = false, Float_t x1 = -1., Float_t x2 = 1., 
		  bool draw = true, 
		  float _MinimumEntriesForGaussFit = 20., float widthToFitGaussAround = 2.,
		  float graph_x1 = -1., float graph_x2 = 1.,
		  float graph_y1 = -1., float graph_y2 = 1.,
		  int col     =  2,
		  int width   =  1,
		  int mark    = 22,
		  float size  =  0.7);
 
  int GetH2AndPlotProfile(TString hname, TString filename, TString dirname, 
			  TGraphErrors **profile_graph, 
			  TGraphErrors **resolution_graph,
			  TPad **profile_can, TPad **resolution_can, 
			  bool doGauss, bool plotSlices, TString title, TString SlicesTag,
			  bool zoom , Float_t x1, Float_t x2,
			  bool draw, bool getH2, TH2D **h2copy, TFile **file,
			  float _MinimumEntriesForGaussFit, float widthToFitGaussAround,
			  float graph_x1, float graph_x2,
			  float graph_y1 = 1, float graph_y2 = -1);

  TH1D* MakeHistoFromGraph(TGraphErrors *graph, TH1D* BinningHisto, TString name);
  TH1D* MakeHalfMaxDiffHisto(std::vector<TH1D*> histos, std::string name);
  std::vector<std::string> TokenizeString(std::string input, std::string delim);
  TH1D* MakeSymErrorBand(TH1D *center, TH1D *ErrUp, TString name, bool reset = false, double center_value = 0.);
  TH1D* MakeSymErrorBandHisto(TH1D *center, TH1D *Err, float Sign, TString name, bool reset, double center_value);

  void ZeroErrorBars(TH1D *h);
  void SquareHisto(TH1D *h);
  void SquareRootHisto(TH1D *h);
  void AddInSquares(TH1D *AddTo, TH1D *what, float OverAllSign);

  void Transpose(TH2 *histo);
  TH1D *MakeDiagonalHisto(TH2D *h2, TString name);

  TLegend *MakeLegend( float x1, float y1, float x2, float y2);
  TPaveText* MakePaveText(TString item, float x1 = 0.50, float y1 = 0.74, float x2 = 0.85, float y2 = 0.90, float size = 0.04);
  // ___________________________________

  TH1D *MakePoissonSmearedHisto(TH1D *h, TRandom3 *rand, std::string tag = "_smeared");
  TH1D* MakeConstHisto(TH1D *h, float CentralConst, float Const, std::string name, bool MultiplyByHisto = false);

  // ___________________________________

  TH1D* GetHisto1DFromFile(TFile *file, std::string dirname, std::string hname);
  TH2D* GetHisto2DFromFile(TFile *file, std::string dirname, std::string hname);
  TObject* GetFromFile(TFile *file, std::string dirname, std::string hname);
  TObject* GetFromFile(TFile *file, std::string hname); 
  TObject* GetFromFile(std::string filename, TFile *file, std::string hname);

  TH1D *TruncateBinsFromHisto(TH1D *hist, TString name, int binout1, int binout2);
  TGraphErrors* MakeGraphFromHisto(TH1D *histo, TString name, bool ignoreYerrors = true);

  // july 2009:
  int AddHistosAver(TH1D *from, TH1D* to, double AverFactor = 2.);
  TGraphErrors* AddGraphsAver(TGraphErrors *from, TGraphErrors *to, double AverFactor = 2., TString MergeTag = "_merge");
  int PrintGraphContent(TGraphErrors *graph);

  int AbsHisto(TH1D *hist);

// ____________________________________________

} // namespace

#endif

 

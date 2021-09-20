#ifndef KINEM_UTILS
#define KINEM_UTILS

#include "TMyLorentzVector.h"

#define sqr(X) pow(X, 2)

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

bool cmpPt(const TMyLorentzVector *a, const TMyLorentzVector *b);
void printJetArray(vector<TMyLorentzVector*> jets, TString tag = "");

struct mindrdetadphi
{
    double mindr;
    double mindeta;
    double mindphi;

};

bool MatchedWithinDR(TMyLorentzVector *vec, vector<TMyLorentzVector*> lvecs, mindrdetadphi &vals, double DRcut = 0.5);

double computeneutrinoz(const TMyLorentzVector lepton, double metx, double mety, bool invertSolutions, double& discriminant);
double computeneutrinoz_from_mt(double mthad, const TMyLorentzVector lepton, const TMyLorentzVector blep, double metx, double mety, bool invertSolutions, double& discriminant);
double mT(TMyLorentzVector l, TMyLorentzVector nu);

#endif

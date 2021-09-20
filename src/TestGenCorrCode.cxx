// jk 6.8.2020

#include "TH2D.h"
#include "TString.h"

#include "HistoHolder.h"

#include <iostream>
#include <map>
using namespace std;

TString GeV = " [GeV]";

class varProps {
public:
    varProps(TString tit, int i, TString un, TString fillstr) {
        this -> title = tit;
        this -> id = i;
        this -> unit = un;
        this -> fillstr = fillstr;
    };
    varProps() {};
    ~varProps() {};
    TString title;
    int id;
    TString unit;
    TString fillstr;
};

std::map<TString,varProps> m_vars;


int main()
{

    int ivar = 0;
    m_vars["DiTopPout"] = varProps("|p_{out}|", ivar++, GeV, "pseudo.Pout[0]");
    m_vars["DiTopDeltaPhi"] = varProps("#Delta#phi", ivar++, "", "pseudo.DeltaPhi");
    m_vars["DiTopMass"] = varProps("m^{t#bar{t}}", ivar++, GeV, "pseudo.pseudottbar.M()");
    m_vars["DiTopPt"] = varProps("p^{t#bar{t}}_{T}", ivar++, GeV, "pseudo.pseudottbar.Pt()");
    m_vars["TopPt"] = varProps("p^{t}_{T}", ivar++, GeV, "pseudo.pseudotophadron.Pt()");
    m_vars["CosThetaStar"] = varProps("|cos#theta*|", ivar++, "", "pseudo.CosThetaStar[0]");
    m_vars["Delta"] = varProps("#delta^{t#bar{t}}", ivar++, "", "pseudo.deltattbar");
    m_vars["Chittbar"] = varProps("#chi^{t#bar{t}}", ivar++, "", "pseudo.Chittbar");
    m_vars["Yboost"] = varProps("y_{boost}^{t#bar{t}}", ivar++, "", "pseudo.Yboost");
    m_vars["Rttbar"] = varProps("R^{t2,t1}", ivar++, "", "pseudo.Rttbar");

    m_vars["jApla"] = varProps("jets Aplanarity", ivar++, "", "GlobalVars.jApla");
    m_vars["jSpher"] = varProps("jets Sphericity", ivar++, "", "GlobalVars.jSpher");
    //m_vars["JApla"] = varProps("Jets Aplanarity", ivar++, "", "GlobalVars.JApla");
    //m_vars["JSpher"] = varProps("Jets Sphericity", ivar++, "", "GlobalVars.JSpher");
    m_vars["HTj"] = varProps("H_{T}^{j}", ivar++, GeV, "GlobalVars.HTj");
    //m_vars["SumMj"] = varProps("#sum m^{j}", ivar++, GeV, "GlobalVars.SumMj");
    //m_vars["HTJ"] = varProps("H_{T}^{J}", ivar++, GeV, "GlobalVars.HTJ");
    m_vars["HTjPlusMet"] = varProps("H_{T}^{j} + E^{miss}_{T}", ivar++, GeV, "GlobalVars.HTj + GlobalVars.met");
    //m_vars["HTJPlusMet"] = varProps("H_{T}^{J} + E^{miss}_{T}", ivar++, GeV, "GlobalVars.HTJ + GlobalVars.met");
    m_vars["SumMJ"] = varProps("#sum m^{J}", ivar++, GeV, "GlobalVars.SumMJ");

    m_vars["Met"] = varProps("E^{miss}_{T}", ivar++, GeV, "GlobalVars.met");
    //m_vars["MjVis"] = varProps("m_{sum j}^{vis}", ivar++, GeV, "GlobalVars.MjVis");
    m_vars["MJVis"] = varProps("m_{sum J}^{vis}", ivar++, GeV, "GlobalVars.MJVis");

    m_vars["DiTopPtRel"] = varProps("p_{T}^{t#bar{t}} / m^{t#bar{t}}", ivar++, "", "pseudo.DiTopPtRel");
    m_vars["DiTopPtGeo"] = varProps("p_{T}^{t#bar{t}} / #sqrt{p_{T}^{t1} p_{T}^{t2}}", ivar++, "", "pseudo.DiTopPtGeo");
    m_vars["DiTopMassGeo"] = varProps("m^{t#bar{t}} / #sqrt{p_{T}^{t1} p_{T}^{t2}}", ivar++, "", "pseudo.DiTopMassGeo");
    m_vars["DiTopPoutRel"] = varProps("|p_{out}| / m^{t#bar{t}}", ivar++, "", "pseudo.DiTopPoutRel[0]");
    m_vars["DiTopPoutGeo"] = varProps("|p_{out}| / #sqrt{p_{T}^{t1} p_{T}^{t2}}", ivar++, "", "pseudo.DiTopPoutGeo[0]");
    m_vars["TopPtRel"] = varProps("p_{T}^{t} / m^{t#bar{t}}", ivar++, "", "pseudo.TopPtRel[0]");



    HistoHolder *hold = new HistoHolder("test");


    map<TString, double*> m_PhysObjBins;
    map<TString, int> m_nPhysObjBins;

    cout << "----------------------------------------------------------------------------" << endl;
    // test code for real usage in HistoMaker.cpp
    // https://stackoverflow.com/questions/6963894/how-to-use-range-based-for-loop-with-stdmap
    for (auto& [var1, varProp1]: m_vars) {
        for (auto& [var2, varProp2]: m_vars) {
            if (var1 == var2)
                continue;
            if (varProp1.id > varProp2.id)
                continue;
            /*
          if (var1.Contains("Apla") && var2.Contains("Apla"))
              continue;
          if (var1.Contains("HT") && var2.Contains("HT"))
              continue;
          if (var1.Contains("Spher") && var2.Contains("Spher"))
              continue;
          */
            // so far fragile but works;-)
            // don't care about correlations between jet and LJet related vars
            //if (var1.Contains("j") && var2.Contains("J"))
            //    continue;

            TString name = var2 + "Vs" + var1;
            TString tag = " + reptag";
            /*
      hold -> AddTH2D(name + tag,
              name + ";" + varProp1.title + varProp1.unit + " ;" + varProp2.title  + varProp2.unit,
              m_nPhysObjBins[var1], m_PhysObjBins[var1], m_nPhysObjBins[var2], m_PhysObjBins[var2]);
              */
        }
    }
    // info
    // cout << var1.Data() << " " << var2.Data() << endl;

    cout << "----------------------------------------------------------------------------" << endl;

    for (auto& [var1, varProp1]: m_vars) {
        for (auto& [var2, varProp2]: m_vars) {
            if (var1 == var2)
                continue;
            if (varProp1.id > varProp2.id)
                continue;
            /* jk 30.4.2021
            if (var1.Contains("j") && var2.Contains("J"))
                continue;
            if (var1.Contains("J") && var2.Contains("j"))
                continue;
            */

            TString name = var2 + "Vs" + var1;
            TString tag = " + reptag";

            // prepare filling code
            TString fillname = "m_level + \"" + var2 + "\" + \"Vs\" + m_level + \"" + var1 + "\"" + tag;

            // JK 23.3.2021 add here weights of 0.5 fro pseudotophadron and pseudotoplepton pTs!
            // otherwise leads to different DiTopMass significances as projected from
            // singly or doubly filled histograms!

            bool var1isTwoFold = varProp1.fillstr.Contains("pseudotophadron") || varProp1.fillstr.Contains("[0]");
            bool var2isTwoFold = varProp2.fillstr.Contains("pseudotophadron") || varProp2.fillstr.Contains("[0]");
            if (var1isTwoFold && var2isTwoFold) {
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << varProp1.fillstr << " , " << varProp2.fillstr <<  ", rweight*0.25);" << endl;
            } else if (var1isTwoFold || var2isTwoFold) {
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << varProp1.fillstr << " , " << varProp2.fillstr <<  ", rweight*0.5);" << endl;
            } else {
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << varProp1.fillstr << " , " << varProp2.fillstr <<  ", rweight);" << endl;
            }

            if (varProp1.fillstr.Contains("pseudotophadron")) {
                TString fillstr1 = varProp1.fillstr;
                TString wstr = "0.5";
                if (var2isTwoFold)
                    wstr = "0.25";
                fillstr1.ReplaceAll("pseudotophadron","pseudotoplepton");
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << fillstr1 << " , " << varProp2.fillstr <<  ", rweight*" << wstr.Data() << "); // add1" << endl;
            }
            if (varProp2.fillstr.Contains("pseudotophadron")) {
                TString fillstr2 = varProp2.fillstr;
                fillstr2.ReplaceAll("pseudotophadron","pseudotoplepton");
                TString wstr = "0.5";
                if (var1isTwoFold)
                    wstr = "0.25";
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << varProp1.fillstr << " , " << fillstr2 <<  ", rweight*" << wstr.Data() << "); // add1" << endl;
            }


            // [0] <-> [1]
            if (varProp1.fillstr.Contains("[0]") && !varProp2.fillstr.Contains("[0]") ) {
                TString fillstr1 = varProp1.fillstr;
                fillstr1.ReplaceAll("[0]","[1]");
                TString wstr = "0.5";
                if (varProp2.fillstr.Contains("pseudotophadron"))
                    wstr = "0.25";
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << fillstr1 << " , " << varProp2.fillstr <<  ", rweight*" << wstr.Data() << "); // add2 " << endl;
                if (varProp2.fillstr.Contains("pseudotophadron")) {
                    TString fillstr2 = varProp2.fillstr;
                    fillstr2.ReplaceAll("pseudotophadron","pseudotoplepton");
                    cout << "hold -> FillTH2D(" << fillname.Data() << ", " << fillstr1 << " , " << fillstr2 <<  ", rweight*" << wstr.Data() << "); // add2 " << endl;
                }
            } else if (varProp2.fillstr.Contains("[0]") && !varProp1.fillstr.Contains("[0]") ) {
                TString fillstr2 = varProp2.fillstr;
                fillstr2.ReplaceAll("[0]","[1]");
                TString wstr = "0.5";
                if (varProp1.fillstr.Contains("pseudotophadron"))
                    wstr = "0.25";
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << varProp1.fillstr << " , " << fillstr2 <<  ", rweight*" << wstr.Data() << "); // add3 " << endl;
                if (varProp1.fillstr.Contains("pseudotophadron")) {
                    TString fillstr1 = varProp1.fillstr;
                    fillstr1.ReplaceAll("pseudotophadron","pseudotoplepton");
                    cout << "hold -> FillTH2D(" << fillname.Data() << ", " << fillstr1 << " , " << fillstr2 <<  ", rweight*" << wstr.Data() << "); // add3 " << endl;
                }
            } else if (varProp1.fillstr.Contains("[0]") && varProp2.fillstr.Contains("[0]") ) {
                TString fillstr1 = varProp1.fillstr;
                TString fillstr2 = varProp2.fillstr;
                fillstr1.ReplaceAll("[0]","[1]");
                fillstr2.ReplaceAll("[0]","[1]");
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << fillstr1 << " , " << fillstr2 <<  ", rweight*0.25); // add4" << endl;
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << varProp1.fillstr << " , " << fillstr2 <<  ", rweight*0.25); // add5" << endl;
                cout << "hold -> FillTH2D(" << fillname.Data() << ", " << fillstr1 << " , " << varProp1.fillstr <<  ", rweight*0.25); // add6" << endl;

            }



        }

    }





    cout << "----------------------------------------------------------------------------" << endl;
    cout << "ivarsDict = {" << endl;
    for (auto& [var, varProp]: m_vars) {
        cout << "  '" << var << "' : " << varProp.id << "," << endl;
    }
    cout << "}" << endl;


    cout << "----------------------------------------------------------------------------" << endl;
    cout << "ivarsLabelsDict = {" << endl;
    for (auto& [var, varProp]: m_vars) {
        cout << "  '" << var << "' : '" << varProp.title.Data() << "'," << endl;
    }
    cout << "}" << endl;


    cout << "----------------------------------------------------------------------------" << endl;
    cout << "hnamesDict = {  " << endl;
    for (auto& [var1, varProp1]: m_vars) {
        for (auto& [var2, varProp2]: m_vars) {
            if (var1 == var2)
                continue;
            if (varProp1.id > varProp2.id)
                continue;
            //if (var1.Contains("j") && var2.Contains("J"))
            //    continue;

            TString name = var2 + "Vs" + var1;
            cout << "  '" << name.Data() << "' : " << "[" << varProp1.id << "," << varProp2.id << ", True]," << endl;
        }
    }
    cout << "}" << endl;

    cout << "----------------------------------------------------------------------------" << endl;

    cout << "allBH2vars = [";
    for (auto& [var1, varProp1]: m_vars) {
        for (auto& [var2, varProp2]: m_vars) {
            if (var1 == var2)
                continue;
            if (varProp1.id > varProp2.id)
                continue;
            //if (var1.Contains("j") && var2.Contains("J"))
            //    continue;

            TString name = var2 + "Vs" + var1;
            cout << "'" << name << "', ";
        }
    }
    cout << " ]" << endl;


} // macro 


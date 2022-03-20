#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <stdio.h>
#include <cstdio>

#include <TCanvas.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TStopwatch.h>

//#include "tools.h"

using namespace std;

vector <TH1F*> make_histo(string name){
    ifstream myfile;
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
    vector <float> t(1030);
    vector <vector<float>> v;
    vector <TH1F*> histos;
    iota(begin(t), end(t), 0);


    if (myfile.is_open()){
        string tp;
        vector <float> w;
        int k=0, m, j=0;
        while(getline(myfile, tp)){ //read data from file object and put it into string.
            if (isalpha(tp[0]) == 0) {
                if (k>0) {
                    k=0;
                    w.clear();
                }
                m = stof(tp);
                w.push_back(m);
                j++;
            }
            else {
                if (j>0) {
                    j=0;
                    v.push_back(w);
                }                
                k++;
            }
        }
        //non prendo l'ultimo evento perchè è incompleto
        myfile.close(); //close the file object.
    }    
    
    vector <float> charge (v.size());
    vector <float> amp (v.size());
    float h, u, min, charge_max=FLT_MIN, amp_max=FLT_MIN, charge_min=FLT_MAX, amp_min=FLT_MAX;
    vector <float> idx_strange;
    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        u=0;
        for (int j=0; j<100; j++){
            h += v[i][j] / 100;
            if (v[i][j]<14600 && u==0){
                //cout << "valore piccolo: "<< v[i][j] <<endl;
                idx_strange.push_back(i);
                u++;
            }
        }
        for (int j=0; j<1000; j++){
            charge[i] += h-v[i][j];
        }
        min =*min_element(v[i].begin(), v[i].end());
        amp[i]=h-min;
        charge_max = (charge[i]>charge_max) ? charge[i] : charge_max;
        amp_max = (amp[i]>amp_max) ? amp[i] : amp_max;
        charge_min = (charge[i]<charge_min) ? charge[i] : charge_min;
        amp_min = (amp[i]<amp_min) ? amp[i] : amp_min;      
    }
    // if pmt1 allora range_max = ....
    
    TH1F *histo_charge=  new TH1F(&(name + "_charge")[0],&(name + "_charge")[0], 1300, 0, /*static_cast<int>(charge_max)-150000*/ 250000);
    TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 1400, 0, /*static_cast<int>(amp_max)-4500*/ 3500);
    for (long unsigned int i=0; i< charge.size(); i++){
        if (find(idx_strange.begin(), idx_strange.end(), i) == idx_strange.end()){
            histo_charge->Fill(charge[i]);
            histo_amp->Fill(amp[i]);
        }
    }
    histo_charge->SetName(&(name + "_charge")[0]);
    histo_amp->SetName(&(name + "_amp")[0]);
    histos.push_back(histo_charge);
    histos.push_back(histo_amp);

    TCanvas *c_charge = new TCanvas(&(name + "_charge")[0] , &(name + "_charge")[0]);
    histo_charge->Draw();
    c_charge->SaveAs(&("plots/" + name + "_charge.png")[0]);

    TCanvas *c_amp = new TCanvas(&(name + "_amp")[0], &(name + "_amp")[0]);
    histo_amp->Draw();
    c_amp->SaveAs(&("plots/" + name + "_amp.png")[0]);

    cout << "numero di idx_strange  "<<idx_strange.size() <<endl;

    TCanvas *c_wave = new TCanvas(&(name + "_wave")[0], &(name + "_wave")[0]);
    TGraph* gr = new TGraph(t.size(), &t[0], &v[0][0]);
    gr->SetNameTitle(&(name + "_wave")[0], &(name + "_wave")[0]);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");
    c_wave->SaveAs(&("waveform/" + name + "_wave.png")[0]);

    return histos;
}

void main_func (){
    gROOT->SetBatch(kFALSE);
    vector <TH1F*> histos;
    TFile *outfile= new TFile("histograms/histograms_RealTimeCalibration.root", "RECREATE"/* "UPDATE"*/);
    list <string> names ={/*"pmt1_NA_e6_100_or_run1",
                        "pmt2_NA_e6_100_or_run1",
                        "pmt1_NA_e6_700_or_run2",
                        "pmt2_NA_e6_700_or_run2",
                        "pmt1_NA_e6_700_or_run3",
                        "pmt2_NA_e6_700_or_run3",
                        "pmt1_NA_e6_ext_run1",
                        "pmt2_NA_e6_ext_run1",
                        "pmt1_NA_e6_ext_run2",
                        "pmt2_NA_e6_ext_run2",
                        "pmt1_NA_e6_ext_run3",
                        "pmt3_NA_e6_ext_run3",
                        "pmt1_NA_e6_100_run1",
                        "pmt2_NA_e6_100_run1",
                        "pmt3_NA_e6_100_run1",
                        "pmt1_NA_e6_100_run2",
                        "pmt2_NA_e6_100_run2",
                        "pmt3_NA_e6_30_run2",
                        "pmt1_co_100", 
                        "pmt2_co_100",
                        "pmt3_co_100_run1",
                        "pmt3_co_30_run2",
                        "pmt1_na_100", 
                        "pmt2_na_100", 
                        "pmt1_cs_100", 
                        "pmt2_cs_100",
                        "pmt3_cs_100_run1",
                        "pmt3_cs_30_run2",
                        "pmt1_bkg_100", 
                        "pmt2_bkg_100", 
                        "pmt1_null", 
                        "pmt2_null",
                        "pmt1_NA_cs2_e6_ext_run1",
                        "pmt2_NA_cs2_e6_ext_run1",
                        "pmt1_NA_cs1_e6_ext_solo_run1",
                        "pmt2_NA_cs1_co1_e6_ext_run1", 
                        "pmt1_NA_cs1_co1_e6_ext_run1", */
                        "pmt1_NA_cs_e6_100_run1", 
                        "pmt2_NA_cs_e6_100_run1",
                        "pmt3_NA_cs_e6_30_run1",
                        "pmt1_NA_cs_co_e6_100_run1",
                        "pmt2_NA_cs_co_e6_100_run1",
                        "pmt3_NA_cs_co_e6_30_run1", 
                        "pmt1_NA_cs_e6_500_run2", 
                        "pmt2_NA_cs_e6_500_run2",
                        "pmt3_NA_cs_e6_100_run2",
                        "pmt1_NA_cs_co_e6_500_run2", 
                        "pmt2_NA_cs_co_e6_500_run2",
                        "pmt3_NA_cs_co_e6_100_run2"};

    TStopwatch time_tot;
    time_tot.Start();                
    for(list<string>::const_iterator name = names.begin(); name != names.end(); ++name){
        TStopwatch time;
        time.Start();
        cout << "Processing: " << *name << endl; 
        histos=make_histo(*name);
        histos[0]->Write();
        histos[1]->Write();
        time.Stop();
        time.Print();
    }
    //getchar();
    outfile->Close();
    time_tot.Stop();
    time_tot.Print();
}

int main() {
    main_func();
}



// execute with:
// $ root -l
// $ .L tools.cpp
// $ .x main_func.cpp


// execute with:
// $ root -l
// $.L tools.cpp+
// $ root -l tools_cpp.so main_func.cpp

// $ g++ -Wall -Wextra -Werror -pedantic -std=c++14 main_func.cpp tools.cpp  `root-config --glibs --cflags` -o main_func
// ./main_func

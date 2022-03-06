#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <float.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>

#include "tools.h"

using namespace std;

vector <TH1F*> make_histo(string path,string name){
    ifstream myfile;
    myfile.open(path, ios::in | ios::out);  
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
    float h, min, charge_max=FLT_MIN, amp_max=FLT_MIN, charge_min=FLT_MAX, amp_min=FLT_MAX;

    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        for (int j=0; j<200; j++){
            h += v[i][j] / 200;
        }
        for (int j=400; j<900; j++){
            charge[i] += h-v[i][j];
        }
        min =*min_element(v[i].begin(), v[i].end());
        amp[i]=h-min;
        charge_max = (charge[i]>charge_max) ? charge[i] : charge_max;
        amp_max = (amp[i]>amp_max) ? amp[i] : amp_max;
        charge_min = (charge[i]<charge_min) ? charge[i] : charge_min;
        amp_min = (amp[i]<amp_min) ? amp[i] : amp_min;        
    }
    TH1F *histo_charge=  new TH1F(&(name + "_charge")[0],&(name + "_charge")[0], 7000, -1000, static_cast<int>(charge_max));
    TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 700, 0, static_cast<int>(amp_max));
    for (long unsigned int i=0; i< charge.size(); i++){
        histo_charge->Fill(charge[i]);
        histo_amp->Fill(amp[i]);
    }
    histo_charge->SetName(&(name + "_charge")[0]);
    histo_amp->SetName(&(name + "_amp")[0]);
    histos.push_back(histo_charge);
    histos.push_back(histo_amp);

    TCanvas *c_wave = new TCanvas(&(name + "_wave")[0], &(name + "_wave")[0]);
    TGraph* gr = new TGraph(t.size(), &t[0], &v[1][0]);
    gr->SetNameTitle(&(name + "_wave")[0], &(name + "_wave")[0]);
    gr->Draw("AP*");
    c_wave->SaveAs(&("waveform/" + name + "_wave.png")[0]);

    TCanvas *c_charge = new TCanvas(&(name + "_charge")[0] , &(name + "_charge")[0]);
    histo_charge->Draw();
    c_charge->SaveAs(&("plots/" + name + "_charge.png")[0]);

    TCanvas *c_amp = new TCanvas(&(name + "_amp")[0], &(name + "_amp")[0]);
    histo_amp->Draw();
    c_amp->SaveAs(&("plots/" + name + "_amp.png")[0]);
    
    return histos;
}





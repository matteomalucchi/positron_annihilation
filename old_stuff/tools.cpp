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
    float h, min,min_idx,min_h,  charge_max=FLT_MIN, amp_max=FLT_MIN, charge_min=FLT_MAX, amp_min=FLT_MAX;
    vector <float> idx_strange;
    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        for (int j=0; j<100; j++){
            h += v[i][j] / 200;
            if (v[i][j]<14600){
                cout << "valore piccolo: "<< v[i][j] <<endl;
                idx_strange.push_back(i);
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
        if (charge[i]<-400000){
            min_idx=i;
            min_h=h;
        }
        amp_min = (amp[i]<amp_min) ? amp[i] : amp_min;      
    }
    TH1F *histo_charge=  new TH1F(&(name + "_charge")[0],&(name + "_charge")[0], 1300, static_cast<int>(charge_min), static_cast<int>(charge_max));
    TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 1400,  static_cast<int>(amp_min), static_cast<int>(amp_max));
    for (long unsigned int i=0; i< charge.size(); i++){
        if (find(idx_strange.begin(), idx_strange.end(), i)== idx_strange.end()){
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

    cout << "charge min "<<charge_min <<endl;
        cout << "charge min h "<<min_h <<endl;

    TCanvas *c_wave = new TCanvas(&(name + "_wave")[0], &(name + "_wave")[0]);
    TGraph* gr = new TGraph(t.size(), &t[0], &v[min_idx][0]);
    gr->SetNameTitle(&(name + "_wave")[0], &(name + "_wave")[0]);
    gr->Draw("AP*");
    c_wave->SaveAs(&("waveform/" + name + "_wave.png")[0]);

    return histos;
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


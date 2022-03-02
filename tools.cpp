#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <TCanvas.h>
#include <TH1F.h>

#include "tools.h"


using namespace std;

vector <TH1F*> make_histo(string path,string name){
    ifstream myfile;
    myfile.open(path, ios::in | ios::out);  
    vector <float> t(1030);
    vector <vector<float>> v;
    vector <TH1F*> histos;
    TH1F *histo_charge=  new TH1F(&(name + "_charge")[0],&(name + "_charge")[0], 800, 0, 200000);
    TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 800, 0, 15000);
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
        //non prendo l'ultimo evetno perchè è incompleto
        myfile.close(); //close the file object.
    }    
    
    vector <float> charge (v.size());
    vector <float> amp (v.size());
    float h, min;

    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        for (int j=0; j<200; j++){
            h += v[i][j] / 200;
        }
        for (int j=400; j<900; j++){
            charge[i] += abs(h-v[i][j]);
        }
        min =*min_element(v[i].begin(), v[i].end());
        amp[i]=abs(h-min);
    }
    for (long unsigned int i=0; i< charge.size(); i++){
        histo_charge->Fill(charge[i]);
        histo_amp->Fill(amp[i]);
    }
    histos.push_back(histo_charge);
    histos.push_back(histo_amp);

    return histos;
}



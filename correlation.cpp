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

vector<float> find_idx_range(string path,string name, vector<float> idx_corr ={}){
    ifstream myfile;
    myfile.open(path, ios::in | ios::out);  
    vector <vector<float>> v;

    if (myfile.is_open()){
        string tp;
        vector <float> w;
        int k=0, m, j=0, bgn=0;
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
    float h;
    vector <float> idx_range;
    float num_events=0;
    if (idx_corr.empty()){
        for( int i = 0; i < v.size(); i++ ){
        idx_corr.push_back( i );        
        }
    } 
    for (float &i : idx_corr){
        h = 0;
        for (int j=0; j<100; j++){
            h += v[i][j] / 100;
        }
        for (int j=0; j<1000; j++){
            charge[i] += h-v[i][j];
        }
        if (name.find("pmt1_NA_e6_ext_run1")<name.length()){
            if (charge[i] > 40000) {
                num_events++;
                if (charge[i]> 73000 && charge[i]< 79000) idx_range.push_back(i);
                //if (charge[i]> 100000 && charge[i]< 200000) idx_range.push_back(i);

            }
        }      
        if (name.find("pmt2_NA_e6_100")<name.length()){
                if (charge[i] > 26000) {
                    num_events++;
                    if (charge[i]> 70000 && charge[i]< 75000) idx_range.push_back(i);
                }
        }
        if (name.find("pmt2_NA_e6_ext_run1")<name.length()){
                if (charge[i] > 26000) {
                    num_events++;
                    if (charge[i]> 65500 && charge[i]< 70500) idx_range.push_back(i);
                }
        }
    }
    idx_range.push_back(num_events);

    return idx_range;
}

void correlation (){
    gROOT->SetBatch(kFALSE);
    vector<vector<float>> idx_ranges;
    ofstream out_file("correlation_prob_1a2a.txt");

    list <string> names ={
                        //"pmt1_NA_e6_100_run1",
                        "pmt2_NA_e6_100_run1",            
                        "pmt1_NA_e6_ext_run1",
                        "pmt2_NA_e6_ext_run1"
                        };

    for(list<string>::const_iterator name = names.begin(); name != names.end(); ++name){
        cout << "Processing: " << *name << endl; 
        idx_ranges.push_back(find_idx_range("data/" + *name + ".txt", *name));
        float a =  (idx_ranges.back().size()-1)/idx_ranges.back().back();
        cout << a <<endl;
        cout << idx_ranges.back().back() <<endl;
        string prob="prob singola picco:     "+ *name + " = " + to_string(a) + "\n";
        string events = "eventi considerati:    " +*name+ " = " + to_string(idx_ranges.back().back()) + "\n\n";
        out_file<< prob;
        out_file << events;
    }
    cout << "Processing: final " << "pmt2_NA_e6_ext_run1" << endl;
    idx_ranges.push_back(find_idx_range("data/pmt2_NA_e6_ext_run1.txt", "pmt2_NA_e6_ext_run1", idx_ranges[1]));
    float a =  (idx_ranges.back().size()-1)/idx_ranges.back().back();
    cout << a <<endl;
    cout << idx_ranges.back().back() <<endl;
    
    string prob= "prob condizionata picco 2 a picco 1 = " + to_string(a) + "\n";
    string events = "eventi considerati = " + to_string(idx_ranges.back().back()) + "\n\n";
    out_file<< prob;
    out_file << events;

}

int main() {
    correlation();
}


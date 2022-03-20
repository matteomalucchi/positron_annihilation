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

vector<float> find_idx_range(string path,string name, vector<float> idx_corr ={}, int &time_bgn, int &time_end){
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
            if (tp.find("Time")<tp.length()){
                if (bgn==0) {
                    time_bgn=stof(tp.substr(20,tp.size()-1));
                    cout << time_bgn <<endl;
                    bgn++;
                }
                time_end=stof(tp.substr(20,tp.size()-1));
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

    cout << time_end <<endl;
    
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
        if (name.find("pmt1_NA_e6_or_run1")<name.length()){
            if (charge[i] > 40000) {
                picco1++;
                if (charge[i]> 73000 && charge[i]< 79000) idx_range.push_back(i);
                //if (charge[i]> 100000 && charge[i]< 200000) idx_range.push_back(i);

            }
        }      
        if (name.find("pmt2_NA_e6_or_run1")<name.length()){
            if (charge[i] > 26000) {
                picco2++;
                if (charge[i]> 70000 && charge[i]< 75000) idx_range.push_back(i);
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
                        //"pmt1_NA_e6_100",
                        "pmt2_NA_e6_100",            
                        "pmt1_NA_e6_ext_run1",
                        "pmt2_NA_e6_ext_run1"
                        };

    int time_bgn=0, time_end=0;
    for(list<string>::const_iterator name = names.begin(); name != names.end(); ++name){
        cout << "Processing: " << *name << endl; 
        idx_ranges.push_back(find_idx_range("data/" + *name + ".txt", *name, time_bgn, time_end));
        float a =  (idx_ranges.back().size()-1)/idx_ranges.back().back();
        cout<< "Frequenza temporale generale: "<< idx_ranges.back().back()/(time_end-time_bgn)*8*pow(10,-9) <<endl;
        cout<< "Frequenza temporale nel picco: "<<  (idx_ranges.back().size()-1)/(time_end-time_bgn)*8*pow(10,-9) <<endl;
        cout << a <<endl;
        cout << idx_ranges.back().back() <<endl;
        string prob= *name + " = " + to_string(a) + "\n";
        string events = "eventi considerati = " + to_string(idx_ranges.back().back()) + "\n\n";
        out_file<< prob;
        out_file << events;
    }
    cout << "Processing: final " << "pmt2_NA_e6_ext_run1" << endl;
    idx_ranges.push_back(find_idx_range("data/pmt2_NA_e6_ext_run1.txt", "pmt2_NA_e6_ext_run1", idx_ranges[1], time_bgn, time_end));
    float a =  (idx_ranges.back().size()-1)/idx_ranges.back().back();
    cout << a <<endl;
    cout << idx_ranges.back().back() <<endl;
    
    string prob= "pmt2_NA_e6_ext_run1 = " + to_string(a) + "\n";
    string events = "eventi considerati = " + to_string(idx_ranges.back().back()) + "\n\n";
    out_file<< prob;
    out_file << events;

}

int main() {
    correlation();
}

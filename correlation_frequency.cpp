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

vector<float> find_idx_range(string path,string name, long int &time_tot){
    ifstream myfile;
    myfile.open(path, ios::in | ios::out);  
    vector <vector<float>> v;
    long int time_before=0, time_now=0, time_bgn=0, time_min=pow(2,32);
    int g=0, a=0;
    if (myfile.is_open()){
        string tp;
        vector <float> w;
        int k=0, j=0, bgn=0, b=0;
        float m;
        while(getline(myfile, tp) /*&& b<1000000*/){ //read data from file object and put it into string.
           // cout <<b <<endl;
            b++;
            if (isalpha(tp[0]) == 0) {
                if (k>0) {
                    k=0;
                    w.clear();
                }
                m = stof(tp);
                w.push_back(m);
                j++;
            }
            if (isalpha(tp[0]) != 0)  {
                if (j>0) {
                    j=0;
                    v.push_back(w);
                }                
                k++;
            }
            if (tp.find("Time")<tp.length()){
                if (bgn==0) {
                    time_bgn=stol(tp.substr(tp.find(":")+1,tp.size()));
                    bgn++;
                }
                time_now=stol(tp.substr(tp.find(":")+1,tp.size()));
                if (time_now<time_before){
                    g++;
                    a++;
                    //cout << time_now << "      " <<time_before<<endl;
                    if (time_now<time_min) time_min=time_now;
                }
                time_before=stol(tp.substr(tp.find(":")+1,tp.size()));
                //if (g==1) cout << time_now <<endl;
            }
            //a++;
        }
        time_tot=(time_now + g*(pow(2, 31)-1) - time_bgn)*8 *pow(10,-9);
        //cout<<time_min<<endl;

        //non prendo l'ultimo evento perchè è incompleto
        myfile.close(); //close the file object.
    }    

    
    vector <float> charge (v.size());
    vector <float> amp (v.size());
    float h;
    vector <float> idx_range;

    for (int i=0; i<v.size();i++){
        h = 0;
        for (int j=0; j<100; j++){
            h += v[i][j] / 100;
        }
        for (int j=0; j<1000; j++){
            charge[i] += h-v[i][j];
        }
        if (name.find("pmt1")<name.length()){
            if (charge[i]> 72000 && charge[i]< 80000) idx_range.push_back(i);
        }      
        if (name.find("pmt2")<name.length()){
            if (charge[i]> 64000 && charge[i]< 75000) idx_range.push_back(i);
        }
    }
    return idx_range;
}

void correlation_frequency (){
    gROOT->SetBatch(kFALSE);
    vector<vector<float>> idx_ranges;
    ofstream out_file("correlation_freq.txt");

    list <string> names ={
                        "pmt1_NA_e6_100_or_run1",
                        "pmt2_NA_e6_100_or_run1"
                        };

    long int time_tot;
    for(list<string>::const_iterator name = names.begin(); name != names.end(); ++name){
        cout << "Processing: " << *name << endl; 
        idx_ranges.push_back(find_idx_range("data/" + *name + ".txt", *name, time_tot));
    }


    int max_i =max(*max_element(idx_ranges[0].begin(), idx_ranges[0].end()),
                   *max_element(idx_ranges[1].begin(), idx_ranges[1].end()));
    int coincidenze=0;
    for (int i=0; i<=max_i; i++){
        if  ((find(idx_ranges[0].begin(), idx_ranges[0].end(), i) != idx_ranges[0].end())
            && (find(idx_ranges[1].begin(), idx_ranges[1].end(), i) != idx_ranges[1].end())){
            coincidenze++;
        }
    }
    out_file << "numero coincidenze:    "<< coincidenze << "\n";
    cout << coincidenze << endl;

    double n_a= idx_ranges[0].size(), n_b= idx_ranges[1].size() ;

    cout <<n_a<<endl;
    out_file << "numero eventi nel picco pmt1:    "<<n_a << "\n";
    cout <<n_b<<endl;
    out_file << "numero eventi nel picco pmt2:    "<<n_b << "\n";



    cout <<n_a<<endl;
    cout <<n_b<<endl;

    cout << time_tot << endl;
    out_file << "tempo totale:    "<< time_tot << "\n";

    double freq_a=n_a/time_tot;
    double freq_b=n_b/time_tot;
    double time_acq= 4*pow(10,-9)*1030;

    double coinc_casuali= freq_a*freq_b*time_acq*time_tot;

    cout << freq_a <<endl;
    cout << freq_b <<endl;
    cout << coinc_casuali <<endl;
    out_file << "freq pmt1 nel picco:    "<< freq_a << "\n";
    out_file << "freq pmt2 nel picco:    "<< freq_b << "\n";
    out_file << "coincidenze casuali stimate:    "<< coinc_casuali << "\n";

    cout << "coincidenze casuali attese= " <<coinc_casuali<< " pm "<< sqrt(n_a*n_b*(n_a+n_b))*time_acq/time_tot <<endl;
    cout << "coincidenze= "<<coincidenze<<" pm "<< sqrt(coincidenze) << endl;

    out_file << "coincidenze casuali attese= " <<coinc_casuali<< " pm "<< sqrt(n_a*n_b*(n_a+n_b))*time_acq/time_tot << "\n";
    out_file << "coincidenze= "<<coincidenze<<" pm "<< sqrt(coincidenze) << "\n";


}

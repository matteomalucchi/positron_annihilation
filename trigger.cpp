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

void trigger (){
    list <string> names ={"pmt1_null", 
                        "pmt2_null"};
    vector <vector<float>> q;
    for(list<string>::const_iterator name = names.begin(); name != names.end(); ++name){
        string path = "data/" + *name + ".txt", *name
        ifstream myfile;
        vector <float> t(1030);
        vector <vector<float>> v;
        myfile.open(path, ios::in | ios::out);  
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
        }
        //non prendo l'ultimo evento perchè è incompleto
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
            charge[i] += h-v[i][j];
        }
        min =*min_element(v[i].begin(), v[i].end());
        amp[i]=h-min;
    }


}
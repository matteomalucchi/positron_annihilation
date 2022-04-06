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
#include <TSpline.h>

using namespace std;


void time_res(){
    ifstream myfile;
    string name="pmt1_NA_e6_100_run2";
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
    vector <vector<float>> v;
    vector <TH1F*> histos;
    iota(begin(t), end(t), 0);

    if (myfile.is_open()){
        string tp;
        vector <float> w;
        int k=0, m, j=0, n=0;
        while(getline(myfile, tp) && n<1000){ //read data from file object and put it into string.
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
                    n++;
                }                
                k++;
            }
        }
        //non prendo l'ultimo evento perchè è incompleto
        myfile.close(); //close the file object.
    }
    int time=25, idx=10;
    vector <float> t(time), q(time), y;
    float min =*min_element(v[idx].begin()+300, v[idx].end()-80), h=0;
    auto bin_min = find(v[i].begin()+300, v[i].end()-80, min);
    int bin_mez=0;
    float diff=1000;
    float tmez=0;
    for (int j=0; j<20; j++){
        h += v[i][j] / 20;
    }
    TCanvas *c_wave = new TCanvas(&(name + "_wave")[0], &(name + "_wave")[0]);
    TGraph* gr = new TGraph(t.size(), &t[0], &v[idx][460]);
    gr->SetNameTitle(&(name + "_wave")[0], &(name + "_wave")[0]);
    gr->SetMarkerStyle(1);
    gr->Draw("AP");

    TSpline3* sp3 = new TSpline3("sp3",gr);
    sp3->Draw("same");

    c_wave->SaveAs(&( name + "_wave.png")[0]);
    float s;
    int max =10;
    for(int i; i<max;i++){
        s=i/max;
        for(int j=0; j<time;j++){
            q[j]=s+t[j];
            y.push_back(sp3->Eval(q[j]));
        }
        for (int j=0; j<time; j++){
        // Algoritmo per ottenere time continuo
            if(abs(y[j]-min/2-h/2)<diff && j<bin_min-y.begin()){
                bin_mez=j;
                diff=abs(v[i][j]-min/2-h/2);
            }
        }       
        float TT[5],VV[5];
        for(int k=0;k<5;k++){
            TT[k] = k;
            VV[k] = v[i][bin_mez-2+k];
        }
        TGraph* gr_fit = new TGraph(5, TT,VV);
        gr_fit->Fit("pol1", "Q");
        tmez=(min/2 + h/2 - gr_fit->GetFunction("pol1")->GetParameter(0))/(gr_fit->GetFunction("pol1")->GetParameter(1));
        time[i]=(tmez+t[bin_mez-2])*4;

    }


}
    
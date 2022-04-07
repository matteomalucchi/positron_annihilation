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

ofstream time_file("triple/time_res.txt");

void time_res(){
    ifstream myfile;
    string name="pmt1_NA_c6_ext_coinc12_merc_metal_run6";
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
    vector <vector<float>> v;

    if (myfile.is_open()){
        string tp;
        vector <float> w;
        int k=0, m, j=0, n=0;
        while(getline(myfile, tp) && n<20){ //read data from file object and put it into string.
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


    int time_tot=15, idx=1, start_time=529;
    vector <float> t(time_tot), q(time_tot), y(time_tot), t_tot(1030);
    iota(begin(t), end(t), 0);
    iota(begin(t_tot), end(t_tot), 0);

    float min =*min_element(v[idx].begin()+400, v[idx].end()-500), h=0;
    auto bin_min = find(v[idx].begin()+400, v[idx].end()-500, min);

    for (int j=0; j<20; j++){
        h += v[idx][j] / 20;
    }

    TCanvas *c_tot = new TCanvas(&(name + "_tot")[0], &(name + "_tot")[0]);
    TGraph* gr_tot = new TGraph(t_tot.size(), &t_tot[0], &v[idx][0]);
    gr_tot->SetNameTitle(&(name + "_tot")[0], &(name + "_tot")[0]);
    gr_tot->SetMarkerStyle(2);
    gr_tot->Draw("AP"); 

    TCanvas *c_wave = new TCanvas(&(name + "_wave")[0], &(name + "_wave")[0]);
    TGraph* gr = new TGraph(t.size(), &t[0], &v[idx][start_time]);
    gr->SetNameTitle(&(name + "_wave")[0], &(name + "_wave")[0]);
    gr->SetMarkerStyle(2);
    gr->Draw("AP");

    TSpline3* sp3 = new TSpline3("sp3",gr);
    sp3->SetLineColor(kRed);
    sp3->Draw("same");
    c_wave->SaveAs(&( name + "_wave.png")[0]);

    float s;
    float max =100;
    int bin_mez=0;
    float diff=1000, time_mes;
    float tmez=0;
    vector<float> times(max);

    for(float i=0; i<max;i++){
        s=i/max;
        for(int j=0; j<time_tot;j++){
            q[j]=s+t[j];
            y[j]=sp3->Eval(q[j]);
            //cout <<q[j]<<"  "<<y[j]<<endl;
        }
        for (int j=0; j<time_tot; j++){
        // Algoritmo per ottenere time_tot continuo
            if(abs(y[j]-min/2-h/2)<diff /*&& j<bin_min-y.begin()*/){
                bin_mez=j;
                diff=abs(y[j]-min/2-h/2);
            }
        }
        cout <<bin_mez<<"  "<<q[bin_mez]<<"  "<< y[bin_mez]<<endl;       
        float TT[5],VV[5];
        for(int k=0;k<5;k++){
            TT[k] = k;
            VV[k] = y[bin_mez-2+k];
        }
        //TCanvas *c_wave = new TCanvas(&(to_string(i))[0], &(to_string(i))[0]);

        TGraph* gr_fit = new TGraph(5, TT,VV);
        gr_fit->Fit("pol1", "Q");
        gStyle->SetOptFit(111);

        //gr_fit->Draw();

        tmez=(min/2 + h/2 - gr_fit->GetFunction("pol1")->GetParameter(0))/(gr_fit->GetFunction("pol1")->GetParameter(1));
        time_mes=(tmez+q[bin_mez-2])*4;
        /*
        cout << "time measured  "<<s<<" =    "<<time_mes <<endl;
        cout << "time measured absolute  "<<s<<" =    "<<(start_time*4+time_mes) <<endl;
        cout <<"--------------------------"<<endl;
        cout<<endl;*/
        time_file <<  "time measured "<<s<<" =    "<<time_mes <<"\n";
        time_file << "time measured absolute  "<<s<<" =    "<<(start_time*4+time_mes) <<endl;
        time_file<<endl;

        times[i]=time_mes;
    }
    float min_times =*min_element(times.begin(), times.end());
    float max_times =*max_element(times.begin(), times.end());
    cout << "resolution=    " << max_times-min_times<<endl;
    time_file << "resolution=    " << max_times-min_times<<endl;

}
    
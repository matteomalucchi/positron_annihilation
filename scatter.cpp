#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <float.h>
#include <tuple>

#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>

#include "tools.h"

using namespace std;

tuple<vector<float>,vector<float>, vector<long unsigned int>>  charge_amp(string name){
    cout << "Processing: " << name << endl; 
    ifstream myfile;     
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
    vector <vector<float>> v;
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
    float h, min, s, u;
    vector<long unsigned int> idx_strange;
    int a=0;
    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        s=0;
        u=0;
        for (int j=0; j<100; j++){
            h += v[i][j] / 100;
            if (v[i][j]<14600 && u==0){
                //cout << "valore piccolo: "<< i <<endl;
                idx_strange.push_back(i);
                u++;
                //break;
            }
        }
        for (int j=0; j<1000; j++){
                s+= h-v[i][j];
        }                
        min =*min_element(v[i].begin(), v[i].end());

        charge.push_back(s);
        amp.push_back(h-min); 

        if (a==0){
            vector <float> t(1030);
            iota(begin(t), end(t), 0);
            TCanvas *c_wave = new TCanvas(&(name + "_wave")[0], &(name + "_wave")[0]);
            TGraph* gr = new TGraph(t.size(), &t[0], &v[idx_strange[idx_strange.size()-1]][0]);
            gr->SetNameTitle(&(name + "_wave")[0], &(name + "_wave")[0]);
            gr->SetMarkerStyle(21);
            gr->Draw("AP");
            a++;
            //c_wave->SaveAs(&("scatter/" + name + "_wave.png")[0]);            
        }    
    }
    return make_tuple(charge, amp, idx_strange);
}




void scatter(){
    map <string, string> pair_names ={
        {"pmt1_NA_e6_ext_run1","pmt2_NA_e6_ext_run1"},
        {"pmt1_NA_e6_ext_run2","pmt2_NA_e6_ext_run2"},
        {"pmt1_NA_e6_ext_run3", "pmt3_NA_e6_ext_run3"},
        {"pmt1_NA_e6_100_or_run1", "pmt2_NA_e6_100_or_run1"},
        {"pmt1_NA_e6_700_or_run2", "pmt2_NA_e6_700_or_run2"},
        {"pmt1_NA_e6_700_or_run3", "pmt2_NA_e6_700_or_run3"}
    };

    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first;
        const auto name_b = pair_name.second;

        string name = name_a.substr(5, name_b.size()-1);
        string pmt_a = name_a.substr(0, 4);
        string pmt_b = name_b.substr(0, 4);

        vector <float> charge_a, amp_a, charge_b, amp_b;
        vector<long unsigned int> idx_strange_a, idx_strange_b;

        tie(charge_a, amp_a, idx_strange_a) = charge_amp(name_a);
        tie(charge_b, amp_b, idx_strange_b) = charge_amp(name_b); 
        
        cout <<"tot = " <<charge_a.size() <<endl;
        int q=0;
        for (long unsigned int i=0; i< charge_a.size(); i++){
            if ((find(idx_strange_a.begin(), idx_strange_a.end(), i) != idx_strange_a.end()) || (find(idx_strange_b.begin(), idx_strange_b.end(), i) != idx_strange_b.end())){
            //if(charge_a[i]<-10000 || charge_b[i]<-10000){
                cout << charge_a[i] << "   " << charge_b[i] <<endl;
                charge_a.erase(charge_a.begin()+i);
                charge_b.erase(charge_b.begin()+i);    
                q++;            
            }
        }
        cout << "erased = " << q <<endl;
        cout << charge_a.size() << "   " << charge_b.size() <<endl;
        

        TCanvas *c_scatter_charge = new TCanvas(&("scatter_" + name + "_charge")[0], &("scatter_" + name + "_charge")[0]);
        TGraph* gr_charge = new TGraph(min(charge_a.size(), charge_b.size()),&charge_a[0],&charge_b[0]);
        gr_charge->SetMarkerStyle(1);
        gr_charge->SetNameTitle(&("scatter_" + name + "_charge")[0],&("scatter_" + name + "_charge")[0]);
        gr_charge->GetXaxis()->SetTitle(&(pmt_a)[0]);
        gr_charge->GetYaxis()->SetTitle(&(pmt_b)[0]);
        gr_charge->Draw("AP");
        c_scatter_charge->SaveAs(&("scatter/scatter_" + name +"_charge"+ ".png")[0]);

        TCanvas *c_scatter_amp = new TCanvas(&("scatter_" + name + "_amp")[0],&("scatter_" + name + "_amp")[0]);
        TGraph* gr_amp = new TGraph(min(amp_a.size(), amp_b.size()),&amp_a[0],&amp_b[0]);
        gr_amp->SetMarkerStyle(1);
        gr_amp->SetNameTitle(&("scatter_" + name + "_amp")[0],&("scatter_" + name + "_amp")[0]);
        gr_amp->GetXaxis()->SetTitle(&(pmt_a)[0]);
        gr_amp->GetYaxis()->SetTitle(&(pmt_b)[0]);
        gr_amp->Draw("AP");
        c_scatter_amp->SaveAs(&("scatter/scatter_" + name +"_amp" +".png")[0]);
    }
}





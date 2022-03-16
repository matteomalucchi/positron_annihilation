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
//#include<bits/stdc++.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>

using namespace std;

vector<vector<float>>  charge_amp(string name){
    cout << "Processing: " << name << endl; 
    ifstream myfile;     
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
    vector <vector<float>> v;
    if (myfile.is_open()){
        string tp;
        vector <float> w;
        int k=0, m, j=0, f=0;
        while(getline(myfile, tp)/* && f<100000000*/){ //read data from file object and put it into string.
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
        f++;
        }
        //non prendo l'ultimo evento perchè è incompleto
        myfile.close(); //close the file object.
    }    
    
    vector <float> charge(v.size());
    vector <float> amp(v.size());
    float h, min, s, u;
    vector<float> idx_strange;
    int a=0;
    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        s=0;
        u=0;
        for (int j=0; j<200; j++){
            h += v[i][j] / 200;
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

        charge[i]=(s);
        amp[i]=(h-min); 

        if (a==0 /*&& (find(idx_strange.begin(), idx_strange.end(), i) == idx_strange.end()) && charge[i]<-300000*/){
            cout << "wave   " << i << "     charge" << charge[i] <<"    amp" << amp[i]<<endl;
            vector <float> t(1030);
            iota(begin(t), end(t), 0);
            TCanvas *c_wave = new TCanvas(&(name + "_wave")[0], &(name + "_wave")[0]);
            TGraph* gr = new TGraph(t.size(), &t[0], &v[i][0]);
            gr->SetNameTitle(&(name + "_wave")[0], &(name + "_wave")[0]);
            gr->SetMarkerStyle(2);
            gr->Draw("AP");
            a++;
            //c_wave->SaveAs(&("scatter/" + name + "_wave.png")[0]);            
        }    
    }
    cout<<idx_strange.size() <<endl;
    vector<vector<float>> vecs;
    vecs.push_back(charge);
    vecs.push_back(amp);
    vecs.push_back(idx_strange);
    return vecs;
}




void scatter(){
    map <string, string> pair_names ={
        //{"pmt1_NA_e6_ext_run1","pmt2_NA_e6_ext_run1"},
        {"pmt1_NA_e6_ext_run2","pmt2_NA_e6_ext_run2"}
        //{"pmt1_NA_e6_ext_run3", "pmt3_NA_e6_ext_run3"},
        /*{"pmt1_NA_e6_100_or_run1", "pmt2_NA_e6_100_or_run1"},
        {"pmt1_NA_e6_700_or_run2", "pmt2_NA_e6_700_or_run2"},
        {"pmt1_NA_e6_700_or_run3", "pmt2_NA_e6_700_or_run3"}*/
    };

    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first;
        const auto name_b = pair_name.second;

        string name = name_a.substr(5, name_b.size()-1);
        string pmt_a = name_a.substr(0, 4);
        string pmt_b = name_b.substr(0, 4);

        vector <float> charge_a, amp_a, charge_b, amp_b, idx_strange_a, idx_strange_b;
        //vector<long unsigned int> idx_strange_a, idx_strange_b;
        vector<vector<float>> vecs_a, vecs_b;

        vecs_a = charge_amp(name_a);
        vecs_b= charge_amp(name_b);

        charge_a=vecs_a[0];
        amp_a=vecs_a[1];
        idx_strange_a=vecs_a[2];
        charge_b=vecs_b[0];
        amp_b=vecs_b[1];
        idx_strange_b=vecs_b[2];


        cout<<"strange a " <<idx_strange_a.size() <<endl;
        for (long unsigned int i=0; i< idx_strange_a.size(); i++){
           // cout<<i << "    " <<idx_strange_a[i] <<endl;
        }
        cout<<"strange b "<<idx_strange_b.size() <<endl;
        for (long unsigned int i=0; i< idx_strange_b.size(); i++){
            //cout<<i << "    " <<idx_strange_b[i] <<endl;
        }
        cout <<charge_a.size()<<endl;
        cout <<charge_b.size()<<endl;
        cout <<amp_a.size()<<endl;
        cout <<amp_b.size()<<endl;
        //for (int e=charge_a.size()-1; e< charge_a.size(); e++){
            //if (charge_a[e]==0){
            int e = 4489;
            cout << "idx  "<<e<<endl;
            cout <<charge_a[e]<<endl;
            cout <<charge_b[e]<<endl;
            cout <<amp_a[e]<<endl;
            cout <<amp_b[e]<<endl;
            //break;
        

        cout <<"tot charge= " <<charge_a.size()<< "   " << charge_b.size() <<endl;
        int q=0;
        vector<float> charge_a_def, charge_b_def, amp_a_def, amp_b_def; 

        for (int i=0; i< charge_a.size(); i++){
            if ((find(idx_strange_a.begin(), idx_strange_a.end(), i) == idx_strange_a.end()) && (find(idx_strange_b.begin(), idx_strange_b.end(), i) == idx_strange_b.end())){
            //if ((count(idx_strange_a.begin(), idx_strange_a.end(), i) ) || (count(idx_strange_b.begin(), idx_strange_b.end(), i) )){
            //if(charge_a[i]<-10000 || charge_b[i]<-10000){
                if (i==4489){
                    cout <<"---"<<endl;
                    cout<< i <<endl;
                    cout << charge_a[i] << "   " << charge_b[i] <<endl;
                    cout << charge_a[i+1] << "   " << charge_b[i+1] <<endl;
                    cout << charge_a[i-1] << "   " << charge_b[i-1] <<endl;  
                    cout <<"---"<<endl;
                  
                }

                charge_a_def.push_back(charge_a[i]);
                charge_b_def.push_back(charge_b[i]);
                amp_a_def.push_back(amp_a[i]);
                amp_b_def.push_back(amp_b[i]);

                //charge_a.erase(charge_a.begin()+i);
                //charge_b.erase(charge_b.begin()+i);  

                q++;            
            }
        }
    
        cout << "erased = " << q <<endl;
        cout << "remaining charge =" << charge_a.size() << "   " << charge_b.size() <<endl;
        cout << charge_a_def[4489]<< endl;


        // Covarianza Campione

        vector <float> picco_a,picco_b;
        for (long unsigned int n=0 ; n<charge_a_def.size() ; n++){

            if(charge_a_def[n]<80000 && charge_a_def[n]>73000 && charge_b_def[n]<74000 && charge_b_def[n]>65000){ // Intervalli tarati su intervalli di picco A dei pmt1 e pmt2 in ext_run2

                picco_a.push_back(charge_a_def[n]);
                picco_b.push_back(charge_b_def[n]);
            } 
        }
        float picco_a_med=accumulate(picco_a.begin(), picco_a.end(), 0)/picco_a.size();
        float picco_b_med=accumulate(picco_b.begin(), picco_b.end(), 0)/picco_b.size();
        float cov_camp=0;

        for(int n=0; n<picco_a.size() ; n++){
            cov_camp+=(picco_a[n]-picco_a_med)*(picco_b[n]-picco_b_med)/picco_a.size();
        }

        cout<<cov_camp<<endl;


        TCanvas *c_scatter_charge = new TCanvas(&("scatter_" + name + "_charge")[0], &("scatter_" + name + "_charge")[0]);
        TGraph* gr_charge = new TGraph(min(charge_a_def.size(), charge_b_def.size()),&charge_a_def[0],&charge_b_def[0]);
        gr_charge->SetMarkerStyle(1);
        gr_charge->SetNameTitle(&("scatter_" + name + "_charge")[0],&("scatter_" + name + "_charge")[0]);
        gr_charge->GetXaxis()->SetTitle(&(pmt_a)[0]);
        gr_charge->GetYaxis()->SetTitle(&(pmt_b)[0]);
        gr_charge->Draw("AP");
        c_scatter_charge->SaveAs(&("scatter/scatter_" + name +"_charge"+ ".png")[0]);

        TCanvas *c_scatter_amp = new TCanvas(&("scatter_" + name + "_amp")[0],&("scatter_" + name + "_amp")[0]);
        TGraph* gr_amp = new TGraph(min(amp_a_def.size(), amp_b_def.size()),&amp_a_def[0],&amp_b_def[0]);
        gr_amp->SetMarkerStyle(1);
        gr_amp->SetNameTitle(&("scatter_" + name + "_amp")[0],&("scatter_" + name + "_amp")[0]);
        gr_amp->GetXaxis()->SetTitle(&(pmt_a)[0]);
        gr_amp->GetYaxis()->SetTitle(&(pmt_b)[0]);
        gr_amp->Draw("AP");
        c_scatter_amp->SaveAs(&("scatter/scatter_" + name +"_amp" +".png")[0]);
    }
}





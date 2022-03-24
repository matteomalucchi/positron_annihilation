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


using namespace std;

vector<vector<float>> energy_time(string name){
    ifstream myfile;
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
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
    vector <float> time (v.size());
    float h, u,r, min;
    vector <float> idx_strange;
    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        u=0;
        r=0;
        for (int j=0; j<20; j++){
            h += v[i][j] / 20;
            if (v[i][j]<14600 && u==0){
                idx_strange.push_back(i);
                u++;
            }
        }
        for (int j=300; j<950; j++){
            charge[i] += h-v[i][j];
            if(v[i][j]<14600 && r==0){
                time[i]=j;
                r++;
            }
        }
        min =*min_element(v[i].begin()+300, v[i].end()-80);
        amp[i]=h-min;
    }
    vector<vector<float>> vecs;
    vecs.push_back(charge);
    vecs.push_back(amp);
    vecs.push_back(time);
    vecs.push_back(idx_strange);

    return vecs;

}
/*
vector <TH1F*> make_histo(string name, vector<float> charge, vector<float> amp, vector<float> idx_strange){
    int range_charge=250000;
    int range_amp=3500; 
    if (name.find("pmt3")<name.length()){
        range_charge=50000;
        range_amp=600;
    }
    TH1F *histo_charge=  new TH1F(&(name + "_charge")[0],&(name + "_charge")[0], 1300, 0, range_charge);
    TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 1400, 0, range_amp);
    for (long unsigned int i=0; i< charge.size(); i++){
        if (find(idx_strange.begin(), idx_strange.end(), i) == idx_strange.end()){
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

    cout << "numero di idx_strange  "<<idx_strange.size() <<endl;


    return histos;
}
*/
void triple_coincidence (){
    gROOT->SetBatch(kFALSE);
    vector <vector<TH1F*>> histos;
    vector <TH1F*> histo_pmt3_12;
    vector<vector <vector<float>>> infos;

    TFile *outfile= new TFile("histograms/histograms_triple_coincidence.root", "RECREATE"/* "UPDATE"*/);
    vector<vector<string>> names ={
        //{"pmt1_NA_e6_ext_triple_90deg_run1", "pmt2_NA_e6_ext_triple_90deg_run1", "pmt3_NA_e6_ext_triple_90deg_run1"},
        {"pmt1_NA_l1_ext_triple_close_run2", "pmt2_NA_l1_ext_triple_close_run2", "pmt3_NA_l1_ext_triple_close_run2"},
    };

    TStopwatch time_tot;
    time_tot.Start();                
    for(int i=0;i<names.size();i++){
        TFile *tree_file= new TFile("triple/prova.root", "RECREATE"/* "UPDATE"*/);

        TStopwatch time;
        time.Start();
        for (int j=0; j<names[i].size(); j++){

            infos.push_back(energy_time(names[i][j]));
            //make_tree(infos, tree_file);


            time.Stop();
            time.Print();
        }

        for (int i=0; i< infos[0][0].size(); i++){
            if ((find(infos[0][3].begin(), infos[0][3].end(), i) != infos[0][3].end())
                || (find(infos[1][3].begin(), infos[1][3].end(), i) != infos[1][3].end())
                || (find(infos[2][3].begin(), infos[2][3].end(), i) != infos[2][3].end())){
                    for (int j=0; j<3; j++){
                        for (int k=0; k<3; k++){
                            infos[j][k].erase(infos[j][k].begin()+i);
                        }
                    }
            }
        }

        for (int n=0 ; n<infos[0][0].size(); n++){
            if(infos[0][0][n]<83000 && infos[0][0][n]>75000 && infos[0][0][n]<83000 && infos[0][0][n]>75000){ // Intervalli tarati su intervalli di picco A dei pmt1 e pmt2 in ext_run2
                picco_a.push_back(charge_a_def[n]);
                picco_b.push_back(charge_b_def[n]);
                //if (a==0) {
                    //cout << picco_a[picco_a.size()-1] <<"  " << picco_b[picco_b.size()-1] << endl;
                 //   a++;}

            } 
        }


        //make_histo(names[i][j], infos.back()[0],infos.back()[1], infos.back()[2]);



    }
    //getchar();
    outfile->Close();
    time_tot.Stop();
    time_tot.Print();
}


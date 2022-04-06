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

string type_of_file = "triple_lin_run8";

ofstream time_file("triple/time_tot.txt");

vector<vector<float>> take_params(string path){
    ifstream myfile;
    list<string> pmts= {"pmt1", "pmt2", "pmt3"};
    int m=0, idx_i, idx_f;
    string number, tp;
    vector<vector<float>> params(pmts.size());
    for(list<string>::const_iterator pmt = pmts.begin(); pmt != pmts.end(); ++pmt){
        myfile.open(path, ios::in | ios::out);  
        if (myfile.is_open()){
            while(getline(myfile, tp)){ 
                if(tp.find(*pmt)<tp.length()){
                    idx_i=tp.find("=");
                    idx_f=tp.find("+");
                    number=tp.substr(idx_i+2, idx_f-idx_i-3);
                    params[m].push_back(stof(number));
                }
            }
        myfile.close(); //close the file object.
        }
    m++;
    }
    return params;
}

vector<vector<vector<float>>> energy_time(string name,long int &time_tot, vector<float> param){
    ifstream myfile;
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
    vector <float> t(1030);
    iota(begin(t), end(t), 0);
    vector <vector<float>> v;
    vector <TH1F*> histos;
    int bgn=0, g=0;
    long int time_before=0, time_now=0, time_bgn=0;

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
            if (tp.find("Time")<tp.length()){
                if (bgn==0) {
                    time_bgn=stol(tp.substr(20,tp.size()-1));
                    bgn++;
                }
                time_now=stol(tp.substr(20,tp.size()-1));
                if (time_now<time_before){
                    g++;
                }
                time_before=stol(tp.substr(20,tp.size()-1));
            }
        }
        time_tot=(time_now + g*(pow(2, 32)-1) - time_bgn)*8 *pow(10,-9);

        //non prendo l'ultimo evento perchè è incompleto
        myfile.close(); //close the file object.
    }    
    
    vector <float> charge (v.size());
    vector <float> amp (v.size());
    vector <float> time (v.size());
    float h, r,u, min;
    vector <float> idx_strange;
    vector <float> mask_strange(v.size(), 1.0);

    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        u=0;
        r=0;
        for (int j=0; j<20; j++){
            h += v[i][j] / 20;
            if (v[i][j]<14600 && u==0){
                idx_strange.push_back(i);
                mask_strange[i]=0;
                u++;
            }
        } 
        min =*min_element(v[i].begin()+300, v[i].end()-80);
        auto bin_min = find(v[i].begin()+300, v[i].end()-80, min);

        if (type_of_file.find("quad")<type_of_file.length()){
            amp[i]=(-param[4]+sqrt(param[4]*param[4]-4*param[5]*(param[3]-(h-min))))/(2*param[5]);
        }
        else if (type_of_file.find("lin")<type_of_file.length()){
            amp[i]=(h-min-param[2])/param[3];
        }
        int bin_mez=0;
        float diff=1000;
        float tmez=0;

        for (int j=300; j<950; j++){
            charge[i] += h-v[i][j];
        // Algoritmo per ottenere time continuo
            if(abs(v[i][j]-min/2-h/2)<diff && j<bin_min-v[i].begin()){
                bin_mez=j;
                diff=abs(v[i][j]-min/2-h/2);
                }
        }
        if (type_of_file.find("quad")<type_of_file.length()){
            charge[i]=(-param[1]+sqrt(param[1]*param[1]-4*param[2]*(param[0]-(charge[i]))))/(2*param[2]);
        }
        else if (type_of_file.find("lin")<type_of_file.length()){ 
            charge[i]=(charge[i]-param[0])/param[1];  
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

    vector<vector<float>> vecs;
    vecs.push_back(charge);
    vecs.push_back(amp);
    vecs.push_back(time);
    vecs.push_back(mask_strange);

    vector<vector<vector<float>>> vecs_final;
    vecs_final.push_back(vecs);
    vecs_final.push_back(v);

    return vecs_final;
}

void triple_coincidence (){
    gROOT->SetBatch(kTRUE);
    ROOT::EnableImplicitMT();
    vector <vector<TH1F*>> histos;
    vector<vector <vector<float>>> infos, waves, tot_vec;
    vector <float> t(1030);
    iota(begin(t), end(t), 0);
    long int time_tot, y;

    TFile *tree_file= new TFile(&("triple/ntuple_"+ type_of_file +".root")[0], "RECREATE"/* "UPDATE"*/);

    //TFile *outfile= new TFile("triple/waves.root", "RECREATE"/* "UPDATE"*/);
    vector<vector<string>> names ={
        /*{"pmt1_NA_e6_ext_triple_90deg_run1", "pmt2_NA_e6_ext_triple_90deg_run1", "pmt3_NA_e6_ext_triple_90deg_run1"},
        {"pmt1_NA_l1_ext_triple_close_run2", "pmt2_NA_l1_ext_triple_close_run2", "pmt3_NA_l1_ext_triple_close_run2"},
        {"pmt1_NA_l1_ext_triple_close_run3", "pmt2_NA_l1_ext_triple_close_run3", "pmt3_NA_l1_ext_triple_close_run3"},
        {"pmt1_NA_c6_ext_triple_merc_aero_run4", "pmt2_NA_c6_ext_triple_merc_aero_run4", "pmt3_NA_c6_ext_triple_merc_aero_run4"},
        {"pmt1_NA_c6_ext_coinc12_merc_metal_run5", "pmt2_NA_c6_ext_coinc12_merc_metal_run5", "pmt3_NA_c6_ext_coinc12_merc_metal_run5"},
        {"pmt1_NA_c6_ext_coinc12_merc_metal_run6", "pmt2_NA_c6_ext_coinc12_merc_metal_run6", "pmt3_NA_c6_ext_coinc12_merc_metal_run6"},
        {"pmt1_NA_c6_ext_coinc12_merc_metal_block_run7", "pmt2_NA_c6_ext_coinc12_merc_metal_block_run7", "pmt3_NA_c6_ext_coinc12_merc_metal_block_run7"},*/
        {"pmt1_NA_c6_ext_coinc12_merc_metal_block_run8", "pmt2_NA_c6_ext_coinc12_merc_metal_block_run8", "pmt3_NA_c6_ext_coinc12_merc_metal_block_run8"},
    };
    vector<vector<float>> params;
    vector<TNtuple*> ntuples;
    if (type_of_file.find("real")<type_of_file.length()){
        if (type_of_file.find("quad")<type_of_file.length()){
            if (type_of_file.find("fix")<type_of_file.length()){
                params= take_params("real_time_calibration/lin_params_low_quadfit_fix.txt"); 
            }          
            else params= take_params("real_time_calibration/lin_params_low_quadfit.txt");        
        } 
        else if (type_of_file.find("lin")<type_of_file.length()){
            params= take_params("real_time_calibration/lin_params_low_new_ranges.txt");
        }
    }

    else if (type_of_file.find("triple_lin")<type_of_file.length()){
        params= take_params("triple_calib/lin_params.txt");
    }

    // loop over various runs
    for(int i=0;i<names.size();i++){
        string run = names[i][0].substr(names[i][0].size()-4, names[i][0].size()-1);
        // loop over various pmt
        for (int j=0; j<names[i].size(); j++){
            string pmt_name = names[i][j].substr(0,4);
            tot_vec=energy_time(names[i][j], time_tot, params[j]);
            infos.push_back(tot_vec[0]);
            waves.push_back(tot_vec[1]);
        }
        cout << "tot time    "<< time_tot <<" s"<<endl;
        time_file<<run <<"  tot time    "<< time_tot <<" s"<< "\n";


        //TTree* tree=new TTree(&("waveforms_"+run)[0],&("waveforms_"+run)[0]);

        y=0;
        // save waveforms
        if (run=="run8"){
            for  (int o=0; o<infos[3*i][0].size(); o++){
                if ( infos[3*i+0][0][o]<0.45 && infos[3*i+1][0][o]<0.45 &&infos[3*i+2][0][o]<0.45 && 
                    abs(infos[3*i+0][2][o]-infos[3*i+1][2][o])<13 &&(abs(infos[3*i+2][2][o]-infos[3*i+1][2][o])<13  || abs(infos[3*i+0][2][o]-infos[3*i+2][2][o])<13)
                    && infos[3*i+0][3][o]*infos[3*i+1][3][o]*infos[3*i+2][3][o]==1
                    /*&&infos[3*i+2][0][o]+infos[3*i+0][0][o]+infos[3*i+1][0][o]>1.1*/ ){
                    TCanvas *c_wave = new TCanvas(&(run +"_wave_"+o)[0], &(run +"_wave_"+o)[0]);
                    //tree->Branch(&(run +"_wave_"+o)[0], &c_wave);

                    TGraph* gr1 = new TGraph(t.size(), &t[0], &waves[i*3][o][0]);
                    gr1->GetYaxis()->SetRangeUser(13000, 14800);
                    gr1->SetNameTitle(&("pmt1_"+run+"_wave_"+infos[3*i+0][0][o])[0], &("pmt1_"+run+ "_wave_"+infos[3*i+0][0][o])[0]);
                    gr1->SetLineColor(kGreen);
                    gr1->Draw();

                    TGraph* gr2 = new TGraph(t.size(), &t[0], &waves[i*3+1][o][0]);
                    gr2->SetNameTitle(&("pmt2_"+run+"_wave_"+infos[3*i+1][0][o])[0], &("pmt2_"+run+ "_wave_"+infos[3*i+1][0][o])[0]);
                    gr2->SetLineColor(kBlue);
                    gr2->Draw("same");

                    TGraph* gr3 = new TGraph(t.size(), &t[0], &waves[i*3+2][o][0]);
                    gr3->SetNameTitle(&("pmt3_"+run+"_wave_"+infos[3*i+2][0][o])[0], &("pmt3_"+run+ "_wave_"+infos[3*i+2][0][o])[0]);
                    gr3->SetLineColor(kRed);
                    gr3->Draw("same");

                    c_wave->BuildLegend();
                    c_wave->SaveAs(&("triple/waves_8/" +run+ "_wave_"+to_string(infos[3*i+2][0][o]+infos[3*i+0][0][o]+infos[3*i+1][0][o])+".png")[0]);
                    c_wave->Write();
                    //tree->Fill();
                    y++;
                }
            }
            //tree->Write();
        }
        
    }

/*
    // eliminate strange events
    int u;
    for (int k=0;k<names.size();k++){
    //cout <<infos[k][0].size()<<endl;
        u=0;
        for (int p=0; p< infos[k*3][0].size(); p++){
            if (infos[k*3][3][p]==0 || infos[k*3+1][3][p]==0 || infos[k*3+2][3][p]==0){
                for (int t=0; t<names[k].size();t++){
                    for(int j=0; j<infos[k*3].size(); j++){
                        //cout << infos[t+k][j].size()<<endl;
                        infos[t+k*3][j].erase(infos[t+k*3][j].begin()+p);
                        //cout << infos[t+k][j].size()<<endl; 
                    }    
                }      
                p--; 
            }
            //if (p>=infos[k*3][0].size()-1) break;
            //cout <<p <<endl;
            u++;
        }
    }
    //cout <<u <<endl;
*/

    //loop over run
    for(int j=0; j<names.size();j++){
        string run = names[j][0].substr(names[j][0].size()-4, names[j][0].size()-1);
        // each branch is charge, amp, time and mask_strange_tot
        TNtuple *ntuple= new TNtuple(&run[0], &run[0], "c1:a1:t1:c2:a2:t2:c3:a3:t3:m");
        for (int m=0; m < infos[3*j][0].size(); m++){
            ntuple->Fill(infos[3*j+0][0][m],infos[3*j+0][1][m],infos[3*j+0][2][m],
                        infos[3*j+1][0][m],infos[3*j+1][1][m],infos[3*j+1][2][m],
                        infos[3*j+2][0][m],infos[3*j+2][1][m],infos[3*j+2][2][m],
                        infos[3*j+0][3][m]*infos[3*j+1][3][m]*infos[3*j+2][3][m]);
        }
        ntuples.push_back(ntuple);
    }

    tree_file->cd();
    for(int j=0; j<ntuples.size();j++){
        ntuples[j]->Write();
    }


    tree_file->Close();
}


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
        while(getline(myfile, tp)/* && f<5000000*/){ //read data from file object and put it into string.
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

    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        s=0;
        u=0;
        for (int j=0; j<200; j++){
            h += v[i][j] / 200;
            if (v[i][j]<14600 && u==0){
                idx_strange.push_back(i);
                u++;
            }
        }
        for (int j=0; j<1000; j++){
                s+= h-v[i][j];
        }                
        min =*min_element(v[i].begin(), v[i].end());

        charge[i]=(s);
        amp[i]=(h-min);   
    }
    vector<vector<float>> vecs;
    vecs.push_back(charge);
    vecs.push_back(amp);
    vecs.push_back(idx_strange);
    return vecs;
}




void scatter(){
    map <string, string> pair_names ={
        {"pmt1_NA_e6_ext_run1","pmt2_NA_e6_ext_run1"},
        //{"pmt1_NA_e6_ext_run2","pmt2_NA_e6_ext_run2"},
        //{"pmt1_NA_e6_ext_run3", "pmt3_NA_e6_ext_run3"},
        /*{"pmt1_NA_e6_100_or_run1", "pmt2_NA_e6_100_or_run1"},
        {"pmt1_NA_e6_700_or_run2", "pmt2_NA_e6_700_or_run2"},
        {"pmt1_NA_e6_700_or_run3", "pmt2_NA_e6_700_or_run3"}*/
    };
    TFile *outfile= new TFile("scatter/scatter_plots_run1_ext.root", "RECREATE");

    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first;
        const auto name_b = pair_name.second;

        string name = name_a.substr(5, name_b.size()-1);
        string pmt_a = name_a.substr(0, 4);
        string pmt_b = name_b.substr(0, 4);

        vector <float> charge_a, amp_a, charge_b, amp_b, idx_strange_a, idx_strange_b;
        vector<vector<float>> vecs_a, vecs_b;

        vecs_a = charge_amp(name_a);
        vecs_b= charge_amp(name_b);

        charge_a=vecs_a[0];
        amp_a=vecs_a[1];
        idx_strange_a=vecs_a[2];
        charge_b=vecs_b[0];
        amp_b=vecs_b[1];
        idx_strange_b=vecs_b[2];

        cout <<"tot charge= " <<charge_a.size()<< "   " << charge_b.size() <<endl;
        int q=0;
        vector<float> charge_a_def, charge_b_def, amp_a_def, amp_b_def; 

        for (int i=0; i< charge_a.size(); i++){
            if ((find(idx_strange_a.begin(), idx_strange_a.end(), i) == idx_strange_a.end())
                && (find(idx_strange_b.begin(), idx_strange_b.end(), i) == idx_strange_b.end())){
                charge_a_def.push_back(charge_a[i]);
                charge_b_def.push_back(charge_b[i]);
                amp_a_def.push_back(amp_a[i]);
                amp_b_def.push_back(amp_b[i]);

            }
        }


        // Covarianza Campione

        vector <float> picco_a,picco_b;
        int a=0;
        for (long unsigned int n=0 ; n<charge_a_def.size() ; n++){

            if(charge_a_def[n]<80000 && charge_a_def[n]>73000 && charge_b_def[n]<74000 && charge_b_def[n]>65000){ // Intervalli tarati su intervalli di picco A dei pmt1 e pmt2 in ext_run2
                picco_a.push_back(charge_a_def[n]);
                picco_b.push_back(charge_b_def[n]);
                //if (a==0) {
                    //cout << picco_a[picco_a.size()-1] <<"  " << picco_b[picco_b.size()-1] << endl;
                 //   a++;}

            } 
        }
   /*     float picco_a_med=accumulate(picco_a.begin(), picco_a.end(), 0)/picco_a.size();
        float picco_b_med=accumulate(picco_b.begin(), picco_b.end(), 0)/picco_b.size();

        cout<<"Media pmt1= "<<picco_a_med<<endl;
        cout<<"Media pmt2= "<<picco_b_med<<endl;


        float cov_camp=0;
        int negatives=0;
        for(int n=0; n<picco_a.size() ; n++){
            cov_camp+=(picco_a[n]-picco_a_med)*(picco_b[n]-picco_b_med)/picco_a.size();
            //cout <<(picco_a[n]-picco_a_med)*(picco_b[n]-picco_b_med)/picco_a.size()<<endl;
            if ((picco_a[n]-picco_a_med)*(picco_b[n]-picco_b_med) < 0){
                negatives++;
            } 
        }

        cout<<"Dimensione campione= "<<picco_a.size()<<endl;
        cout<<"Negatives = "<<negatives<<endl;
        cout<<"Covarianza campione= "<<cov_camp<<endl;
*/



        // Fit di scatter per trovare coefficiente di correlazione quando si trovano vicino al picco    

        TCanvas *c_fit_lin = new TCanvas(&("scatter_" + name + "_charge_correlationfit")[0], &("scatter_" + name + "_charge_correlationfit")[0]);
        TGraphErrors* gr = new TGraphErrors(picco_a.size(),&picco_a[0],&picco_b[0],nullptr,nullptr);
/*        c_fit_lin->Update();
        TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
        st->SetX1NDC(0.7); //new x start position
        st->SetX2NDC(0.9); //new x end position
        st->SetY1NDC(0.7); //new x start position
        st->SetY2NDC(0.9); //new x end position
        gr->SetMarkerStyle(1);*/
        gr->SetNameTitle(&("scatter plot (pmt1 vs pmt2) charge")[0],&("scatter plot (pmt1 vs pmt2) charge")[0]);
        gr->GetXaxis()->SetTitle(&(pmt_a+" charge [a.u.]")[0]);
        gr->GetYaxis()->SetTitle(&(pmt_b+" charge [a.u.]")[0]);
        gr->Draw("APE");
        TF1 *   linear  = new TF1("linear","[0]+[1]*x",-1, 1.4,"W");
        linear->SetParNames ("Off-set","Correlation factor");
        gr->Fit("linear");
        linear->Draw("SAME");
        gStyle->SetOptFit(111);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        c_fit_lin->Write();

        c_fit_lin->SaveAs(&("scatter/scatter_linear" + name +"_charge"+ ".pdf")[0]);



        TCanvas *c_scatter_charge = new TCanvas(&("scatter_" + name + "_charge")[0], &("scatter_" + name + "_charge")[0]);
        TGraph* gr_charge = new TGraph(min(charge_a_def.size(), charge_b_def.size()),&charge_a_def[0],&charge_b_def[0]);
        gr_charge->SetMarkerStyle(1);
        gr_charge->SetNameTitle(&("scatter plot (pmt1 vs pmt2) charge")[0],&("scatter plot (pmt1 vs pmt2) charge")[0]);
        gr_charge->GetXaxis()->SetTitle(&(pmt_a+" charge [a.u.]")[0]);
        gr_charge->GetYaxis()->SetTitle(&(pmt_b+" charge [a.u.]")[0]);
        gPad->SetLeftMargin(0.15);
        gr_charge->Draw("AP");
        c_scatter_charge->SaveAs(&("scatter/scatter_" + name +"_charge"+ ".pdf")[0]);
        c_scatter_charge->Write();

        TCanvas *c_scatter_amp = new TCanvas(&("scatter_" + name + "_amp")[0],&("scatter_" + name + "_amp")[0]);
        TGraph* gr_amp = new TGraph(min(amp_a_def.size(), amp_b_def.size()),&amp_a_def[0],&amp_b_def[0]);
        gr_amp->SetMarkerStyle(1);
        gr_amp->SetNameTitle(&("scatter plot (pmt1 vs pmt2) amplitude")[0], &("scatter plot (pmt1 vs pmt2) amplitude")[0]);
        gr_amp->GetXaxis()->SetTitle(&(pmt_a+" amplitude [a.u.]")[0]);
        gr_amp->GetYaxis()->SetTitle(&(pmt_b+" amplitude [a.u.]")[0]);
        gPad->SetLeftMargin(0.15);
        gr_amp->Draw("AP");
        c_scatter_amp->SaveAs(&("scatter/scatter_" + name +"_amp" +".pdf")[0]);
        c_scatter_amp->Write();
    }
    outfile->Close();
}





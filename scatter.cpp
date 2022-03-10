#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <float.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>

#include "tools.h"

using namespace std;

void scatter(){
    ifstream myfile;

    list <string> names ={"pmt1_NA_e6_ext_run1",
                        "pmt2_NA_e6_ext_run1",};
    vector <vector<float>> charge (2);
    vector <vector<float>> amp (2);
    int n=0, u=0;
    for (list<string>::const_iterator name = names.begin(); name != names.end(); ++name){                    
        myfile.open("data/" + *name + ".txt", ios::in | ios::out);  
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
        

        float h, min, s;
        long unsigned int idx_strange=-1;

        for (long unsigned int i=0; i<v.size(); i++){
            h = 0;
            s=0;
            for (int j=0; j<200; j++){
                h += v[i][j] / 200;
                if (v[i][j]<14600){
                    cout << "valore piccolo: "<< i <<endl;
                    idx_strange=i;
                    u++;
                    break;
                }
            }
            for (int j=0; j<1000; j++){
                 s+= h-v[i][j];
            }                
            min =*min_element(v[i].begin(), v[i].end());

            if (i != idx_strange){
                charge[n].push_back(s);
                amp[n].push_back(h-min);      
            }
        }
        n++;
    }
    cout << u << endl;
    TCanvas *c_scatter_charge = new TCanvas("scatter plot charge", "scatter plot charge");
    TGraph* gr_charge = new TGraph(charge[0].size(),&charge[0][0],&charge[1][0]);
    gr_charge->SetMarkerStyle(1);
    gr_charge->SetNameTitle("scatter plot charge", "scatter plot charge");
    gr_charge->Draw("AP");
    c_scatter_charge->SaveAs("images/scatter_plot_ext_charge.png");

    TCanvas *c_scatter_amp = new TCanvas("scatter plot amp", "scatter plot amp");
    TGraph* gr_amp = new TGraph(amp[0].size(),&amp[0][0],&amp[1][0]);
    gr_amp->SetMarkerStyle(1);
    gr_amp->SetNameTitle("scatter plot amp", "scatter plot amp");
    gr_amp->Draw("AP");
    c_scatter_amp->SaveAs("images/scatter_plot_ext_amp.png");
}





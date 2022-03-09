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

    list <string> names ={"coincidence_pmt1_ext",
                        "coincidence_pmt2_ext",};
    vector <vector<float>> charge (2);
    vector <vector<float>> amp (2);
    int n=0;
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

        for (long unsigned int i=0; i<v.size(); i++){
            h = 0;
            s=0;
            for (int j=0; j<200; j++){
                h += v[i][j] / 200;
            }
            for (int j=400; j<900; j++){
                 s+= h-v[i][j];
            }
            charge[n].push_back(s);
            min =*min_element(v[i].begin(), v[i].end());
            amp[n].push_back(h-min);      
        }
        n++;
    }
    TCanvas *c_scatter = new TCanvas("scatter plot", "scatter plot");
    TGraph* gr = new TGraph(charge[0].size(),&charge[0][0],&charge[1][0]);
    gr->SetMarkerStyle(1);
    gr->Draw("AP");
    c_scatter->SaveAs("scatter_plot_ext.png");
}





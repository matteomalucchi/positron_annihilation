#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
    
using namespace std;

void plot (){
    ifstream myfile;
    myfile.open("data/wave0_co_1.txt", ios::in | ios::out);  
    vector <double> t(1030);
    vector <vector<double>> v;
    iota(begin(t), end(t), 0);

    if (myfile.is_open()){
        string tp;
        vector <double> w;
        int k=0, m, j=0;
        while(getline(myfile, tp)){ //read data from file object and put it into string.
            if (isalpha(tp[0]) == 0) {
                if (k>0) {
                    k=0;
                    w.clear();
                }
                m = stod(tp);
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
        //non prendo l'ultimo evetno perchè è incompleto
        myfile.close(); //close the file object.
        //cout<< v[v.size()-1][v[v.size()-1].size()-1] <<endl;

    }    
    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example");
    TGraph* gr = new TGraph(t.size(), &t[0], &v[1][0]);
    //gr->Draw("AC*");
    vector <double> charge (v.size());
    double h;
    TCanvas *c2 = new TCanvas("c2","histo charge");

    for (int i=0; i<v.size(); i++){
        h = 0;
        for (int j=0; j<200; j++){
            h += v[i][j] / 200;
        }
        for (int j=400; j<900; j++){
            charge[i] += abs(h-v[i][j]);
        }
    }
    TH1F *isto_co_1=  new TH1F("isto_co_1","co_1", 1000, 200000, 600000);
    for (int i=0; i< charge.size(); i++){
        isto_co_1->Fill(charge[i]);
    }
    isto_co_1->Draw();

}


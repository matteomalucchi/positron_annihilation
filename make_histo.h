
auto make_histo(string path, TH1F* histo){
    ifstream myfile;
    myfile.open(path, ios::in | ios::out);  
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

    }    
    vector <double> charge (v.size());
    double h;

    for (int i=0; i<v.size(); i++){
        h = 0;
        for (int j=0; j<200; j++){
            h += v[i][j] / 200;
        }
        for (int j=400; j<900; j++){
            charge[i] += abs(h-v[i][j]);
        }
    }
    for (int i=0; i< charge.size(); i++){
        histo->Fill(charge[i]);
    }
    return histo;
}



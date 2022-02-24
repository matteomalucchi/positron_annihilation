#include <iostream>
#include <fstream>
#include <string>
#include <vector>
    
using namespace std;

int main (){
    ifstream myfile;
    myfile.open("data/wave0_co_1", ios::in | ios::out);  
    string tp;
    vector <vector<int>> v;
    int i=0, k=0;
    while(getline(myfile, tp)){ //read data from file object and put it into string.
        if (isalpha(tp[0]) == 0) {
            if (k>0) {
                k=0;
                i++;
            }
            v[i].push_back(stoi(tp));
            cout<< tp <<endl;
        }
        else k++;
    }
    myfile.close(); //close the file object.
    //cout<< v[0][0] <<endl;
    vector <int> t(1030);
    iota(begin(t), end(t), 0);

    


}


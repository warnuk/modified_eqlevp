#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdio.h>
#include <array>

using namespace std;

string lower_to_upper(string input) {
  for (string::size_type i=0; i<input.length(); ++i)
         input[i] = toupper(input[i]);
  return input;
}

string upper_to_lower(string input) {
  for (string::size_type i=0; i<input.length(); ++i)
         input[i] = tolower(input[i]);
  return input;
}

int main() {
    string label, system, units;
    double temp, dens, ph, pco2, dil, max_sal, pkmol, pkeq;
    double na, k, li, ca, mg, cl, so4, no3, b, si, alk;
    int output_step, print_step, verbose, output;
    vector<string> add_min; vector<string> rem_min;

    ifstream file;
    string line;
    
    file.open("input.dat");

    while (file.good()) {
        getline(file, line);
        string name, value;
        stringstream linestream(line);
        

        getline(linestream, name, ',');

        if (name == "label") {
            getline(linestream, value, ',');
            label = value;
        } else if (name == "temp") {
            getline(linestream, value, ',');
            temp = stod(value);
        } else if (name == "dens") {
            getline(linestream, value, ',');
            dens = stod(value);
        } else if (name == "ph") {
            getline(linestream, value, ',');
            ph = stod(value);
        } else if (name == "na") {
            getline(linestream, value, ',');
            na = stod(value);
        } else if (name == "k") {
            getline(linestream, value, ',');
            k = stod(value);
        } else if (name == "li") {
            getline(linestream, value, ',');
            li = stod(value);
        } else if (name == "ca") {
            getline(linestream, value, ',');
            ca = stod(value);
        } else if (name == "mg") {
            getline(linestream, value, ',');
            mg = stod(value);
        } else if (name == "cl") {
            getline(linestream, value, ',');
            cl = stod(value);
        } else if (name == "so4") {
            getline(linestream, value, ',');
            so4 = stod(value);
        } else if (name == "no3") {
            getline(linestream, value, ',');
            no3 = stod(value);
        } else if (name == "b") {
            getline(linestream, value, ',');
            b = stod(value);
        } else if (name == "si") {
            getline(linestream, value, ',');
            si = stod(value);
        } else if (name == "alk") {
            getline(linestream, value, ',');
            alk = stod(value);
        } else if (name == "pco2") {
            getline(linestream, value, ',');
            pco2 = stod(value);
        } else if (name == "system") {
            getline(linestream, value, ',');
            system = value;
        } else if (name == "units") {
            getline(linestream, value, ',');
            units = value;
        } else if (name == "dil") {
            getline(linestream, value, ',');
            dil = stod(value);
        } else if (name == "add") {
            while (linestream.good()) {
                string substr;
                getline(linestream, substr, ',');
                add_min.push_back(lower_to_upper(substr));
            }     
        } else if (name == "remove") {
            while (linestream.good()) {
                string substr;
                getline(linestream, substr, ',');
                rem_min.push_back(lower_to_upper(substr));
            } 
        } else if (name == "max_sal") {
            getline(linestream, value, ',');
            max_sal = stod(value);
        } else if (name == "pkmol") {
            getline(linestream, value, ',');
            pkmol = stod(value);
        } else if (name == "pkeq") {
            getline(linestream, value, ',');
                pkeq = stod(value);
        } else if (name == "print_step") {
            getline(linestream, value, ',');
            print_step = stoi(value);
        } else if (name == "output_step") {
            getline(linestream, value, ',');
            output_step = stoi(value);
        } else if (name == "verbose") {
            getline(linestream, value, ',');
            verbose = stoi(value);
        } else if (name == "output") {
            getline(linestream, value, ',');
            output = stoi(value);
        }
    }
    file.close();

    cout << "add:" << endl;
    for (int i=0; i<add_min.size(); i++) {
        cout << "  " << add_min[i] << endl;
    }

    cout << "remove:" << endl;
    for (int i=0; i<rem_min.size(); i++) {
        cout << "  " << rem_min[i] << endl;
    }

    return 0;
}
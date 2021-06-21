#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <stdio.h>

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

vector<string> comma_split(string input) {
  vector<string> output;
  stringstream stream(input);
    while(stream.good()) {
      string substr;
      getline(stream, substr, ',');
      output.push_back(substr);
    };
  return output;
}

class Simulation {
  public:

  // attributes
  int n = 25; int ntot = 12; int ncomplex = 14; int n0 = 10;
  int max_cat = 5; int max_an = 5; int max_constit = 30;
  int max_conv = 100; double epsilon = 1e-8; double eps = 1e-12;
  double pk = 0.1; double pk0 = 0.1; double pkf = 0.0001;
  double pkstd = 0.1; double dph = 0.2; double mh2o = 55.51;
  int ndepact = 0; int nconv = 2;
  double diltot = 1; double poa = 0.0;
  double tinit; int nminer; string units;
  int nc, na, nm, nt, nm0;
  double sc, sa, cmax, amax, dca, delta;
  int icat, iani;
  string min_S;
  double aw, fi;
  int simulation_loop;
  bool error;
  double a, b, c, xu, eq;
  double stdi, ee, s;
  
  string label; 
  double temp, dens, ph;
  double pco2, max_salinity;
  string system; string output_units;
  double dilute;
  string add_minerals, rem_minerals;
  vector<string> add_min; vector<string> rem_min;
  double pkmol, pkeq;
  int print_step, output_step, verbose, output;
  string log_file, event_file, chem_file, min_file;
  ifstream file; string line;
  

  // Initialize arrays of zeros
  double psc[15] = { 0 }; double tot[13] = { 0 };
  double tot0[13] = { 0 }; double totinit[13] = { 0 };
  double molal[26] = { 0 }; double act[26] = { 0 };
  double gact[26] = { 0 }; double atom[26] = { 0 };
  double xx[26] = { 0 }; double z[26][26] = { 0 };
  double zz[26] = { 0 }; double cat[6] = { 0 };
  double ani[6] = { 0 };  int nch[26] = { 0 };
  int kmat[26][26] = { 0 }; int ica[26] = { 0 };
  int nchcat[6] = { 0 }; int nchani[6] = { 0 };

  vector<string> constituants;
  vector<string> aq_S;

  // constructor
  Simulation(string input_file) {
    // Declare temporary input chemistry variables
    double na, k, li, ca, mg, cl, so4, no3, b, si, alk;

    // Open the input file into the file stream
    file.open(input_file);

    // Read the input data into simulation attributes
    file >> label; file >> temp; file >> dens; 
    file >> ph; file >> na; file >> k;
    file >> li; file >> ca; file >> mg;
    file >> cl; file >> so4; file >> no3;
    file >> b; file >> si; file >> alk;
    file >> pco2; file >> system; file >> output_units;
    file >> dilute; file >> add_minerals; 
    file >> rem_minerals; file >> max_salinity;
    file >> pkmol; file >> pkeq;
    file >> print_step; file >> output_step;
    file >> verbose; file >> output;

    // Convert input chemistry from mmoles to moles, pH to [H+]
    tot[1] = na / 1000; tot[2] = k / 1000; tot[3] = li / 1000;
    tot[4] = ca / 1000; tot[5] = mg / 1000; tot[6] = cl / 1000;
    tot[7] = so4 / 1000; tot[8] = no3 / 1000; tot[9] = b / 1000;
    tot[10] = si / 1000; tot[11] = pow(10, -ph); tot[12] = alk / 1000;

    // Close the input file stream
    file.close();
  }

  // primary simulation methods
  void run_eql() {
    
    // initialize the simulation

      // files
    log_file = label + ".log";
    event_file = label + ".j" + system + "@";
    chem_file = label + ".j" + system + "&";
    min_file = label + ".j" + system + "%";

    // Print the greeting and current date/time
    if (verbose == 1) {
        time_t now;
        time (&now);

        cout << "\nThis is EQL..............\n";
        cout << asctime(localtime(&now)) << "\n";
    };

    add_min = comma_split(add_minerals);
    rem_min = comma_split(rem_minerals);
    
    tinit = temp;
    nminer = 0;

    for (int i = 1; i <= ntot; i++) {
      if (tot[i] > 0) nminer ++;
    }

    if (dens == 0) dens = 1;
    
    if (dens == 1) {
      units = "molal";
    } else {
      units = "molar";
    }

    // Open and read "aqu.dat"
    file.open("aqu.dat");
    for (int i = 0; i <= n; i++) {
      string text; 
      vector<string> line;
      file >> text;
      line = comma_split(text);
      
      aq_S.push_back(line[0]);
      atom[i] = stod(line[1]);
      nch[i] = stod(line[2]);
    }
    file.close();

    for (int i=1; i<=n; i++) ica[i] = 1;
    
    // Open and read "matrice1"
    //ifstream matrice1;
    //matrice1.open("matrice1");
    //string line;
    file.open("matrice1");
    for (int i=1; i<=n; i++) {
      getline(file, line);
      stringstream linestream(line);
      string value;
      for (int j=1; j<=n; j++) {
        getline(linestream, value, ',');
        kmat[i][j] = stoi(value);
      }
    }
    file.close();

    kmat[13][0] = -1; kmat[15][0] = -1; kmat[20][0] = 3; kmat[21][0] = 5;

    // Open and read "complex 3"
    file.open("complex3");
    for (int i=1; i<=14; i++) {
      getline(file, line);
      stringstream linestream(line);
      string value;
      double at, bt, ct, dt, et;
      
      getline(linestream, value, ',');
      getline(linestream, value, ',');
      at = stod(value);
      getline(linestream, value, ',');
      bt = stod(value);
      getline(linestream, value, ',');
      ct = stod(value);
      getline(linestream, value, ',');
      dt = stod(value);
      getline(linestream, value, ',');
      et = stod(value);

      psc[i] = pow(10, (at +
                        bt/300*temp +
                        ct/30000*pow(temp, 2) +
                        dt/3000000*pow(temp,3) +
                        et/300000000*pow(temp,4)));
    }
    file.close();

    // Open and read murtf2
    {
      file.open("murtf2");
      getline(file, line);
      stringstream linestream(line);
      string value;
      getline(linestream, value, ',');
      nc = stoi(value);
      getline(linestream, value, ',');
      na = stoi(value);
      getline(linestream, value, ',');
      nm = stoi(value);
      nt = nc + na;
      file.close();
    }

    // Initialize arrays of zeros
    vector<string> mineral_S, mineral0_S;

    double wmin[nm+1][nt+1] = { 0 };
    double mu[nt+1] = { 0 };
    double mum[nm+1] = { 0 };
    double psol[nm+1] = { 0 };
    double pai[nm+1] = { 0 };
    int kinvar[nt+1] = { 0 };
    int lmin[nm+1] = { 0 };
    int nwmin[nm+1] = { 0 };

    {
      // Open and read murtf2
      file.open("murtf2");
      
      // Ignore the first line
      getline(file, line);

      // Iterate through the chemical species
      for (int i=0; i<=14; i++) {
        getline(file, line);
        stringstream linestream(line);
        string value, x_S;
        double at, bt, ct, dt, et;
        
        getline(linestream, value, ',');
        x_S = upper_to_lower(value);
        getline(linestream, value, ',');
        at = stod(value);
        getline(linestream, value, ',');
        bt = stod(value);
        getline(linestream, value, ',');
        ct = stod(value);
        getline(linestream, value, ',');
        dt = stod(value);
        getline(linestream, value, ',');
        et = stod(value);

        for (int j=0; j<=nt; j++) {

          if (x_S == aq_S.at(j)) {

            mu[j] = (at + 
                    bt / 300 * temp + 
                    ct / 30000 * pow(temp, 2) +
                    dt / 3000000 * pow(temp, 3) + 
                    et / 300000000 * pow(temp, 4));
          }
        }
      }
      string x = "n/a";
      mineral_S.push_back(x);
      for (int k=1; k<=nm; k++) {
        getline(file, line);
        stringstream linestream(line);
        string value, x_S, nom_ion_S;
        double at, bt, ct, dt, et;
        int ncomp;
        double c_ion;

        getline(linestream, value, ',');
        mineral_S.push_back(value);
        getline(linestream, value, ',');
        ncomp = stoi(value);

        for (int i=1; i<=ncomp; i++) {
          getline(linestream, value, ',');
          c_ion = stod(value);
          getline(linestream, value, ',');
          nom_ion_S = value;
          x_S = upper_to_lower(nom_ion_S);

          for (int j=0; j<=nt; j++) {
            if (x_S == aq_S.at(j)) {
              wmin[k][j] = c_ion;
            }
          }
        }

        getline(linestream, value, ',');
        at = stod(value);
        getline(linestream, value, ',');
        bt = stod(value);
        getline(linestream, value, ',');
        ct = stod(value);
        getline(linestream, value, ',');
        dt = stod(value);
        getline(linestream, value, ',');
        et = stod(value);

        mum[k] = (at +
                  bt / 300 * temp +
                  ct / 30000 * pow(temp, 2) +
                  dt / 3000000 * pow(temp, 3) + 
                  et / 300000000 * pow(temp, 4));
      }
      file.close();
    }

    for (int k=1; k<=nm; k++) {
      double u = mum[k];
      for (int i=0; i<=nt; i++) {
        u = u - wmin[k][i] * mu[i];
      }
      psol[k] = exp(u);
    }

    sc = 0; cmax = 0;
    for (int i=1; i<=ntot; i++) {
      if (nch[i] > 0 and i != 11) {
        sc = sc + tot[i] * nch[i];
        if (tot[i] * nch[i] >= cmax) {
          cmax = tot[i] * nch[i];
          icat = i;
        }
      }
    }

    sa = 0; amax=0;
    for (int i=1; i<=ntot; i++) {
      if (nch[i] < 0) {
        sa = sa + tot[i] * -nch[i];
        if (tot[i] * (-nch[i]) >= amax) {
          amax = tot[i] * (-nch[i]);
          iani = i;
        }
      }
    }

    if (sc + sa != 0) dca = 200 * abs(sc - sa) / (sc + sa);
    delta = sc - sa;
    cout << "sum cations = " << sc << endl;
    cout << "sum of anions = " << sa << endl;
    cout << "Electrical balance = " << (dca * 100 + 0.5) / 100 << "%" << endl;
    tot[icat] = tot[icat] - delta / 2 / nch[icat];
    tot[iani] = tot[iani] + delta / 2 / (-nch[iani]);
    for (int i=1; i<=12; i++) {
      tot0[i] = tot[i];
    }

    {
      string x = "n/a";
      mineral0_S.push_back(x);

      file.open("murtf3");
      getline(file, line);
      stringstream linestream(line);
      string value;
      getline(linestream, value, ',');
      nc = stoi(value);
      getline(linestream, value, ',');
      na = stoi(value);
      getline(linestream, value, ',');
      nm0 = stoi(value);
      
      for (int i=1; i<=(nc+na+1); i++) {
        getline(file, line);
      }
      for (int k=1; k<=nm0; k++) {
        getline(file, line);
        stringstream linestream(line);
        string value;
        getline(linestream, value,  ',');
        mineral0_S.push_back(value);
      }
      for (int k=1; k<=nm; k++) {
        nwmin[k] = 0;
      }
      for (int k=1; k<=nm; k++) {
        for (int l=1; l<=nm0; l++) {
          if (mineral0_S.at(l) == mineral_S.at(k)) {
            nwmin[k] = 1;
          }
        }
      }
      file.close();
    }

    min_S = "murtf3";

    // run the main eql loop
    // ...loop conditions:
    // ...500: calculate molalities
    // ...200: iterate activities
    // ...400: screen output
    // ...error=true: exit/error
    error = false;
    simulation_loop = 500;
    while (!error) {
      if (simulation_loop == 500) {
        molal[1] = tot[1]; molal[2] = tot[2]; molal[3] = tot[3];
        molal[6] = tot[6]; molal[8] = tot[8]; molal[11] = tot[11];
        molal[13] = psc[1] / molal[11];
        if (tot[9] > 0) {
          a = 1 + molal[13] / psc[7];
          b = 3 * psc[8] * molal[13];
          c = 4 * psc[9] * pow(molal[13], 2);
          xu = tot[9] / 2;
          double u = xu;
          //do while...
          while (true) {
            eq = a * xu + b * pow(xu, 3) + c * pow(xu, 4);
            if (200 * abs(eq - tot[9]) / (eq + tot[9]) < pk) break;
            
            u = u / 2;
            if (eq > tot[9]) {
              xu -= u;
            } else {
              xu += u;
            }
          }
          molal[9] = xu;
          molal[19] = molal[9] * molal[13] / psc[7];
          molal[20] = molal[13] * pow(molal[9], 3) * psc[8];
          molal[21] = pow(molal[13], 2) * pow(molal[9], 4) * psc[9];
        }
        molal[14] = ((tot[12] + molal[11] - molal[13] - 
                      molal[19] - molal[20] - 2 * molal[21]) / 
                     (2 + molal[11] / psc[2]));
        molal[12] = ((tot[12] + molal[11] - molal[13] -
                      molal[19] - molal[20] - 2 * molal[21]) /
                     (1 + 2 * psc[2] / molal[11]));
        molal[15] = (molal[12] * molal[11] / psc[3]);
        molal[4] = (tot[4] / (1 + molal[14] / psc[4] + molal[19] / psc[10]));
        molal[16] = (molal[4] * molal[14] / psc[4]);
        molal[22] = (molal[4] * molal[19] / psc[10]);
        molal[5] = tot[5] / (1 + molal[14] / psc[5] + 
                             molal[13] / psc[6] + 
                             molal[19] / psc[11]);
        molal[17] = molal[5] * molal[14] / psc[5];
        molal[18] = molal[5] * molal[13] / psc[6];
        molal[23] = molal[5] * molal[19] / psc[11];
        molal[10] = tot[10] / (1 + psc[12] / molal[11]);
        molal[24] = tot[10] / (1 + molal[11] / psc[12]);
        molal[7] = tot[7] / (1 + molal[11] / psc[13]);
        molal[25] = molal[7] * molal[11] / psc[13];

        // charge balance
        sc = 0; cmax = 0;
        for (int i=1; i<=n; i++) {
          if (nch[i] > 0) {
            sc += molal[i] * nch[i];
            if (molal[i] * nch[i] > cmax) {
              cmax = molal[i] * nch[i];
              icat = i;
            }
          }
        }
        sa = 0; amax = 0;
        for (int i=1; i<=n; i++) {
          if (nch[i] < 0) {
            sa += molal[i] * (-nch[i]);
            if (molal[i] * (-nch[i]) > amax) {
              amax = molal[i] * (-nch[i]);
              iani = i;
            }
          }
        }
        delta = sc - sa;
        molal[icat] -= delta / 2 / nch[icat];
        molal[iani] += delta / 2 / (-nch[iani]);

        sc = (molal[1] + molal[2] + molal[3] + molal[4] * 2 + 
              molal[5] * 2 + molal[11] + molal[18] + molal[22] + 
              molal[23]);
        sa = (molal[6] + molal[7] * 2 + molal[8] + molal[12] + 
              molal[13] + molal[14] * 2 + molal[19] + molal[20] + 
              molal[21] * 2 + molal[24] + molal[25]);

        cout << endl;
        cout << "Sum of cations = " << sc << " corrected for " << aq_S[icat] << endl;
        cout << "Sum of anions = " << sa << " corrected for " << aq_S[iani] << endl;
        cout << endl;

        float s = 0;
        for (int i=1; i<=n; i++) {
          s += molal[i] * atom[i];
        }
        if (units == "molar") {
          ee = 1000 / (1000 * dens - s);
          for (int i=1; i<=ntot; i++) {
            if (i != 11) {
              tot[i] = tot[i] * ee;
            }
          }
          for (int i=1; i<=n; i++) {
            if (i != 11) {
              molal[i] = molal[i] * ee;
            }
          }
        } else if (units == "molal") {
          ee = 1;
        }
        stdi = s * ee;

        simulation_loop = 200;

      }
      if (simulation_loop == 200) {
        while (true) {
          int nu = 1;
          int ncompt = 0;
          while (nu != 0) {
            eql_actp();
            for (int i=1; i<=n; i++) {
              act[i] = molal[i] * gact[i];
            }
            act[0] = aw;
            tot[11] = pow(10, -ph) / gact[11];
            act[11] = pow(10, -ph);
            for (int i=1; i<=12; i++) {
              for (int j=1; j<=n; j++) {
                if (molal[j] !=0 ) z[i][j] = kmat[i][j];
              }
              double u = 0;
              for (int j=1; j<=n; j++) {
                u += kmat[i][j] * molal[j];
              }
              zz[i] = tot[i] - u;
            }

            for (int i=13; i<=n; i++) {
              for (int j=1; j<=n;+) {
                if (molal[j] != 0) {
                  z[i][j] = kmat[i][j] / molal[j];
                } else if (molal[j] == 0) {
                  z[i][j] = 0;
                }
              }
              double u = 0;
              for (int j=0; j<=n; j++) {
                if (act[j] > 0) {
                  u += kmat[i][j] * log(act[j]);
                }
              }
              zz[i] = log(psc[i-12]) - u;
            }
            // line 409
          }
        }
      }
      simulation_loop = 500;
      error = true;
    }
    

  } 

  // seconrdary simulation methods
  void eql_actp() {
  
  }

  void density() {

  }

  void invar() {

  }

  double g0() {
    return 0;
  }

  double g1() {
    return 0;
  }

  double j0() {
    return 0;
  }

  double j1() {
    return 0;
  }

  double temperature(double at, double bt, double ct, double dt, double et) {
    return (at + bt * temp + ct * pow(temp, 2) + dt * pow(temp,3) + et * pow(temp, 4));
  }

};





int main() {
  
  Simulation test("example.txt");
  test.run_eql();

  // Successful exit
  return 0;
}
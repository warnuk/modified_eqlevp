#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <regex>

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

std::string ltrim(const std::string &s) {
    return std::regex_replace(s, std::regex("^\\s+"), std::string(""));
}
 
std::string rtrim(const std::string &s) {
    return std::regex_replace(s, std::regex("\\s+$"), std::string(""));
}
 
std::string trim(const std::string &s) {
    return ltrim(rtrim(s));
}

double g0(double x) {
    return 2 * (1 - (1+x) * exp(-x)) / pow(x, 2);
}

double g1(double x) {
    return -2 * (1 - (1 + x + pow(x, 2) / 2) * exp(-x)) / pow(x, 2);
}

double j0(double x) {
    double ya, yb, yc, yd;
    ya = 4.581; yb = -0.7237; yc = -0.012; yd = 0.528;
    return x / (4 + ya * pow(x, yb) * exp(yc * pow(x, yd)));
}

double j1(double x) {
    double ya, yb, yc, yd;
    ya = 4.581; yb = -0.7237; yc = -0.012; yd = 0.528;
    return ((4 + ya * pow(x, yb) * (1 - yb - yc * yd * pow(x, yd)) * 
            exp(yc * pow(x, yd))) / pow((4 + ya * pow(x, yb) * exp(yc * pow(x, yd))), 2));
}

double temperature(double at, double bt, double ct, double dt, double et, double temp) {
    return (at + bt * temp + ct * pow(temp, 2) + dt * pow(temp,3) + et * pow(temp, 4));
}

// define constants
const int n = 25; const int ntot = 12; const int ncomplex = 14; 
const int n0 = 10; const int max_cat = 5; const int max_an = 5; 
const int max_constit = 30; const int max_conv = 100; 
const double epsilon = 1e-8; const double eps = 1e-12;
const double pk0 = 0.1; const double pkf = 0.0001;
const double pkstd = 0.1;


class Simulation {
    public:

    // declare arrays and vectors
    array<string,max_constit+1> constit_S; array<string,n+1> aq_S;
    double tot[ntot+1] = {0}; double tot0[ntot+1] = {0}; 
    double totinit[ntot+1] = {0}; double molal[n+1] = {0}; 
    double act[n+1] = {0}; double gact[n+1] = {0}; 
    double atom[n+1] = {0}; double xx[n+1];
    double zz[n+1] = {0}; double z[n+1][n+1] = {0};
    double psc[15] = {0}; int nch[n+1] = {0};
    double cat[max_cat] = {0}; double ani[max_an] = {0};
    int nchcat[max_cat+1] = {0}; int nchani[max_an] = {0};
    int kmat[n+1][n+1] = {0}; int ica[n+1] = {0};

    // declare variables
    int i, j, k, l, ii, nminer, nm, nt, nc, na, ncomp, icat, iani, nconstit;
    int nm0, nu, ni, nj, num, nwmp, nwm, kinvariant, nw, nmneuf, npasfich, npasimp, ncompt, npasecran;
    string num_S, unit_S, x_S, e_S, repinit_S, min_S, pco_S, rep_S;
    string zone_S, m_S, p_S, y_S, fichmin_S, molemin_S, fich_S, prnt_S, constituant_S;
    string at_S, bt_S, ct_S, dt_S, et_S;
    double temp, dens, ph, at, bt, ct, dt, et, pkmol, pkeq;
    double ap0, bp0;
    double a, b, c, u, sc, cmax, sa, amax, dca, delta, xu, eq, s, ee, stdi, aw, fi, std, ef, tinit, stdmax;
    double po, po0, poinit, phinit, ph0, pc, alcar, albor, alsil, aloh, alh, altest, d, v, xi, eps;

    string label, tra_file, log_file, event_file, chem_file, min_file;
    string pco2, syst_S; string output_units, units;
    int print_step, output_step, verbose, output;
    double max_sal, dilute;
    
    ifstream file; 
    ofstream outfile;
    string line;

    int nconv = 2; int ndepact = 0; string pco2_S = "";
    double pk = 0.1; double dph = 0.2; double mh2o = 55.51;
    double dil = 0; double diltot = 1; double poa = 0;

    // Declare vectors
    vector<string> add_min; vector<string> rem_min;
    vector<string> nom_ion_S, mineral_S, mineral0_S;
    vector<double> c_ion, mu, mum, psol, pai;
    vector<int> kinvar, lmin, nwmin;
    vector<vector<double>> wmin;
    
    
    

    // Declare vectors for activity calculation
    vector<int> nzc, nza;
    vector<vector<double>> b0, b1, b2, c0;
    vector<vector<double>>  tcp, tap, lcp, lap;
    vector<vector<vector<double>>> scp, sap, xip;

    Simulation(string fname) {
        run_eql(fname);
    }

    private:

        void run_eql(string fname) {
            // read user inputs
            read_inputs(fname);

            if (verbose == 1) {
                time_t now;
                time (&now);
                cout << "\nThis is EQL..............\n";
                cout << asctime(localtime(&now)) << "\n";
            }
            
            // files
            tra_file = label + ".tra"; log_file = label + ".log";
            event_file = label + ".j" + syst_S + "@";
            chem_file = label + ".j" + syst_S + "&";
            min_file = label + ".j" + syst_S + "%";

            // set initial values
            tinit = temp;

            // Set nminer based on nonzero input chemistry
            nminer = 1;
            for (int i = 1; i <= ntot; i++) {
                if (tot[i] > 0) nminer += 1;
            }

            // Interpret units based on reported density
            if (dens == 0) dens = 1.e0;
            if (dens == 1.e0) {
                units = "molal";
            } else {
                units = "molar";
            }

            // Open and read "aqu.dat"
            file.open("aqu.dat");
            for (int i = 0; i <= n; i++) {
                getline(file, line);
                stringstream linestream(line); string value;
                getline(linestream, value, ','); aq_S[i] = trim(value);
                getline(linestream, value, ','); atom[i] = stod(value);
                getline(linestream, value, ','); nch[i] = stoi(value);
            }
            file.close();

            // Set default ica to 1
            for (int i=1; i<=n; i++) {
                ica[i] = 1;
            }

            // Open and read "matrice1"
            file.open("matrice1");
            for (int i=1; i<=n; i++) {
                getline(file, line);
                stringstream linestream(line); string value;
                for (int j=1; j<=n; j++) {
                    getline(linestream, value, ','); kmat[i][j] = stoi(value);
                }
            }
            file.close();
            kmat[13][0] = -1; kmat[15][0] = -1; kmat[20][0] = 3; kmat[21][0] = 5;

            // Open and read "complex 3"
            file.open("complex3");
            for (int i=1; i<=14; i++) {
                getline(file, line);
                stringstream linestream(line); string value;
                //double at, bt, ct, dt, et;
                
                getline(linestream, value, ','); x_S = trim(value);
                getline(linestream, value, ','); at = stod(value);
                getline(linestream, value, ','); bt = stod(value);
                getline(linestream, value, ','); ct = stod(value);
                getline(linestream, value, ','); dt = stod(value);
                getline(linestream, value, ','); et = stod(value);

                psc[i] = pow(10, (at + bt/300*temp + ct/30000*pow(temp, 2) + dt/3000000*pow(temp,3) + et/300000000*pow(temp,4)));
            }
            file.close();

            {
                // Open and read murtf2
                file.open("murtf2");
                getline(file, line);
                stringstream linestream(line); string value;
                getline(linestream, value, ','); nc = stoi(value);
                getline(linestream, value, ','); na = stoi(value);
                getline(linestream, value, ','); nm = stoi(value);
                nt = nc + na;
                file.close();
            }

            // Initialize arrays of zeros + empty string arrays
            for (int i=0; i<=nt; i++) {
                nom_ion_S.push_back(""); c_ion.push_back(0); mu.push_back(0); kinvar.push_back(0);
            }
            for (int i=0; i<=nm; i++) {
                mineral_S.push_back(""); mineral0_S.push_back("");
                mum.push_back(0); psol.push_back(0); pai.push_back(0);
                lmin.push_back(0); nwmin.push_back(0);

                vector<double> temp1d;
                for (int k=0; k<=nt; k++) temp1d.push_back(0);
                wmin.push_back(temp1d);
            }

            {
                // Open and read murtf2
                file.open("murtf2");

                // Ignore the first line
                getline(file, line);
            
                // Iterate through the chemical species
                for (int i=0; i<=14; i++) {
                    getline(file, line);
                    stringstream linestream(line); string value, x_S;
                    double at, bt, ct, dt, et;
                    
                    getline(linestream, value, ','); x_S = upper_to_lower(trim(value));
                    getline(linestream, value, ','); at = stod(value);
                    getline(linestream, value, ','); bt = stod(value);
                    getline(linestream, value, ','); ct = stod(value);
                    getline(linestream, value, ','); dt = stod(value);
                    getline(linestream, value, ','); et = stod(value);
            
                    for (int j=0; j<=nt; j++) {
                        if (x_S == aq_S[j]) {
                            mu[j] = (at + bt / 300 * temp + ct / 30000 * pow(temp, 2) + dt / 3000000 * pow(temp, 3) + et / 300000000 * pow(temp, 4));
                        }
                    }
                }

                for (int k=1; k<=nm; k++) {
                    getline(file, line);
                    stringstream linestream(line); string value, x_S;
                    //int ncomp; double c_ion;
            
                    getline(linestream, value, ','); mineral_S[k] = trim(value);
                    getline(linestream, value, ','); ncomp = stoi(value);
            
                    for (int i=1; i<=ncomp; i++) {
                        getline(linestream, value, ','); c_ion[i] = stod(value);
                        getline(linestream, value, ','); nom_ion_S[i] = trim(value);
                        x_S = upper_to_lower(nom_ion_S[i]);

                        for (int j=0; j<=nt; j++) {
                            if (x_S == aq_S[j]) {
                                wmin[k][j] = c_ion[i];
                            }
                        }
                    }

                    getline(linestream, value, ','); at = stod(value);
                    getline(linestream, value, ','); bt = stod(value);
                    getline(linestream, value, ','); ct = stod(value);
                    getline(linestream, value, ','); dt = stod(value);
                    getline(linestream, value, ','); et = stod(value);

                    mum[k] = (at + bt / 300 * temp + ct / 30000 * pow(temp, 2) + dt / 3000000 * pow(temp, 3) +  et / 300000000 * pow(temp, 4));
                } 
                file.close();
            }

            // Set psol
            for (int k=1; k<=nm; k++) {
                u = mum[k];
                for (int i=0; i<=nt; i++) {
                    u = u - wmin[k][i] * mu[i];
                }
                psol[k] = exp(u);
            }

            // Charge balance
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
            // Open and read mineral database
            {
                file.open("murtf3");
                getline(file, line);
                stringstream linestream(line); string value;
                getline(linestream, value, ','); nc = stoi(value);
                getline(linestream, value, ','); na = stoi(value);
                getline(linestream, value, ','); nm0 = stoi(value);

                for (int i=1; i<=(nc+na+1); i++) {
                    getline(file, line);
                }

                for (int k=1; k<=nm0; k++) {
                    getline(file, line);
                    stringstream linestream(line); string value;
                    getline(linestream, value, ','); mineral0_S[k] = trim(value);
                }
                for (int k=1; k<=nm; k++) {
                    nwmin[k] = 0;
                }
                for (int k=1; k<=nm; k++) {
                    for (int l=1; l<=nm0; l++) {
                        if (mineral0_S[l] == mineral_S[k]) nwmin[k] = 1;
                    }
                }
                file.close();
            }
            min_S = "murtf3";

            LOOP500:
                molal[1] = tot[1]; molal[2] = tot[2]; molal[3] = tot[3];
                molal[6] = tot[6]; molal[8] = tot[8]; 
                tot[11] = pow(10, -ph);
                molal[11] = tot[11];
                molal[13] = psc[1] / molal[11];
                if (tot[9] > 0) {
                    a = 1 + molal[13] / psc[7];
                    b = 3 * psc[8] * molal[13];
                    c = 4 * psc[9] * pow(molal[13], 2);
                    xu = tot[9] / 2;
                    u = xu;
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
                molal[14] = (tot[12] + molal[11] - molal[13] - molal[19] - molal[20] - 2 * molal[21]) / (2 + molal[11] / psc[2]);
                molal[12] = (tot[12] + molal[11] - molal[13] - molal[19] - molal[20] - 2 * molal[21]) / (1 + 2 * psc[2] / molal[11]);
                molal[15] = molal[12] * molal[11] / psc[3];
                molal[4] = tot[4] / (1 + molal[14] / psc[4] + molal[19] / psc[10]);
                molal[16] = molal[4] * molal[14] / psc[4];
                molal[22] = molal[4] * molal[19] / psc[10];
                molal[5] = tot[5] / (1 + molal[14] / psc[5] + molal[13] / psc[6] +  molal[19] / psc[11]);
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

                s = 0;
                for (int i=1; i<=n; i++) {
                    s += molal[i] * atom[i];
                }
                if (units == "molar") {
                    ee = 1000 / (1000 * dens - s);
                    for (int i=1; i<=ntot; i++) {
                        if (i != 11) tot[i] = tot[i] * ee;
                    }
                    for (int i=1; i<=n; i++) {
                        if (i != 11) molal[i] = molal[i] * ee;
                    }
                } else if (units == "molal") {
                    ee = 1;
                }
                stdi = s * ee;
                
            LOOP200:
                while (true) {
                    nu = 1;
                    ncompt = 0;
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
                                if (molal[j] != 0 ) z[i][j] = kmat[i][j];
                            }
                            u = 0;
                            for (int j=1; j<=n; j++) {
                                u += kmat[i][j] * molal[j];
                            }
                            zz[i] = tot[i] - u;
                        }

                        for (int i=13; i<=n; i++) {
                            for (int j=1; j<=n; j++) {
                                if (molal[j] != 0) {
                                    z[i][j] = kmat[i][j] / molal[j];
                                } else if (molal[j] == 0) {
                                    z[i][j] = 0;
                                }
                            }
                            u = 0;
                            for (int j=0; j<=n; j++) {
                                if (act[j] > 0) {
                                    u += kmat[i][j] * log(act[j]);
                                }
                            }
                            zz[i] = log(psc[i-12]) - u;
                        }

                        for (int k=1; k<=n0; k++) {
                            if (tot[k] == 0 and k != 12) {
                                ica[k] = 0;
                                for (int i=k+1; i<=n; i++) {
                                    if (kmat[i][k] != 0) ica[i] = 0;
                                }
                                for (int j=k+1; j<=n; j++) {
                                    if (kmat[k][j] != 0) ica[j] = 0;
                                }
                            }
                        }

                        ni = n; nj = n;
                        for (int k=n; k>=1; --k) {
                            if (ica[k] == 0) {
                                for (int i=k; i<=ni-1; i++) {
                                    for (int j=1; j<=nj; j++) {
                                        z[i][j] = z[i+1][j];
                                    }
                                    zz[i] = zz[i+1];
                                }
                                ni -= 1;
                                for (int j=k; j<=nj-1; j++) {
                                    for (int i=1; i<=ni; i++) {
                                        z[i][j] = z[i][j+1];
                                    }
                                }
                                nj -= 1;
                            }
                        }

                        for (int k=2; k<=ni; k++) {
                            for (int i=k; i<=ni; i++) {
                                if (z[i][k-1] != 0) {
                                    u = z[k-1][k-1] / z[i][k-1];
                                    for (int j=k; j<=ni; j++) {
                                        z[i][j] = z[k-1][j] - z[i][j] * u;
                                    }
                                    zz[i] = zz[k-1] - zz[i] * u;
                                }
                            }
                        }

                        xx[ni] = zz[ni] / z[ni][ni];
                        for (int i=ni-1; i>=1; --i) {
                            s = 0;
                            for (int j=i+1; j<=ni; j++) {
                                s += z[i][j] * xx[j];
                            }
                            xx[i] = (zz[i] - s) / z[i][i];
                        }

                        for (int k=1; k<=n; k++) {
                            if (ica[k] == 0) {
                                for (int i=ni; i>= k; --i) {
                                    xx[i+1] = xx[i];
                                }
                                xx[k] = 0;
                                ni += 1;
                            }
                        }

                        ncompt += 1;
                        cout << "iteration molalities " << ncompt << endl;
                        if (ncompt >= 100) {
                            for (int i=1; i<=n; i++) {
                                if (molal[i] + xx[i] / nconv < 0) {
                                    cout << "the equation set diverges: end of program" << endl;
                                    stop_simulation();
                                }
                            }
                        }
                        for (int i=1; i<=n; i++) {
                            if (molal[i] + xx[i] / nconv < 0) {
                                molal[i] = eps;
                            } else {
                                molal[i] += xx[i] / nconv;
                            }
                        }

                        nu = 0;
                        for (int i=1; i<=n; i++) {
                            if (ica[i] == 1) {
                                if (200 * abs(xx[i] / nconv / 
                                    (2 * molal[i] - xx[i] / nconv)) > pk) {
                                        nu = 1;
                                }
                            }
                        }
                    }

                    std = 0;
                    for (int i=0; i<=n; i++) {
                        std += molal[i] * atom[i];
                    }
                    cout << "tdsi = " << stdi << endl;
                    cout << "tds = " << std << endl;

                    if (abs(std - stdi) / (std + stdi) * 200 < pkstd) {
                        break;
                    } else {
                    if (units == "molar") {
                        ef = (1000 + std) / dens / 1000;
                        for (int i=1; i<=ntot; i++) {
                            if (i != 11) tot[i] = tot[i] / ee * ef;
                        }
                        for (int i=0; i<=n; i++) {
                            if (i != 11) molal[i] = molal[i] / ee * ef;
                        }
                        ee = ef;
                    }
                    stdi = std;
                    }
                }

                if (units == "molal" and dil == 0) {
                    for (int i=1; i<=5; i++) {
                        cat[i] = tot[i];
                        nchcat[i] = nch[i];
                    }
                    for (int j=1; j<=3; j++) {
                        ani[j] = tot[j+5];
                        nchani[j] = -nch[j+5];
                    }
                    ani[4] = molal[12]; ani[5] = molal[14] + molal[16] + molal[17];
                    nchani[4] = -nch[12]; nchani[5] = -nch[14];
                    eql_density();
                }

                po = log10(act[15] / psc[14]);
                if (pco2_S == "") {
                    if (diltot == 1) {
                        poinit = po;
                        phinit = ph;
                    }
                    if (verbose == 1) {
                        cout << "LOG PCO2 = " << po << endl;
                    }

                    if (pco2 != "na" or pco2 != "") {
                        pc = stod(pco2);
                        po0 = po; ph0 = ph;
                        pco2_S = "y";
                    } else if (pco2 == "na" or pco2 == "") {
                        pco2_S = "n";
                    }
                }

                if (pco2_S == "y") {
                    if (verbose == 1) {
                        cout << endl;
                        cout << "Log(PCO2) selected   = " << pc << endl;
                        cout << "Log(PCO2) calculated = " << po << endl;
                        cout << endl;
                    }
                    if (abs(po - pc) > 0.01) {
                        if (po < pc and poa < pc) ph = ph - dph;
                        if (po < pc and poa > pc) {
                            dph = dph / 2; ph = ph - dph;
                        }
                        if (po > pc and poa > pc) ph = ph + dph;
                        if ((po > pc) and (poa < pc)) {
                            dph = dph / 2; ph = ph + dph;
                        }
                        poa = po;
                        goto LOOP200;
                    }
                }
                if (pk > pkf) {
                    pk = pkf;
                    if (verbose == 1) {
                        cout << "last iteration" << endl;
                    }
                    goto LOOP200;
                }

            LOOP400:
                for (int i=1; i<=n0; i++) {
                    totinit[i] = 0;
                    for (int j=1; j<=n; j++) {
                        totinit[i] += molal[j] * kmat[i][j];
                    }
                }
                if (verbose == 1) {
                    cout << endl;
                    cout << label << endl;
                    cout << setw(10) << " " << setw(15) << "MOLALITY" << setw(15) << "ACT COEFF" << setw(15) << "ACTIVITY" << setw(15) << "MOLAL TOT" << endl;
                    for (int i=1; i<=n0; i++) {
                        if (molal[i] != 0) {
                            cout << setw(10) << aq_S[i] << setw(15) << molal[i] << setw(15) << gact[i] << setw(15) << act[i] << setw(15) << totinit[i] << endl;
                        }
                    }
                    for (int i=n0+1; i<=n; i++) {
                        if (molal[i] != 0) {
                            cout << setw(10) << aq_S[i] << setw(15) << molal[i] << setw(15) << gact[i] << setw(15) << act[i] << endl;
                        }
                    }
                    cout << endl;
                }

                alcar = molal[12] + 2 * (molal[14] + molal[16] + molal[17]);
                albor = molal[19] + molal[20] + 2 * molal[21] + molal[22] + molal[23];
                alsil = molal[24];
                aloh = molal[13] + molal[18];
                alh = -molal[11] - molal[25];
                altest = alcar + albor + alsil + aloh + alh;

                if (verbose == 1) {
                    cout << "ELECTRICAL BALANCE     = " << dca << "% corrected on " << aq_S[icat] << " and " << aq_S[iani] << endl;
                    cout << "TOTAL DISSOLVED SOLIDS = " << std << "g/kg(H2O)" << endl;
                    cout << "MOLAL/MOLAR FACTOR     = " << (1000 * dens / (1000 + std)) << endl;
                    cout << "DENSITY                = " << dens << endl;
                    cout << "IONIC STRENGTH         = " << fi << endl;
                    cout << "WATER ACTIVITY         = " << aw << endl;
                    if (diltot > 1) {
                        cout << "DILUTION               = " << diltot << endl;
                    }
                    if ((pco2_S == "") || (pco2_S == "n")) {
                        cout << "pH                     = " << ph << endl;
                        cout << "LOG PCO2               = " << po << endl;
                    } else if (pco2_S == "y") {
                        cout << "INITIAL LOG PCO2       = " << poinit << endl;
                        cout << "INITIAL pH             = " << phinit << endl;
                        cout << "CALCULATED LOG PCO2    = " << po << endl;
                        cout << "CALCULATE pH           = " << ph << endl;
                    }
                    cout << "CARBONATE ALKALINITY   = " << alcar << endl;
                    cout << "BORATE ALKALINITY      = " << albor << endl;
                    cout << "SILICATE ALKALINITY    = " << alsil << endl;
                    cout << "OH ALKALINITY          = " << aloh << endl;
                    cout << "H ALKALINITY           = " << alh << endl;
                    cout << endl;
                    cout << "TOTAL ALKALINITY       = " << altest << endl;
                    cout << "             init. alk = " << tot[12] << endl;

                    //LINE 638
                    cout << endl;
                    cout << setw(18) << left << "TESTS" << setw(18) << "SUM OF SPECIES" << setw(18) << "INIT. CONC." << setw(18) << "BALANCE %" << endl;
                    cout << endl;
                    for (int i=1; i<=12; i++) {
                        if (ica[i] == 1) {
                            u = 0;
                            for (int j=1; j<=n; j++) {
                                u += molal[j] * kmat[i][j];
                            }
                            if (u + tot[i] != 0) d = 200 * abs(u - tot[i]) / (u + tot[i]);
                            if (i == 12) {
                                if (u + tot[i] != 0 and tot[i] != 0) {
                                    cout << setw(18) << left << "ALK" << setw(18) << u << setw(18) << tot[i] << setw(18) << d << endl;
                                } else {
                                    cout << setw(18) << left << "ALK" << setw(18) << u << setw(18) << tot[i] << endl;
                                }
                            } else {
                                zone_S = lower_to_upper(aq_S[i]);
                                cout << setw(18) << left << zone_S << setw(18) << u << setw(18) << tot[i] << setw(18) << d << endl;
                            }
                        }
                    }
                    cout << endl << endl;
                    cout << setw(18) << " " << setw(18) << "LOG(IAP)" << setw(18) << "LOG(K)" << setw(18) << "BALANCE %" << endl;
                    cout << endl;
                }

                for (int i=13; i<=n; i++) {
                    u = 0;
                    if (ica[i] == 1) {
                        for (int j=0; j<=n; j++) {
                            if (act[j] != 0) u += kmat[i][j] * log10(act[j]);
                        }
                        v = log10(psc[i-12]);
                        d = 200 * abs(u - v) / (u + v);
                        if (i == 13) {
                            cout << setw(18) << "H2O" << setw(18) << u << setw(18) << v << setw(18) << d << endl;
                        } else {
                            zone_S = lower_to_upper(aq_S[i]);
                            cout << setw(18) << zone_S << setw(18) << u << setw(18) << v << setw(18) << d << endl;
                        }
                    }
                }

                for (int k=1; k<=nm; k++) {
                    pai[k] = 1;
                    for (int i=0; i<=nt; i++) {
                        pai[k] = pai[k] * pow(act[i], wmin[k][i]); 
                    }
                }

                if (verbose == 1) {
                    cout << endl;
                    cout << setw(18) << " " << setw(18) << "SOLUB. PROD." << setw(18) << "ION. ACT. PROD." << setw(18) << "SATUR. RATIO" << endl;
                    cout << endl;
                    nwm = 0; nwmp = 0;
                    for (int k=1; k<=nm; k++) {
                        if (pai[k] != 0) {
                            if (nwmin[k] == 0) zone_S = " " + mineral_S[k];
                            if (nwmin[k] == 1) zone_S = "*" + mineral_S[k];
                            x_S = " ";
                            if (pai[k] / psol[k] >= 1 and nwmin[k] == 1) {
                                nwm += 1;
                                x_S = "*";
                            } else if (pai[k] / psol[k] >= 0.9 and pai[k] / psol[k] < 1 and nwmin[k] == 1) {
                                nwmp += 1;
                            }
                            cout << setw(18) << zone_S << setw(18) << psol[k] << setw(18) << pai[k] << setw(18) << pai[k] / psol[k] << left << x_S << endl;
                        }
                    }
                    cout << endl;
                    
                }

                if (output == 1) {
                    if (verbose == 1) cout << "LOG FILE IS " << log_file << endl;

                    outfile.open(log_file, ios::trunc);
                    outfile << left << setw(15) << " " << left << setw(20) << "MOLALITY" << left << setw(20) << "ACT. COEFF." << left << setw(20) << "ACTIVITY" << endl;
                    outfile << endl;

                    for (int i=1; i<=n; i++) {
                        if (molal[i] != 0) {
                            outfile << left << setw(15) << aq_S[i] << left << setw(20) << molal[i] << left << setw(20) << gact[i] << left << setw(20) << act[i] << endl;
                        }
                    }

                    outfile << endl;

                    outfile << "ELECTRICAL BALANCE     = " << dca << "% corrected on " << aq_S[icat] << " and " << aq_S[iani] << endl;
                    outfile << "TOTAL DISSOLVED SOLITS = " << std << "g/kg(H2O)" << endl;
                    outfile << "MOLAL/MOLAR FACTOR     = " << (1000 * dens / (1000 + std)) << endl;
                    outfile << "DENSITY                = " << dens << endl;
                    outfile << "IONIC STRENGTH         = " << fi << endl;
                    outfile << "WATER ACTIVITY         = " << aw << endl;
                    if (diltot > 1) {
                        outfile << "DILUTION               = " << diltot << endl;
                    }
                    if (pco2_S == "" or pco2_S == "n") {
                        outfile << "pH                     = " << ph << endl;
                        outfile << "LOG PCO2               = " << po << endl;
                    } else if (pco2_S == "y") {
                        outfile << "INITIAL LOG PCO2       = " << poinit << endl;
                        outfile << "INITIAL pH             = " << phinit << endl;
                        outfile << "CALCULATED LOG PCO2    = " << po << endl;
                        outfile << "CALCULATED pH          = " << ph << endl;
                    }
                    outfile << "CARBONATE ALKALINITY   = " << alcar << endl;
                    outfile << "BORATE ALKALINITY      = " << albor << endl;
                    outfile << "SILICATE ALKALINITY    = " << alsil << endl;
                    outfile << "OH ALKALINITY          = " << aloh << endl;
                    outfile << "H ALKALINITY           = " << alh << endl;
                    outfile << "TOTAL ALKALINITY       = " << setw(25) << altest << "init. alk. = " << tot[12] << endl;

                    outfile << endl;
                    outfile << setw(18) << " " << setw(18) << "SOLUB. PROD." << setw(18) << "ION. ACT. PROD." << setw(18) << "SATUR. RATIO" << endl;
                    outfile << endl;
                    for (int k=1; k<=nm; k++) {
                        if (pai[k] > 0) {
                            outfile << setw(18) << mineral_S[k] << setw(18) << psol[k] << setw(18) << pai[k] << setw(18) << pai[k] / psol[k] << endl;
                        }
                    
                    }
                    outfile.close();
                }

                for (int k=1; k<=nm; k++) {
                    lmin[k] = 0;
                }
                if (verbose == 1) {
                    cout << endl;
                    cout << "The initial solution is oversaturated in " << nwm << " mineral(s) of the data base MURTF3:" << endl;
                }
                for (int k=1; k<=nm; k++) {
                    if (pai[k] / psol[k] >= 1 and nwmin[k] == 1) {
                        if (verbose == 1) {
                            cout << setw(20) << mineral_S[k] << pai[k] / psol[k] << endl;
                        }
                        lmin[k] = 1;
                    }
                }
                if (nwm > nminer) {
                    if (verbose == 1) {
                        cout << endl;
                        cout << "VIOLATION OF THE PHASE RULE :" << endl;
                        cout << "The maximum number of minerals allowed is: " << nminer << endl;
                        cout << "The evaporation program cannot start with this paragenesis" << endl;
                    }
                } else if (nwm <= nminer) {
                    if (verbose == 1) cout << endl;
                    eql_invar();
                    if (kinvariant > 0) {
                        if (verbose == 1) {
                            cout << "The activity of water is constrained by: " << endl;
                            for (int k=1; k<=kinvariant; k++) {
                                cout << mineral_S[kinvar[k]] << endl;
                            }
                            cout << endl;
                        }
                    } else if (kinvariant == -1) {
                        if (verbose == 1) {
                            cout << "System in thermodynamic desequilibrium" << endl;
                            cout << "The activity of water is constrained at different values" << endl;
                            cout << "by more than one mineral assemblage" << endl;
                        }
                    } else if (kinvariant == -2) {
                        if (verbose == 1) {
                            cout << "System in thermodynamic desequilibrium: inconsistent mineral assemblage" << endl;
                        }
                    } else if (kinvariant == 0) {
                        if (verbose == 1) {
                            cout << "No invariant paragenesis detected" << endl;
                        }
                    }
                

                    if (kinvariant != 0) {
                        if (verbose == 1) {
                            cout << "The evaporation program cannot start with this paragenesis" << endl;
                        }
                    }

                    if (kinvariant == 0 and nwmp > 0) {
                        if (verbose == 1) {
                            cout << endl;
                            cout << "The solution is close to saturation in " << nwmp << " mineral(s) of the data base MURTF3:" << endl;
                        }

                        for (int k=1; k<=nm; k++) {
                            if (pai[k] / psol[k] >= 0.9 and pai[k] / psol[k] < 1 and nwmin[k] == 1) {
                                if (verbose == 1) {
                                    cout << setw(20) << mineral_S[k] << pai[k] / psol[k] << endl;
                                }
                                lmin[k] = 1;
                            }
                        }

                        eql_invar();
                        if (verbose == 1) cout << endl;
                        if (kinvariant > 0) {
                            if (verbose == 1) {
                                cout << "At the start of evaporation, the activity of water may be constrained by: " << endl;
                                for (int k=1; k<= kinvariant; k++) {
                                    cout << mineral_S[kinvar[k]] << endl;
                                }
                                cout << endl;
                            }
                        } else if (kinvariant == -1) {
                            if (verbose == 1) {
                                cout << "System in thermodynamic desequilibrium" << endl;
                                cout << "The activity of water is constrained at different values" << endl;
                                cout << "by more than one mineral assemblage" << endl;
                            }
                        } else if (kinvariant == -2) {
                            if (verbose == 1) {
                                cout << "System in thermodynamic desequilibrium: inconsisten mineral assemblage" << endl;
                            }
                        } else if (kinvariant == 0) {
                            if (verbose == 1) {
                                cout << "No invariant paragenesis detected" << endl;
                            }
                        }

                        if (kinvariant != 0) {
                            if (verbose == 1) {
                                cout << "If the evaporation program does not start" << endl;
                                cout << "dilute the solution slightly" << endl;
                            }
                        }
                    }
                }

                if (verbose == 1) {
                    cout << endl;
                    cout << "Dilute the solution:" << endl;
                    cout << "  dilution = 1/" << dilute << endl;
                }

                if (dilute > 1) {
                    dil = dilute; dilute = 1;
                    diltot = diltot * dil;

                    pco2_S = ""; pk = pk0; dph = 0.2;
                    for (int i=1; i<=12; i++) {
                        if (i != 11) tot[i] = tot[i] / dil;
                    }
                    for (int i=1; i<=5; i++) {
                        cat[i] = tot[i];
                        nchcat[i] = nch[i];
                    }
                    for (int j=1; j<=3; j++) {
                        ani[j] = tot[j+5];
                        nchani[j] = -nch[j+5];
                    }
                    ani[4] = molal[12] / dil;
                    ani[5] = (molal[14] + molal[17] + molal[17]) / dil;
                    nchani[4] = -nch[12];
                    nchani[5] = -nch[14];
                    unit_S = "molal";
                    eql_density();
                    dil = 1;
                    goto LOOP500;
                }

                if (diltot > 1) {
                    if (verbose == 1) {
                        cout << "The initial solution has been diluted " << diltot << " times" << endl;
                    }
                }
                if (verbose == 1) cout << endl;

                nmneuf = 0;

                if (add_min.size() > 0 or rem_min.size() > 0) {
                    // add_mineral
                    for (int i=0; i<add_min.size(); i++) {
                        string m_S = add_min[i];
                        for (int k=1; k<=nm; k++) {
                            if (mineral_S[k] == m_S) {
                                if (nwmin[k] == 0) {
                                    nmneuf += 1;
                                    nwmin[k] = 1;
                                }
                                if (verbose == 1) {
                                    cout << add_min[i] << " added" << endl;
                                }
                                break;
                            } else if (k == nm and mineral_S[k] != m_S) {
                                if (verbose == 1) {
                                    cout << add_min[i] << " not found" << endl;
                                }
                            }
                        }
                    }

                    for (int i=0; i<rem_min.size(); i++) {
                        string m_S = rem_min[i];
                        for (int k=1; k<=nm; k++) {
                            if (mineral_S[k] == m_S) {
                                if (nwmin[k] == 1) {
                                    nmneuf -= 1;
                                    nwmin[k] = 0;
                                }
                                if (verbose == 1) {
                                    cout << rem_min[i] << " removed" << endl;
                                }
                                break;
                            } else if (k == nm and mineral_S[k] != m_S) {
                                if (verbose == 1) {
                                    cout << rem_min[i] << " not found" << endl;
                                }
                            }
                        }
                    }

                    // Create murtf0 and transfer data LINE 955
                    // clear add_min and rem_min
                    ifstream murtf2;
                    ofstream murtf0;
                    murtf2.open("murtf2");
                    murtf0.open("murtf0", ios::trunc);
                    {
                        getline(murtf2, line); 
                        string value; stringstream linestream(line);

                        getline(linestream, value); nc = stoi(value);
                        getline(linestream, value); na = stoi(value);

                        murtf0 << nc << "," << na << "," << nm0+nmneuf << endl;
                    }

                    for (int i=1; i<=nc+na+1; i++) {
                        getline(murtf2, line);
                        murtf0 << line << endl;
                    }

                    for (int k=1; k<=nm; k++) {
                        if (nwmin[k] == 1) {
                            getline(murtf2, line);
                            murtf0 << line << endl;
                        } else if (nwmin[k] == 0) {
                            getline(murtf2, line);
                        }
                    }
                    murtf2.close();
                    murtf0.close();

                    min_S = "murtf0";

                    add_min.clear();
                    rem_min.clear();
                    goto LOOP400;
                }

                // Transfer file
                {
                    constituant_S = "label,fc,ah2o,sal,eva,ds,ph,alk";
                    for (int i=1; i<=8; i++) {
                        if (tot[i] > 0) {
                            zone_S = upper_to_lower(aq_S[i]);
                            constituant_S += "," + zone_S;
                        }
                    }
                    if (tot[9] != 0) constituant_S += ",b";
                    if (tot[10] != 0) constituant_S += ",si";
                    constituant_S += ",tds";
                    
                    outfile.open("stockage", ios::trunc);
                    outfile << tra_file;
                    outfile.close();

                    outfile.open(tra_file, ios::trunc);
                    outfile << temp << endl;
                    outfile << tinit << endl;
                    outfile << ph << endl;
                    outfile << phinit << endl;
                    outfile << po << endl;
                    outfile << poinit << endl;
                    outfile << diltot << endl;
                    outfile << constituant_S << endl;
                    for (int i=1; i<=10; i++) {
                        outfile << tot[i] << endl;
                    }
                    outfile << molal[15] << endl;
                    for (int i=1; i<=10; i++) {
                        outfile << molal[i] << endl;
                    }
                    outfile << mh2o << endl;
                    outfile << molal[13] << endl;
                    outfile << molal[11] << endl;
                    outfile << molal[14] << endl;
                    outfile << molal[12] << endl;
                    for (int i=16; i<=25; i++) {
                        outfile << molal[i] << endl;
                    }
                    outfile << syst_S << endl;
                    outfile << xi << endl;
                    outfile << output_step << endl;
                    outfile << output << endl;
                    outfile << print_step << endl;
                    outfile << verbose << endl;
                    outfile << output_units << endl;
                    outfile << log_file << endl;
                    outfile << chem_file << endl;
                    outfile << event_file << endl;
                    outfile << min_file << endl;
                    outfile << min_S << endl;
                    outfile << max_sal << endl;
                    outfile << pkmol << endl;
                    outfile << pkeq << endl;

                    outfile.close();
                }

                int result = system("./evp");


            STOP: 
                stop_simulation();
                
        }


        void read_inputs(string fname) {
            // Declare temporary input chemistry variables
            double na, k, li, ca, mg, cl, so4, no3, b, si, alk;

            // Set default values to 0
            na = 0; k = 0; li = 0; ca = 0; mg = 0; cl = 0; 
            so4 = 0; no3 = 0; b = 0; si = 0; alk = 0;

            // Set default parameters
            label = ""; temp = 25; dens = 1.0;; ph = 7; pco2 = ""; max_sal = 0;
            pkmol = 0.001; pkeq = 0.0000000000001; dilute = 1.0; xi = 0;
            print_step = 1; output_step = 1; verbose = 1; output = 1;

            // Open the input file into the file stream
            file.open(fname);

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
                    pco2 = value;
                } else if (name == "system") {
                    getline(linestream, value, ',');
                    syst_S = value;
                } else if (name == "units") {
                    getline(linestream, value, ',');
                    output_units = value;
                } else if (name == "dil") {
                    getline(linestream, value, ',');
                    dilute = stod(value);
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
                } else if (name == "increment") {
                    getline(linestream, value, ',');
                    xi = stod(value) / 100;
                }
            }
            file.close();

            // Convert input chemistry from mmoles to moles, pH to [H+]
            tot[1] = na / 1000; tot[2] = k / 1000; tot[3] = li / 1000;
            tot[4] = ca / 1000; tot[5] = mg / 1000; tot[6] = cl / 1000;
            tot[7] = so4 / 1000; tot[8] = no3 / 1000; tot[9] = b / 1000;
            tot[10] = si / 1000; tot[11] = pow(10, -ph); tot[12] = alk / 1000;
        }

        void eql_actp() {
            double c[10]; double a[12]; double h[4];
            double at, bt, ct, dt, et;
            double u, z, w, v, f, fj, co, s;
            int nc, na, nn;

            c[1] = molal[1]; c[2] = molal[2]; c[3] = molal[3];
            c[4] = molal[4]; c[5] = molal[22]; c[6] = molal[5];
            c[7] = molal[18]; c[8] = molal[23]; c[9] = molal[11];
            a[1] = molal[6]; a[2] = molal[7]; a[3] = molal[25];
            a[4] = molal[12]; a[5] = molal[14]; a[6] = molal[13];
            a[7] = molal[24]; a[8] = molal[19]; a[9] = molal[20];
            a[10] = molal[21]; a[11] = molal[8];
            h[1] = molal[10]; h[2] = molal[9]; h[3] = molal[15];

            {
            file.open("coefft4");
            
            getline(file, line); stringstream linestream(line); string value;
            getline(linestream, value, ','); nc = stoi(value);
            getline(linestream, value, ','); na = stoi(value);
            getline(linestream, value, ','); nn = stoi(value);
            
            if (ndepact == 0) {
                // initialize 1d, 2d, and 3d vectors to 0
                for (int i=0; i<=nc; i++) nzc.push_back(0);
                for (int i=0; i<=na; i++) nza.push_back(0);
                for (int i=0; i<=nc; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=na; k++) temp1d.push_back(0); b0.push_back(temp1d);
                }
                for (int i=0; i<=nc; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=na; k++) temp1d.push_back(0); b1.push_back(temp1d);
                }
                for (int i=0; i<=nc; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=na; k++) temp1d.push_back(0); b2.push_back(temp1d);
                }
                for (int i=0; i<=nc; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=na; k++) temp1d.push_back(0); c0.push_back(temp1d);
                }
                for (int i=0; i<=nc; i++) {
                    vector<vector<double>> temp2d;
                    for (int k=0; k<=nc; k++) {
                        vector<double> temp1d;
                        for (int j=0; j<=na; j++) temp1d.push_back(0); temp2d.push_back(temp1d);
                    }
                    scp.push_back(temp2d);
                }
                for (int i=0; i<=na; i++) {
                    vector<vector<double>> temp2d;
                    for (int k=0; k<=na; k++) {
                        vector<double> temp1d;
                        for (int j=0; j<=nc; j++) temp1d.push_back(0); temp2d.push_back(temp1d);
                    }
                    sap.push_back(temp2d);
                }
                for (int i=0; i<=nc; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=nc; k++) temp1d.push_back(0); tcp.push_back(temp1d);
                }
                for (int i=0; i<=na; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=na; k++) temp1d.push_back(0); tap.push_back(temp1d);
                }
                for (int i=0; i<=nn; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=nc; k++) temp1d.push_back(0); lcp.push_back(temp1d);
                }
                for (int i=0; i<=nn; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=na; k++) temp1d.push_back(0); lap.push_back(temp1d);
                }
                for (int i=0; i<=nn; i++) {
                    vector<vector<double>> temp2d;
                    for (int k=0; k<=nc; k++) {
                        vector<double> temp1d;
                        for (int j=0; j<=na; j++) temp1d.push_back(0); temp2d.push_back(temp1d);
                    }
                    xip.push_back(temp2d);
                }

                // Read data into matrices/vectors
                for (int i=1; i<=nc; i++) {
                    getline(file, line); stringstream linestream(line); string value;

                    getline(linestream, value, ','); string x_S = trim(value);
                    getline(linestream, value, ','); nzc[i] = stoi(value);
                }
                for (int i=1; i<=na; i++) {
                    getline(file, line); stringstream linestream(line); string value;

                    getline(linestream, value, ','); string x_S = trim(value);
                    getline(linestream, value, ','); nza[i] = stoi(value);
                }
                {
                    getline(file, line); stringstream linestream(line); string value;

                    getline(linestream, value, ','); string x_S = trim(value);
                    getline(linestream, value, ','); at = stod(value);
                    getline(linestream, value, ','); bt = stod(value);
                    getline(linestream, value, ','); ct = stod(value);
                    getline(linestream, value, ','); dt = stod(value);
                    getline(linestream, value, ','); et = stod(value);

                    ap0 = temperature(at, bt, ct, dt, et, temp);
                }

                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        {
                            getline(file, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = trim(value);
                            getline(linestream, value, ','); at = stod(value);
                            getline(linestream, value, ','); bt = stod(value);
                            getline(linestream, value, ','); ct = stod(value);
                            getline(linestream, value, ','); dt = stod(value);
                            getline(linestream, value, ','); et = stod(value);

                            b0[i][j] = temperature(at, bt, ct, dt, et, temp);
                        }
                        {
                            getline(file, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = trim(value);
                            getline(linestream, value, ','); at = stod(value);
                            getline(linestream, value, ','); bt = stod(value);
                            getline(linestream, value, ','); ct = stod(value);
                            getline(linestream, value, ','); dt = stod(value);
                            getline(linestream, value, ','); et = stod(value);

                            b1[i][j] = temperature(at, bt, ct, dt, et, temp);
                        }
                        {
                            getline(file, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = trim(value);
                            getline(linestream, value, ','); at = stod(value);
                            getline(linestream, value, ','); bt = stod(value);
                            getline(linestream, value, ','); ct = stod(value);
                            getline(linestream, value, ','); dt = stod(value);
                            getline(linestream, value, ','); et = stod(value);

                            b2[i][j] = temperature(at, bt, ct, dt, et, temp);
                        }
                        {
                            getline(file, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = trim(value);
                            getline(linestream, value, ','); at = stod(value);
                            getline(linestream, value, ','); bt = stod(value);
                            getline(linestream, value, ','); ct = stod(value);
                            getline(linestream, value, ','); dt = stod(value);
                            getline(linestream, value, ','); et = stod(value);

                            c0[i][j] = temperature(at, bt, ct, dt, et, temp);
                        }
                    }
                }
                for (int i=1; i<=nc-1; i++) {
                    for (int j=i+1; j<=nc; j++) {
                        getline(file, line); stringstream linestream(line); string value;

                        getline(linestream, value, ','); string x_S = trim(value);
                        getline(linestream, value, ','); at = stod(value);
                        getline(linestream, value, ','); bt = stod(value);
                        getline(linestream, value, ','); ct = stod(value);
                        getline(linestream, value, ','); dt = stod(value);
                        getline(linestream, value, ','); et = stod(value);

                        tcp[i][j] = temperature(at, bt, ct, dt, et, temp);
                        tcp[j][i] = tcp[i][j];
                    }
                }

                for (int i=1; i<=na-1; i++) {
                    for (int j=i+1; j<=na; j++) {
                        getline(file, line); stringstream linestream(line); string value;

                        getline(linestream, value, ','); string x_S = trim(value);
                        getline(linestream, value, ','); at = stod(value);
                        getline(linestream, value, ','); bt = stod(value);
                        getline(linestream, value, ','); ct = stod(value);
                        getline(linestream, value, ','); dt = stod(value);
                        getline(linestream, value, ','); et = stod(value);

                        tap[i][j] = temperature(at, bt, ct, dt, et, temp);
                        tap[j][i] = tap[i][j];
                    }
                }

                for (int k=1; k<=nc-1; k++) {
                    for (int i=k+1; i<=nc; i++) {
                        for (int j=1; j<=na; j++) {
                            getline(file, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = trim(value);
                            getline(linestream, value, ','); at = stod(value);
                            getline(linestream, value, ','); bt = stod(value);
                            getline(linestream, value, ','); ct = stod(value);
                            getline(linestream, value, ','); dt = stod(value);
                            getline(linestream, value, ','); et = stod(value);

                            scp[k][i][j] = temperature(at, bt, ct, dt, et, temp);
                            scp[i][k][j] = scp[k][i][j];
                        }
                    }
                }

                for (int k=1; k<=na-1; k++) {
                    for (int i=k+1; i<=na; i++) {
                        for (int j=1; j<=nc; j++) {
                            getline(file, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = trim(value);
                            getline(linestream, value, ','); at = stod(value);
                            getline(linestream, value, ','); bt = stod(value);
                            getline(linestream, value, ','); ct = stod(value);
                            getline(linestream, value, ','); dt = stod(value);
                            getline(linestream, value, ','); et = stod(value);

                            sap[k][i][j] = temperature(at, bt, ct, dt, et, temp);
                            sap[i][k][j] = sap[k][i][j];
                        }
                    }
                }

                for (int i=1; i<=nn; i++) {
                    for (int j=1; j<=nc; j++) {
                        getline(file, line); stringstream linestream(line); string value;

                        getline(linestream, value, ','); string x_S = trim(value);
                        getline(linestream, value, ','); at = stod(value);
                        getline(linestream, value, ','); bt = stod(value);
                        getline(linestream, value, ','); ct = stod(value);
                        getline(linestream, value, ','); dt = stod(value);
                        getline(linestream, value, ','); et = stod(value);

                        lcp[i][j] = temperature(at, bt, ct, dt, et, temp);
                    }
                }

                for (int i=1; i<=nn; i++) {
                    for (int j=1; j<=na; j++) {
                        getline(file, line); stringstream linestream(line); string value;

                        getline(linestream, value, ','); string x_S = trim(value);
                        getline(linestream, value, ','); at = stod(value);
                        getline(linestream, value, ','); bt = stod(value);
                        getline(linestream, value, ','); ct = stod(value);
                        getline(linestream, value, ','); dt = stod(value);
                        getline(linestream, value, ','); et = stod(value);

                        lap[i][j] = temperature(at, bt, ct, dt, et, temp);
                    }
                }

                for (int k=1; k<=nn; k++) {
                    for (int i=1; i<=nc; i++) {
                        for (int j=1; j<=na; j++) {
                            xip[k][i][j] = 0;
                        }
                    }
                }
                xip[2][9][1] = -0.0102e0;
                xip[2][1][2] = 0.046e0;
            }
            
            file.close();
            }

            // declare and initialize temporary matrices to 0
            double ec[(nc+1)][(nc+1)] = { 0 }; double ea[(na+1)][(na+1)] = { 0 };
            double fc[(nc+1)][(nc+1)] = { 0 }; double fa[(na+1)][(na+1)] = { 0 };
            double xc[(nc+1)][(nc+1)] = { 0 }; double xa[(na+1)][(na+1)] = { 0 };
            double pp[(nc+1)][(nc+1)] = { 0 }; double qp[(na+1)][(na+1)] = { 0 };
            double p[(nc+1)][(nc+1)] = { 0 }; double q[(na+1)][(na+1)] = { 0 };
            double pf[(nc+1)][(nc+1)] = { 0 }; double qf[(na+1)][(na+1)] = { 0 };
            double cc[(nc+1)][(na+1)] = { 0 }; double bf[(nc+1)][(na+1)] = { 0 };
            double b[(nc+1)][(na+1)] = { 0 }; double bp[(nc+1)][(na+1)] = { 0 };
            double gc[(nc+1)] = { 0 }; double ga[(na+1)] = { 0 }; double gn[(nn+1)] = { 0 };

            bp0 = 1.2e0; mh2o = 55.51e0;

            u = 0; z = 0;
            for (int i=1; i<=nc; i++) {
                u = u + c[i] * pow(nzc[i], 2); z = z + c[i] * nzc[i];
            }
            for (int j=1; j<=na; j++) {
                u = u + a[j] * pow(nza[j], 2); z = z + a[j] * nza[j];
            }
            fi = u / 2; fj = sqrt(fi);
            u = 6 * ap0 * fj;
            for (int i=1; i<=nc-1; i++) {
                for (int j=i+1; j<=nc; j++) {
                    if (nzc[i] == nzc[j]) {
                        ec[i][j] = 0; fc[i][j] = 0;
                    } else {
                        xc[i][j] = 2 * u;
                        xc[i][i] = pow(nzc[i], 2) * u;
                        xc[j][j] = pow(nzc[j], 2) * u;
                        ec[i][j] = (j0(xc[i][j]) - j0(xc[i][i]) / 2 - j0(xc[j][j]) / 2) / fi / 2;
                        fc[i][j] = ((xc[i][j] * j1(xc[i][j]) - xc[i][i] * j1(xc[i][i]) / 2 - 
                                    xc[j][j] * j1(xc[j][j]) / 2) / pow(fi, 2) / 4 - ec[i][j] / fi);
                        ec[j][i] = ec[i][j]; 
                        fc[j][i] = fc[i][j];
                    }
                }
            }
            for (int i=1; i<=na-1; i++) {
                for (int j=i+1; j<=na; j++) {
                    if (nza[i] == nza[j]) {
                        ea[i][j] = 0; fa[i][j] = 0;
                    } else {
                        xa[i][j] = 2 * u; 
                        xa[i][i] = pow(nza[i], 2) * u;
                        xa[j][j] = pow(nza[j], 2) * u;
                        ea[i][j] = (j0(xa[i][j]) - j0(xa[i][i]) / 2 - j0(xa[j][j]) / 2) / fi / 2;
                        fa[i][j] = ((xa[i][j] * j1(xa[i][j]) - xa[i][i] * j1(xa[i][i]) / 
                                        2 - xa[j][j] * j1(xa[j][j]) / 2) / pow(fi, 2) / 4 - ea[i][j] / fi);
                        ea[j][i] = ea[i][j];
                        fa[j][i] = fa[i][j];
                    }
                }
            }
            for (int i=1; i<=nc-1; i++) {
                for (int j=i+1; j<=nc; j++) {
                    pp[i][j] = fc[i][j]; p[i][j] = tcp[i][j] + ec[i][j]; pf[i][j] = p[i][j] + pp[i][j] * fi;
                    pp[j][i] = pp[i][j]; p[j][i] = p[i][j]; pf[j][i] = pf[i][j];
                }
            }
            for (int i=1; i<=na-1; i++) {
                for (int j=i+1; j<=na; j++) {
                    qp[i][j] = fa[i][j]; q[i][j] = tap[i][j] + ea[i][j]; qf[i][j] = q[i][j] + qp[i][j] * fi;
                    qp[j][i] = qp[i][j]; q[j][i] = q[i][j]; qf[j][i] = qf[i][j];
                }
            }
            w = fj * 12;
            for (int i=1; i<=nc; i++) {
                for (int j=1; j<=na; j++) {
                    cc[i][j] = c0[i][j] / sqrt(nzc[i] * nza[j]) / 2;
                    if (nzc[i] == 2 and nza[j] == 2) v = fj * 1.4e0;
                    if (nzc[i] == 1 or nza[j] == 1) v = fj * 2;
                    bf[i][j] = (b0[i][j] + b1[i][j] * exp(-v) + b2[i][j] * exp(-w));
                    b[i][j] = (b0[i][j] + b1[i][j] * (g0(v)) + b2[i][j] * (g0(w)));
                    bp[i][j] = (b1[i][j] * (g1(v)) / fi + b2[i][j] * (g1(w)) / fi);
                }
            }
            f = -ap0 * (fj / (1 + bp0 * fj) + 2 / bp0 * log(1 + bp0 * fj));
            for (int i=1; i<=nc; i++) {
                for (int j=1; j<=na; j++) {
                    f += c[i] * a[j] * bp[i][j];
                }
            }
            for (int i=1; i<=(nc-1); i++) {
                for (int j=(i+1); j<=nc; j++) {
                    f += c[i] * c[j] * pp[i][j];
                }
            }
            for (int i=1; i<=(na-1); i++) {
                for (int j=(i+1); j<=na; j++) {
                    f += a[i] * a[j] * qp[i][j];
                }
            }
            for (int ii=1; ii<=nc; ii++) {
                u = pow(nzc[ii], 2) * f;
                for (int j=1; j<=na; j++) {
                    u += a[j] * (b[ii][j] * 2 + z * cc[ii][j]);
                } 
                for (int i=1; i<=nc; i++) {
                    if (i != ii) {
                        v = 0;
                        for (int j=1; j<=na; j++) {
                            v += a[j] * scp[ii][i][j];
                        }
                        u += c[i] * (p[ii][i] * 2 + v);
                    }
                }
                for (int i=1; i<=(na-1); i++) {
                    for (int j=(i+1); j<=na; j++) {
                        u += a[i] * a[j] * sap[i][j][ii];
                    }
                }
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += c[i] * a[j] * cc[i][j] * nzc[ii];
                    }
                }
                for (int i=1; i<=nn; i++) {
                    u += h[i] * lcp[i][ii] * 2;
                }
                for (int k=1; k<=nn; k++) {
                    for (int j=1; j<=na; j++) {
                        u += h[k] * a[j] * xip[k][ii][j];
                    }
                }
                gc[ii] = exp(u);
            }
            for (int jj=1; jj<=na; jj++) {
                u = pow(nza[jj], 2) * f;
                for (int i=1; i<=nc; i++) u += c[i] * (b[i][jj] * 2 + z * cc[i][jj]);
                for (int i=1; i<=na; i++) {
                    if (i != jj) {
                        v = 0;
                        for (int j=1; j<=nc; j++) v += c[j] * sap[jj][i][j];
                        u += a[i] * (q[jj][i] * 2 + v);
                    }
                }
                for (int i=1; i<=nc-1; i++) {
                    for (int j=i+1; j<=nc; j++) {
                        u += c[i] * c[j] * scp[i][j][jj];
                    }
                }
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += c[i] * a[j] * cc[i][j] * nza[jj];
                    }
                }
                for (int j=1; j<=nn; j++) {
                    u += h[j] * lap[j][jj];
                }
                for (int k=1; k<=nn; k++) {
                    for (int i=1; i<=nc; i++) {
                        u += h[k] * c[i] * xip[k][i][jj];
                    }
                }
                ga[jj] = exp(u);
            }
            for (int k=1; k<=nn; k++) {
                u = 0;
                for (int i=1; i<=nc; i++) {
                    u += c[i] * lcp[k][i] * 2;
                }
                for (int j=1; j<=na; j++) {
                    u += a[j] * lap[k][j] * 2;
                }
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += c[i] * a[j] * xip[k][i][j];
                    }
                }
                gn[k] = exp(u);
            }
            u = -ap0 * pow(fi, 1.5e0) / (1 + bp0 * fj);
            for (int i=1; i<=nc; i++) {
                for (int j=1; j<=na; j++) {
                    u += c[i] * a[j] * (bf[i][j] + z * cc[i][j]);
                }
            }
            for (int i=1; i<=nc-1; i++) {
                for (int j=i+1; j<=nc; j++) {
                    v = 0;
                    for (int k=1; k<=na; k++) {
                        v += a[k] * scp[i][j][k];
                    }
                    u += c[i] * c[j] * (pf[i][j] + v);
                }
            }
            for (int i=1; i<=na-1; i++) {
                for (int j=i+1; j<=na; j++) {
                    v = 0;
                    for (int k=1; k<=nc; k++) v += c[k] * sap[i][j][k];
                    u += a[i] * a[j] * (qf[i][j] + v);
                }
            }
            for (int k=1; k<=nn; k++) {
                for (int i=1; i<=nc; i++) {
                    u += h[k] * c[i] * lcp[k][i];
                }
            }
            for (int k=1; k<=nn; k++) {
                for (int j=1; j<=na; j++) {
                    u += h[k] * a[j] * lap[k][j];
                }
            }
            for (int k=1; k<=nn; k++) {
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += h[k] * c[i] * a[j] * xip[k][i][j];
                    }
                }
            }

            s = 0;
            for (int i=1; i<=nc; i++) s += c[i];
            for (int j=1; j<=na; j++) s += a[j];

            co = 1 + 2 * u / s;
            aw = exp(-s * co / mh2o);
            gact[1] = gc[1]; gact[2] = gc[2]; gact[3] = gc[3];
            gact[4] = gc[4]; gact[22] = gc[5]; gact[5] = gc[6];
            gact[18] = gc[7]; gact[23] = gc[8]; gact[11] = gc[9];
            gact[6] = ga[1]; gact[7] = ga[2]; gact[25] = ga[3];
            gact[12] = ga[4]; gact[14] = ga[5]; gact[13] = ga[6];
            gact[24] = ga[7]; gact[19] = ga[8]; gact[20] = ga[9];
            gact[21] = ga[10]; gact[8] = ga[11];
            gact[10] = aw * aw * pow(gn[1], log(10.e0));
            gact[9] = gn[2]; gact[15] = gn[3];
            gact[16] = 1.; gact[17] = 1.;
            ndepact = 1;
        }

        void eql_density() {
            int nc, na;

            nc = 5; na = 5;

            double s[nc+1][na+1] = {0};
            double ao[nc+1][na+1] = {0}; double bo[nc+1][na+1] = {0};
            double au[nc+1][na+1] = {0}; double bu[nc+1][na+1] = {0};

            file.open("densite");
            for (int i=1; i<=5; i++) {
                for (int j=1; j<=5; j++) {
                    getline(file, line);
                    stringstream linestream(line);
                    string value;

                    getline(linestream, value, ','); x_S = trim(value);
                    getline(linestream, value, ','); ao[i][j] = stod(value);
                    getline(linestream, value, ','); bo[i][j] = stod(value);
                    getline(linestream, value, ','); au[i][j] = stod(value);
                    getline(linestream, value, ','); bu[i][j] = stod(value);
                }
            }
            file.close();

            if (units == "molar") {
                dens = 1; u = 0;
                for (int j=1; j<=na; j++) u += ani[j];
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        s[i][j] = (nchcat[i] + nchani[j]) / 2 * cat[i] * ani[j] / nchcat[i] / nchani[j] / u;
                        dens += ao[i][j] * s[i][j] + bo[i][j] * pow(s[i][j], 2);
                    }
                }
            } else if (units == "molal") {
                dens = 1; u = 0;
                for (int j=1; j<=na; j++) u += ani[j] * nchani[j];
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        s[i][j] = (nchcat[i] + nchani[j]) / 2 * cat[i] * ani[j] / u;
                        dens += au[i][j] * s[i][j] + bu[i][j] * pow(s[i][j], 2);
                    }
                }
            }
        }

        void eql_invar() {
            vector<string> minv_S, minvar_S;
            int n1, n2, n3, n4, ncond, nbmin, ninvar, det, det1, z, ncm, i1, i2;
            double swap, ah2o;
            int ii, jj, kk;

            kinvariant = 0; ncm = 14;
            nbmin = 10;
            ninvar = nbmin + 3;

            for (int i=0; i<=ninvar; i++) minv_S.push_back("");
            for (int i=0; i<=ninvar; i++) minvar_S.push_back("");
            int kinv[ninvar+1] = {0}; double psminv[ninvar+1] = {0};
            double winv[ninvar+1][ncm+1] = {0};
            double psminvar[ninvar+1] = {0}; 
            double t0[ninvar+1][ninvar+1] = {0}; double t1[ninvar+1][ninvar+1] = {0};
            double t2[ninvar+1][ninvar+1] = {0}; double t3[ninvar+1][ncm+1] = {0};
            double t4[ncm+1][ncm+1] = {0}; double tt4[ncm+1] = {0};

            for (int k=1; k<=3; k++) {
                psminv[k] = log10(psc[k]);
            }
            winv[1][11] = 1; winv[1][13] = 1; winv[1][0] = -1;
            winv[2][11] = 1; winv[2][12] = -1; winv[2][14] = 1;
            winv[3][11] = 1; winv[3][12] = 1; winv[3][0] = -1;

            n1 = 3;
            for (int k=1; k<=nm; k++) {
                if (lmin[k] == 1) {
                    n1 += 1;
                    kinv[n1] = k;
                    minv_S[n1] = mineral_S[k];
                    psminv[n1] = log10(psol[k]);
                    for (int j=0; j<=ncm; j++) {
                        winv[n1][j] = wmin[k][j];
                    }
                }
            }
            for (int i=1; i<=n1; i++) {
                swap = winv[i][0];
                winv[i][0] = winv[i][14];
                winv[i][14] = swap;
            }
            for (int i=1; i<=n1; i++) {
                for (int j=i; j<=n1; j++) {
                    t1[i][j] = 0;
                    for (int k=0; k<=ncm-1; k++) {
                        t1[i][j] += winv[i][k] * winv[j][k];
                        t1[j][i] = t1[i][j];
                        t0[i][j] = t1[i][j];
                        t0[j][i] = t0[i][j];
                    }
                }
            }
            for (int k=2; k<=n1; k++) {
                for (int i=k; i<=n1; i++) {
                    if (abs(t1[i][k-1]) > epsilon) {
                        u = t1[k-1][k-1] / t1[i][k-1];
                        for (int j=k; j<=n1; j++) {
                            t1[i][j] = t1[k-1][j] - t1[i][j] * u;
                            if (abs(t1[i][j]) < epsilon) t1[i][j] = 0;
                        }
                    }
                }
            }
            det = 1;
            for (int i=1; i<=n1; i++) {
                if (abs(t1[i][i]) < epsilon) {
                    det = 0;
                    break;
                }
            }

            if (det == 0) {
                n3 = 0;
                n2 = n1-1;
                for (int kk=1; kk<=n1; kk++) {
                    ii = 0;
                    for (int i=1; i<=n1; i++) {
                        if (i != kk) {
                            ii += 1;
                            jj = 0;
                            for (int j=1; j<=n1; j++) {
                                if (j != kk) {
                                    jj += 1;
                                    t2[ii][jj] = t0[i][j];
                                }
                            }
                        }
                    }

                    for (int k=2; k<=n2; k++) {
                        for (int i=k; i<=n2; i++) {
                            if (abs(t2[i][k-1]) > epsilon) {
                                u = t2[k-1][k-1] / t2[i][k-1];
                                for (int j=k; j<=n2; j++) {
                                    t2[i][j] = t2[k-1][j] - t2[i][j] * u;
                                    if (abs(t2[i][j]) < epsilon) t2[i][j] = 0;
                                }
                            }
                        }
                    }

                    det1 = 1;
                    for (int i=1; i<=n2; i++) {
                        if (abs(t2[i][i]) < epsilon) {
                            det1 = 0;
                            break;
                        }
                    }
                    if (det1 == 1) {
                        n3 += 1;
                        kinvar[n3] = kinv[kk];
                        minvar_S[n3] = minv_S[kk];
                        psminvar[n3] = psminv[kk];
                        for (int j=0; j<=ncm; j++) {
                            t3[n3][j] = winv[kk][j];
                        }
                    }
                }

                if (n3 == 0) {
                    kinvariant = -1;
                } else if (n3 > 0) {
                    n4 = ncm;
                    for (int j=ncm; j>=1; --j) {
                        u = 0;
                        for (int i=1; i<=n3; i++) {
                            u += pow(t3[i][j], 2);
                        }
                        if (u < epsilon) {
                            for (int k=j+1; k<=n4; k++) {
                                for (int i=1; i<=n3; i++) {
                                    t3[i][k-1] = t3[i][k];
                                }
                            }
                            n4 -= 1;
                        }
                    }

                    for (int i=1; i<=n4; i++) {
                        for (int j=1; j<=n4; j++) {
                            t4[i][j] = 0;
                            for (int k=1; k<=n3; k++) {
                                t4[i][j] += t3[k][i] * t3[k][j];
                                t4[j][i] = t4[i][j];
                            }
                        }
                    }
                    for (int i=1; i<=n4; i++) {
                        tt4[i] = 0;
                        for (int k=1; k<=n3; k++) {
                            tt4[i] += t3[k][i] * psminvar[k];
                        }
                    }

                    for (int k=2; k<=n4; k++) {
                        for (int i=k; i<=n4; i++) {
                            if (abs(t4[i][k-1]) > epsilon) {
                                u = t4[k-1][k-1] / t4[i][k-1];
                                for (int j=k; j<=n4; j++) {
                                    t4[i][j] = t4[k-1][j] - t4[i][j] * u;
                                    if (abs(t4[i][j]) < epsilon) {
                                        t4[i][j] = 0;
                                    }
                                }
                                tt4[i] = tt4[k-1] - tt4[i] * u;
                            }
                        }
                    }

                    if (abs(t4[n4][n4]) > epsilon) {
                        ah2o = pow(10, (tt4[n4] / t4[n4][n4]));
                        if (ah2o > 1 or ah2o <= 0.01) {
                            kinvariant = -2;
                        } else {
                            kinvariant = n3;
                            for (int i=kinvariant; i>=1; --i) {
                                if (kinvar[i] == 0) {
                                    for (int k=1; k<=kinvariant-1; k++) {
                                        kinvar[k] = kinvar[k+1];
                                    }
                                    kinvariant -= 1;
                                }
                            }
                        }
                    } else if(abs(t4[n4][n4]) <= epsilon) {
                        kinvariant = -2;
                    }
                }
            }
        }

        void stop_simulation() {

            exit(0);
        }
};

int main() {
    Simulation test("input.dat");

    return 0;
}
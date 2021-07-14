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
const int max_conv = 100; const double epsilon = 1e-8;

class Simulation {
    public:

        // declare variable
        double fc, temp, tinit, ph, phinit, po, poinit, diltot, xi, stdmax, pkmol, pkeq, xi0;
        double at, bt, ct, dt, et, u, sc, cmax, sa, amax, dca, delta, ctot0;
        double xinv, gmin, m, p, g, psc3, fi, fi0, dens, std, ee, hco3, co3, ctot, s, pco, mwev, mwev0;
        double mh2o, alk, ap0, bp0, aw;
        int nc, na, ncm, nm, nminer, nminer0, nbmin, ndepact;
        int kinvariant, ncpt, ncomp, ncmpt, nw, nu;
        int npasi, npasf, kneuf, ii, jj, l, ncmptinv, ninv;
        int output_step, output, print_step, verbose;
        int ic, ia, ix, iy, initdeseq, ki, ksupprim, nperitec;
        string mineraux_S, constituant_S, inc_S, q1_S, q2_S;
        string a_S, x_S, y_S, zone_S, m0_S, q_S, q0_S, my_S, my0_S;
        string tra_file, log_file, chem_file, event_file, min_file;
        string system, units, miner_S;
        string line, value;
        ifstream infile; ofstream outfile;

        // declare arrays and vectors
        array<string,n+1> aq_S;
        double tot[ntot+1] = {0}; double totinit[ntot+1] = {0};
        double tot0[ntot+1] = {0}; double totest[ntot+1] = {0};
        double psc[ncomplex+1]; double gact[n+1] = {0}; 
        double gact0[n+1] = {0}; double gact1[n+1] = {0}; 
        double mol[n+1] = {0}; double mol0[n+1] = {0}; 
        double mol1[n+1] = {0}; double molal[n+1] = {0}; 
        double molal0[n+1] = {0}; double act[n+1] = {0}; 
        double act0[n+1] = {0};  double atom[n+1] = {0}; 
        double kmat[n+1][n+1] = {0}; int nch[n+1] = {0};


        // Declare vectors
        vector<int> ica, kinvar, linvar;
        vector<int> lmin, lmin0, lmin1;
        vector<vector<double>> wmin;
        vector<double> mu, mum, psol, psol0, pai, pai0;
        vector<double> min, min0, minp, minp0, c_ion;
        vector<string> nom_ion_S, mineral_S;

        // Declare vectors for density calculation
        vector<vector<double>> au, bu;

        // Declare vectors for activity calculation
        vector<int> nzc, nza;
        vector<vector<double>> b0, b1, b2, c0;
        vector<vector<double>>  tcp, tap, lcp, lap;
        vector<vector<vector<double>>> scp, sap, xip;

        Simulation() {
            run_evp();
        }

    private:

        void run_evp() {

            // set variables
            ncpt = 0; mwev = 0; fc = 1; q0_S = "";
            mh2o = 55.51; ndepact = 0; ksupprim = 0;
            kinvariant = 0;
            npasi = 0; npasf = 0;
            nperitec = 0; initdeseq = 0; ninv = 0;
            xinv = 0;

            read_transfer();

            if (verbose == 1) {
                time_t now;
                time (&now);
                cout << "\nThis is EVP.............." << endl;;
                cout << endl;
                cout << "STARTING THE EVAPORATION PROGRAM" << endl << endl;
                cout << asctime(localtime(&now)) << endl << endl;
            }

            // read chemical species thermodynamic data
            infile.open("resources/aquv.dat");
            for (int i=0; i<=n; i++) {
                getline(infile, line); stringstream linestream(line);
                getline(linestream, value, ','); aq_S[i] = trim(value);
                getline(linestream, value, ','); atom[i] = stod(value);
                getline(linestream, value, ','); nch[i] = stoi(value);   
            }
            infile.close();

            // read coefficients of partial derivatives of EVP equation set
            infile.open("resources/matrice2");
            for (int i=1; i<=n; i++) {
                getline(infile, line); stringstream linestream(line);
                for (int j=1; j<=n; j++) {
                    getline(linestream, value, ','); kmat[i][j] = stod(value);
                }
            }
            infile.close();

            // set the increment mode based on user input
            if (xi == 0) {
                inc_S = "auto";
            } else {
                inc_S = "manu";
            }
            xi0 = xi;

            // set ica to 1 for species 1-25
            for (int i=1; i<=n; i++) ica[i] = 1;

            if (output == 1) {
                // add heading to the chemistry output file
                outfile.open(chem_file, ios::trunc);
                    outfile << constituant_S << endl;
                outfile.close();

                // add initial events to event file
                outfile.open(event_file, ios::trunc);
                    outfile << "Temperature of solution = " << tinit << " Deg C.   Temperature of simulation = " << temp << " Deg C." << endl;
                    if (diltot > 1) {
                        outfile << "The initial solution has been diluted " << diltot << " times" << endl;
                    }
                    if (ph != phinit) {
                        outfile << "Initial Log(pco2) = " << poinit << "    Selected Log(pco2 = " << po << endl;
                        outfile << "Initial pH = " << phinit << "    Calculated pH = " << ph << endl;
                    }
                outfile.close();

                // open and clear the mineral output file
                outfile.open(min_file, ios::trunc);
                outfile.close();
            }

            infile.open("resources/complex3");
            for (int i=1; i<=ncomplex; i++) {
                getline(infile, line);
                stringstream linestream(line);

                getline(linestream, value, ','); 
                getline(linestream, value, ','); at = stod(value);
                getline(linestream, value, ','); bt = stod(value);
                getline(linestream, value, ','); ct = stod(value);
                getline(linestream, value, ','); dt = stod(value);
                getline(linestream, value, ','); et = stod(value);

                psc[i] = pow(10, (at + bt / 300 * temp + ct / 30000 * pow(temp, 2) + dt / 3000000 * pow(temp, 3) + et / 300000000 * pow(temp, 4)));
            }
            infile.close();
            psc3 = psc[3];
            psc[3] = psc[3] * psc[14] * pow(10, po);
            
            // Read mineral data base
            infile.open(miner_S);
            {
                getline(infile, line);
                stringstream linestream(line);
                getline(linestream, value, ','); nc = stoi(value);
                getline(linestream, value, ','); na = stoi(value);
                getline(linestream, value, ','); nm = stoi(value);

                ncm = nc + na + 1;
                
                // allocate vectors and arrays
                for (int i=0; i<=nm; i++) {
                        vector<double> temp1d;
                        for (int k=0; k<=ncm; k++) temp1d.push_back(0);
                        wmin.push_back(temp1d);
                }
                for (int i=0; i<=ncm; i++) {
                    mu.push_back(0); 
                    nom_ion_S.push_back("");
                    c_ion.push_back(0);
                }
                for (int i=0; i<=nm; i++) {
                    mineral_S.push_back(""); linvar.push_back(0);
                    mum.push_back(0); psol.push_back(0);
                    psol0.push_back(0); pai.push_back(0);
                    pai0.push_back(0); lmin.push_back(0);
                    lmin0.push_back(0); lmin1.push_back(0);
                    min.push_back(0); min0.push_back(0);
                    minp.push_back(0); minp0.push_back(0);
                }
                for (int i=1; i<=ncm; i++) {
                    getline(infile, line); stringstream linestream(line);
                    getline(linestream, value, ','); a_S = upper_to_lower(trim(value));
                    getline(linestream, value, ','); at = stod(value);
                    getline(linestream, value, ','); bt = stod(value);
                    getline(linestream, value, ','); ct = stod(value);
                    getline(linestream, value, ','); dt = stod(value);
                    getline(linestream, value, ','); et = stod(value);

                    for (int j=1; j<=ncm; j++) {
                        if (a_S == aq_S[j]) {
                            mu[j] = at + bt / 300 * temp + ct / 30000 * pow(temp, 2) + dt / 3000000 * pow(temp, 3) + et / 300000000 * pow(temp, 4);
                        }
                    }
                }
                for (int k=1; k<=nm; k++) {
                    getline(infile, line); stringstream linestream(line);
                    getline(linestream, value, ','); mineral_S[k] = trim(value);
                    getline(linestream, value, ','); ncomp = stoi(value);

                    for (int i=1; i<=ncomp; i++) {
                        getline(linestream, value, ','); c_ion[i] = stod(value);
                        getline(linestream, value, ','); nom_ion_S[i] = trim(value);
                    }

                    getline(linestream, value, ','); at = stod(value);
                    getline(linestream, value, ','); bt = stod(value);
                    getline(linestream, value, ','); ct = stod(value);
                    getline(linestream, value, ','); dt = stod(value);
                    getline(linestream, value, ','); et = stod(value);

                    for (int i=1; i<=ncomp; i++) {
                        x_S = upper_to_lower(nom_ion_S[i]);
                        for (int j=0; j<=ncm; j++) {
                            if (x_S == aq_S[j]) {
                                wmin[k][j] = c_ion[i];
                            }
                        }
                    }
                    mum[k] = (at + bt / 300 * temp + ct / 30000 * pow(temp, 2) + dt / 3000000 * pow(temp, 3) + et / 300000000 * pow(temp, 4));
                }
            }
            infile.close();

            for (int k=1; k<=nm; k++) {
                u = mum[k];
                for (int i=0; i<=ncm; i++) {
                    u = u - wmin[k][i] * mu[i];
                }
                psol[k] = exp(u);
                psol0[k] = psol[k];
            }

            if (output == 1) {
                outfile.open(min_file, ios::app);
                outfile << "fc,";
                for (int k=1; k<=nm-1; k++) {
                    zone_S = mineral_S[k] + ",";
                    outfile << zone_S;
                }
                zone_S = mineral_S[nm];
                outfile << zone_S;
                outfile << endl;
                outfile.close();
            }
            sc = 0; cmax = 0;
            for (int i=1; i<=n; i++) {
                if (nch[i] > 0) {
                    sc += mol[i] * nch[i];
                    if (mol[i] * nch[i] > cmax) {
                        cmax = mol[i] * nch[i];
                        ic = i;
                    }
                }
            }
            sa = 0; amax = 0;
            for (int i=1; i<=n; i++) {
                if (nch[i] < 0) {
                    sa += mol[i] * (-nch[i]);
                    if (mol[i] * (-nch[i]) > amax) {
                        amax = mol[i] * (-nch[i]);
                        ia = i;
                    }
                }
            }
            dca = 200 * abs(sc - sa) / (sc + sa);
            delta = sc - sa;
            mol[ic] = mol[ic] - delta / 2 / nch[ic];
            mol[ia] = mol[ia] + delta / 2 / (-nch[ia]);
            sc = 0; sa = 0;
            for (int i=1; i<=n; i++) {
                if (nch[i] > 0) sc += mol[i] * nch[i];
                if (nch[i] < 0) sa += mol[i] * (-nch[i]);
            }
            if (verbose == 1) {
                cout << "Sum of cations = " << sc << endl;
                cout << "Sum of anions  = " << sa << endl;
                cout << endl;
            }

            for (int i=1; i<=11; i++) {
                totinit[i] = 0;
                for (int j=1; j<=n; j++) {
                    totinit[i] += kmat[i][j] * mol[j];
                    tot[i] = totinit[i];
                    tot0[i] = totinit[i];
                }
            }
            tot[12] = 0;
            ctot0 = mol[0] + mol[14] + mol[15] + mol[16] + mol[17];
            
            evp_actp();
            for (int i=0; i<=n; i++) {
                molal[i] = mol[i] * mh2o / mol[11];
                act[i] = molal[i] * gact[i];
            }
            for (int k=1; k<=nm; k++) {
                pai[k] = 1;
                for (int i=1; i<=ncm; i++) {
                    pai[k] = pai[k] * pow(act[i], wmin[k][i]); 
                }
                if (pai[k] >= psol[k]) lmin[k] = 1;
            }

            if (output == 1) {
                outfile.open(event_file, ios::app);
                for (int k=1; k<=nm; k++) {
                    if (lmin[k] == 1) {
                        outfile << "Initial solution oversaturated in " << mineral_S[k] << endl;
                    }
                }
                outfile.close();
            }

            for (int k=1; k<=nm; k++) {
                if (lmin[k] == 1) {
                    psol[k] = pai[k];
                }
            }

            LOOP500:
            while (true) {
                ncmpt = 0;
                ix = 1; iy = 2;
                initdeseq = 0;
                if (mwev == 0) {
                    for (int k=1; k<=nm; k++) {
                        if (psol[k] > psol0[k]) initdeseq = 1;
                    }
                    if (initdeseq == 1) {
                        for (int k=1; k<=nm; k++) {
                            if (psol[k] * 0.95 > psol0[k]) {
                                psol[k] = psol[k] * 0.95;
                                if (verbose == 1) {
                                    cout << setw(20) << mineral_S[k] << setw(15) << psol[k] << setw(15) << psol0[k] << endl;
                                }
                            } else if (psol[k] * 0.95 <= psol0[k]) {
                                psol[k] = psol0[k];
                            }
                        }
                    }
                }

                nw = 1;
                while (nw != 0) {
                    ncmpt += 1;
                    m0_S = "";
                    for (int k=1; k<=nm; k++) {
                        if (lmin[k] == 1) {
                            m0_S += "_" + mineral_S[k];
                        }
                    }
                    for (int i=0; i<=n; i++) gact1[i] = gact[i];
                    evp_actp();
                    if (kinvariant == 0) {
                        for (int i=0; i<=n; i++) gact[i] = (gact[i] + gact1[i] * ix) / iy;
                    }
                    for (int i=0; i<=n; i++) {
                        molal[i] = mol[i] * mh2o / mol[11];
                        act[i] = molal[i] * gact[i];
                    }
                    for (int k=1; k<=nm; k++) {
                        pai[k] = 1;
                        for (int i=1; i<=ncm; i++) {
                            pai[k] = pai[k] * pow(act[i], wmin[k][i]);
                        }
                        if (pai[k] >= psol[k]) {
                            if (min[k] >= 0) {
                                lmin[k] = 1;
                            } else if (min[k] < 0) {
                                lmin[k] = 0;
                                min[k] = 0;
                            }
                        } else if (pai[k] < psol[k]) {
                            if (min[k] <= 0) {
                                lmin[k] = 0;
                                min[k] = 0;
                            } else if (min[k] > 0) {
                                lmin[k] = 1;
                            }
                        }
                    }

                    for (int k=1; k<=nm; k++) {
                        if (psol[k] == 1.e+50) {
                            if (pai[k] < psol0[k] * 0.9) {
                                linvar[k] = 0;
                            } else if (pai[k] >= psol0[k]) {
                                linvar[k] = 1;
                            }
                        }
                    }

                    mineraux_S = "";
                    nminer = 0;
                    for (int k=1; k<=nm; k++) {
                        if (lmin[k] == 1) {
                            nminer += 1;
                            mineraux_S += "_" + mineral_S[k];
                        }
                    }
                    if (ncpt == 1 or ncpt % output_step == 0) {
                        if (verbose == 1) {
                            if (nminer == 0) {
                                cout << ncmpt << " No_minerals" << endl;
                            } else {
                                cout << setw(4) << left << ncmpt << mineraux_S << endl;
                            }
                        }
                    }

                    if (mwev > 0 and fc != 1 and nminer - nminer0 >= 2) {
                        xi = xi / 2;
                        if (xi < epsilon) {
                            if (verbose == 1) {
                                cout << endl;
                                cout << "Program unstable" << endl;
                                cout << "Restart the initialization program (EQL...)" << endl;
                                cout << "and lower the limits of convergence" << endl;
                                cout << endl;
                            }
                            if (output == 1) {
                                outfile.open(event_file, ios::app);
                                outfile << "Program unstable" << endl;
                                outfile << "Restart the initialization program (EQL...)" << endl;
                                outfile << "and lower the limits of convergence" << endl;
                                outfile.close();
                            }
                            stop_simulation();
                        }
                        if (verbose == 1) {
                            cout << "reduction of increment at " << xi << endl;
                        }
                        for (int i=0; i<=n; i++) {
                            mol[i] = mol0[i];
                        }
                        for (int i=1; i<=ntot; i++) {
                            tot[i] = tot0[i];
                        }
                        for (int k=1; k<=nm; k++) {
                            lmin[k] = lmin0[k];
                            min[k] = min0[k];
                            minp[k] = minp0[k];
                        }
                        mwev = mwev0;
                        nminer = nminer0;
                        mwev = mwev + mol[11] * xi;
                        tot[11] = tot[11] - 2 * mol[11] * xi;
                        goto LOOP500;
                    }

                    if (nminer > 1 and mineraux_S != m0_S) {
                        ix = 2; iy = 3;
                        evp_invar();
                        if (kinvariant > 0) {
                            if (system == "o") {
                                for (int i=1; i<=kinvariant; i++) {
                                    if (minp[kinvar[i]] == 0 and min[kinvar[i]] == 0) kneuf = kinvar[i];
                                }
                                goto LOOP2000;
                            } else if (system == "c") {
                                for (int i=0; i<=n; i++) {
                                    mol[i] = mol0[i];
                                    molal[i] = molal0[i];
                                    gact[i] = gact0[i];
                                    act[i] = act0[i];
                                }
                                for (int i=1; i<=ntot; i++) {
                                    tot[i] = tot0[i];
                                }
                                for (int k=1; k<=nm; k++) {
                                    pai[k] = pai0[k];
                                    min[k] = min0[k];
                                }
                                mwev = mwev0;
                                nminer = nminer0;
                                fi = fi0;
                            }
                        }
                    }

                    for (int i=0; i<=n; i++) mol1[i] = mol[i];
                    // LINE 513
                    if (kinvariant == 0) {
                        reseq();
                    } else if (kinvariant > 0) {
                        reseqinv();
                        mwev += xinv/2;
                        tot[11] = tot[11] - xinv;
                    }
                    mol[0] = mol[15] * gact[15] * mol[13] * gact[13] / mol[11] / gact[11] / psc3 / gact[0];
                    nw = 0;
                    for (int i=1; i<=n; i++) {
                        if (mol[i] > 0) {
                            if (200 * abs(mol[i] - mol1[i]) / (mol[1] + mol1[i]) > pkmol) nw = 1;
                        }
                    }
                    ki = kinvariant;
                    if (kinvariant > 0) {
                        for (int k=1; k<=kinvariant; k++) {
                            if (min[kinvar[k]] <= 0) {
                                lmin[kinvar[k]] = 0;
                                min[kinvar[k]] = 0;
                                psol[kinvar[k]] = 1.e+50;
                                mwev = mwev + mol[11] * xi;
                                tot[11] -= 2 * mol[11] * xi;
                                ki = 0; nw = 1;
                            }
                        }
                    }
                    kinvariant = ki;
                    if (nw == 1) {
                        for (int i=0; i<=n; i++) mol[i] = (mol[i] + mol1[i]) / 2;
                    }
                    if (ncmpt == 500) {
                        if (verbose == 1) {
                            cout << endl;
                            cout << "Program unstable" << endl;
                            cout << "Restart the initialization program (EQL...)" << endl;
                            cout << "and lower the limits of convergence." << endl;
                            cout << "Set the increment in manual mode at a value lower than 0.5" << endl;
                            cout << endl;
                        }
                        if (output == 1) {
                            outfile.open(event_file, ios::app);
                            outfile << endl;
                            outfile << "Program unstable" << endl;
                            outfile << "Restart the initialization program (EQL...)" << endl;
                            outfile << "and lower the limits of convergence." << endl;
                            outfile << "Set the increment in manual mode at a value lower than 0.5" << endl;
                            outfile.close();
                        }
                        stop_simulation();
                    }
                }

                for (int k=1; k<=nm; k++) {
                    if (psol[k] == 1.e+50 and linvar[k] == 0) {
                        psol[k] = psol0[k];
                        if (verbose == 1) {
                            cout << "resetting: " << mineral_S[k] << endl;
                        }
                    }
                    linvar[k] = 0;
                }

                npasi += 1;
                npasf += 1;
                ncpt += 1;

                if (system == "o") {
                    for (int k=1; k<=nm; k++) {
                        minp[k] = minp[k] + min[k];
                    }
                    for (int i=1; i<=10; i++) {
                        for (int k=1; k<=nm; k++) {
                            tot[i] = tot[i] - wmin[k][i] * min[k];
                        }
                    }
                    for (int j=1; j<=ncm; j++) {
                        for (int k=1; k<=nm; k++) {
                            tot[11] = tot[11] - wmin[k][j] * kmat[11][j] * min[k];
                        }
                    }
                }
                for (int i=1; i<=ntot-1; i++) {
                    totest[i] = 0;
                    for (int j=1; j<=n; j++) {
                        totest[i] = totest[i] + kmat[i][j] * mol[j];
                    }
                }
                totinit[12] = 0; totest[12] = 0;
                for (int j=1; j<=n; j++) {
                    if (kmat[12][j] > 0) {
                        totest[12] = totest[12] + kmat[12][j] * mol[j];
                    } else if (kmat[12][j] < 0) {
                        totinit[12] = totinit[12] + kmat[12][j] * mol[j];
                    }
                }
                totinit[12] = -totinit[12];
                for (int i=1; i<=10; i++) {
                    for (int k=1; k<=nm; k++) {
                        if (system == "c") totest[i] = totest[i] + min[k] * wmin[k][i];
                        if (system == "o") totest[i] = totest[i] + minp[k] * wmin[k][i];
                    }
                }
                for (int j=1; j<=ncm; j++) {
                    for (int k=1; k<=nm; k++) {
                        if (system == "c") totest[11] = totest[11] + wmin[k][j] * kmat[11][j] * min[k];
                        if (system == "o") totest[11] = totest[11] + wmin[k][j] * kmat[11][j] * minp[k];
                    }
                }
                totest[11] = totest[11] + mwev * 2;
                evp_density();

                fc = mh2o / mol[11];
                alk = molal[12] - molal[13] + molal[14] * 2 + molal[15] + (molal[16] + molal[17]) * 2;
                alk = alk + molal[18] + molal[19] + molal[20] + molal[21] * 2;
                alk = alk + molal[22] + molal[23] + molal[24] - molal[25];
                std = -molal[11] * atom[11];
                for (int i=1; i<=n; i++) {
                    std = std + molal[i] * atom[i];
                }
                ee = 1000000 * dens / (1000 + std);
                hco3 = 0; co3 = 0;
                for (int k=1; k<=nm; k++) {
                    if (system == "c") {
                        hco3 = hco3 + wmin[k][15] * min[k];
                        co3 = co3 + wmin[k][14] * min[k];
                    } else if (system == "o") {
                        hco3 = hco3 + wmin[k][15] * minp[k];
                        co3 = co3 + wmin[k][14] * minp[k];
                    }
                }
                ctot = mol[0] + mol[14] + mol[15] + mol[16] + mol[17] + hco3 + co3;
                my_S = "";
                for (int i=1; i<=nm; i++) {
                    if (lmin[i] == 1 or min[i] != 0) {
                        my_S += "_" + mineral_S[i];
                    }
                }

            LOOP600:
                if (ncpt == 1 or ncpt % output_step == 0) {
                    if (verbose == 1) {
                        cout << endl;
                        q_S = "";
                        if (nminer > 0) {
                            if (system == "c") {
                                cout << setw(24) << " " << setw(42) << "MOLES PREC" << setw(14) << "TESTS" << endl;
                            } else if (system == "o") {
                                cout << setw(24) << " " << setw(14) << "MOLES 1 STEP" << setw(28) << "MOLES TOT" << setw(14) << "TESTS" << endl;
                            }
                        }
                        cout << endl;
                        for (int i=1; i<=nm; i++) {
                            if (lmin[i] == 1 or min[i] != 0) {
                                if (system == "c") {
                                    if (min[i] > min0[i]) {
                                        q_S = "P";
                                    } else if (min[i] < min0[i]) {
                                        q_S = "D";
                                    } else if (min[i] == min0[i]) {
                                        q_S = "=";
                                    }
                                }
                                x_S = mineral_S[i] + " " + q_S;
                                u = 200 * abs(pai[i] - psol[i]) / (pai[i] + psol[i]);

                                if (system == "o") {
                                    cout << setw(24) << x_S << setw(14) << min[i] << setw(28) << minp[i] << setw(14) << u << endl;
                                } else {
                                    cout << setw(24) << x_S << setw(42) << min[i] << setw(14) << u << endl;
                                }
                            }
                        }
                        if (my_S == "") cout << "No_minerals" << endl;
                        cout << endl;

                        cout << setw(10) << " " << setw(14) << "MOLES" << setw(14) << "MOLALITIES" << setw(14) << "ACT COEFF" << setw(14) << "MOLAL TOT" << setw(14) << "TESTS" << endl;
                        cout << endl;

                        for (int i=1; i<=ntot; i++) {
                            if (tot[i] > 0 or i == 12) {
                                u = 200 * abs(totest[i] - totinit[i]) / (totest[i] + totinit[i]);
                                if (i <= 10) {
                                    s = 0;
                                    for (int j=1; j<=n; j++) {
                                        s += molal[j] * kmat[i][j];
                                    }
                                    cout << setw(10) << aq_S[i] << setw(14) << mol[i] << setw(14) << molal[i] << setw(14) << gact[i] << setw(14) << s << setw(14) << u << endl;
                                } else if (i > 10) {
                                    cout << setw(10) << aq_S[i] << setw(14) << mol[i] << setw(14) << molal[i] << setw(14) << gact[i] << setw(14) << " " << setw(14) << u << endl;
                                }
                            }
                        }
                        for (int i=ntot+1; i<=n; i++) {
                            if (mol[i] > 0) {
                                p = 1;
                                for (int j=1; j<=n; j++) {
                                    p = p * pow(act[j], kmat[i][j]);
                                }
                                u = 200 * abs(p - psc[i-12]) / (p + psc[i-12]);
                                cout << setw(10) << aq_S[i] << setw(14) << mol[i] << setw(14) << molal[i] << setw(14) << gact[i] << setw(14) << " " << setw(14) << u << endl;
                            }
                        }
                        pco = log10(act[15] * act[13] / act[11] / psc3 / psc[14]);
                        u = 200 * abs(pco - po) / (pco + po);
                        cout << setw(10) << aq_S[0] << setw(14) << mol[0] << setw(14) << molal[0] << setw(14) << gact[0] << setw(14) << " " << setw(14) << u << endl;
                        cout << endl;
                        cout << setw(40) << left << tra_file << "concentration factor = " << fc << endl;
                        cout << setw(22) << "ionic strength      = " << setw(18) << fi << setw(23) << "salinity (%)         = " << setw(17) << std/(std+1000) * 100 << endl;
                        cout << setw(22) << "activity of water   = " << setw(18) << act[11] << setw(23) << "water evapor. (mol)  = " << setw(17) << mwev << endl;
                        cout << setw(22) << "pH                  = " << setw(18) << -log10(act[13]) << setw(23) << "CO2 exchanged (mol)  = " << setw(17) << ctot - ctot0 << endl;
                        cout << setw(22) << "alkalinity (eq/kg)  = " << setw(18) << alk << setw(23) << "Log PCO2             = " << setw(17) << pco << endl;
                        cout << setw(22) << "number of steps     = " << setw(18) << ncpt << setw(23) << "molal/molar factor   = " << setw(17) << ee/1000 << endl;
                        if (kinvariant == 0) {
                            cout << setw(22) << "increment (%)       = " << setw(18) << xi * 100 << setw(23) << "density              = " << setw(17) << dens << endl;
                        } else {
                            cout << setw(22) << "increment (moles)   = " << setw(18) << xinv << setw(23) << "density              = " << setw(17) << dens << endl;
                        }
                        cout << endl;
                    }
                }

                if (output == 1 and (ncpt == 1 or my_S != my0_S or ncpt % print_step == 0)) {
                    outfile.open(log_file, ios::app);
                    my_S = ""; q_S = "";
                    outfile << endl;
                    if (system == "c") {
                        outfile << setw(24) << " " << setw(42) << "MOLES PREC" << setw(14) << "TESTS" << endl;
                    } else if (system == "o") {
                        outfile << setw(24) << " " << setw(14) << "MOLES 1 STEP" << setw(28) << "MOLES TOT" << setw(14) << "TESTS" << endl;
                    }
                    outfile << endl;

                    for (int i=1; i<=nm; i++) {
                        if (lmin[i] == 1 or min[i] != 0) {
                            if (system == "c") {
                                if (min[i] > min0[i]) {
                                    q_S = "(p)";
                                } else if (min[i] < min0[i]) {
                                    q_S = "(d)";
                                } else if (min[i] == min0[i]) {
                                    q_S = "(=)";
                                }
                            }

                            x_S = mineral_S[i] + " " + q_S;
                            u = 200 * abs(pai[i] - psol[i]) / (pai[i] + psol[i]);

                            if (system == "o") {
                                outfile << setw(24) << x_S << setw(14) << min[i] << setw(28) << minp[i] << setw(14) << u << endl;
                            } else {
                                outfile << setw(24) << x_S << setw(42) << min[i] << setw(14) << u << endl;
                            }
                            my_S += "_" + mineral_S[i];
                        }
                    }
                    if (my_S == "") outfile << "No_minerals" << endl;
                    outfile << endl;

                    if (ncpt == 1) {
                        outfile << setw(10) << " " << setw(14) << "MOLES" << setw(14) << "MOLALITIES" << setw(14) << "ACT COEFF" << setw(14) << "MOLAL TOT" << setw(14) << "TESTS" << endl;
                        outfile << endl;
                    }

                    for (int i=1; i<=ntot; i++) {
                        if (tot[i] > 0 or i==12) {
                            u = 200 * abs(totest[i] - totinit[i]) / (totest[i] + totinit[i]);
                            if (i <= 10) {
                                s = 0;
                                for (int j=1; j<=n; j++) {
                                    s += molal[j] * kmat[i][j];
                                }
                                outfile << setw(10) << aq_S[i] << setw(14) << mol[i] << setw(14) << molal[i] << setw(14) << gact[i] << setw(14) << s << setw(14) << u << endl;
                            } else if (i>10) {
                                outfile << setw(10) << aq_S[i] << setw(14) << mol[i] << setw(14) << molal[i] << setw(14) << gact[i] << setw(14) << " " << setw(14) << u << endl;
                            }
                        }
                    }

                    for (int i=ntot+1; i<=n; i++) {
                        if (mol[i] > 0) {
                            p = 1;
                            for (int j=1; j<=n; j++) {
                                p = p * pow(act[j], kmat[i][j]);
                            }
                            u = 200 * abs(p - psc[i-12]) / (p + psc[i-12]);
                            outfile << setw(10) << aq_S[i] << setw(14) << mol[i] << setw(14) << molal[i] << setw(14) << gact[i] << setw(14) << " " << setw(14) << u << endl;
                        }
                    }
                    pco = log10(act[15] * act[13] / act[11] / psc3 / psc[14]);
                    u = 200 * abs(pco - po) / (pco + po);
                    outfile << setw(10) << aq_S[0] << setw(14) << mol[0] << setw(14) << molal[0] << setw(14) << gact[0] << setw(14) << " " << setw(14) << u << endl;
                    outfile << endl;
                    outfile << tra_file << endl;
                    outfile << "concentration factor     = " << fc << endl;
                    outfile << "ionic strength           = " << fi << endl;
                    outfile << "salinity                 = " << std/(std+1000) << " g/kg" << endl;
                    outfile << "activity of water        = " << act[11] << endl;
                    outfile << "pH                       = " << -log10(act[13]) << endl;
                    outfile << "alkalinity               = " << alk << " eq/kg" << endl;
                    outfile << "water evaporated         = " << mwev << " moles" << endl;
                    outfile << "CO2 exchanged            = " << ctot-ctot0 << " moles" << endl;
                    outfile << "Log PCO2                 = " << pco << endl;
                    outfile << "density                  = " << dens << endl;
                    outfile << "molal/molar factor       = " << ee/1000 << endl;
                    if (kinvariant == 0) {
                        outfile << "increment (%)            = " << xi * 100 << endl;
                    } else {
                        outfile << "increment (moles)        = " << xinv << endl;
                    }
                    outfile << "number of steps          = " << ncpt << endl;
                    outfile << endl;
                    outfile.close();
                }

                if (output == 1) {
                    q_S = "";
                    y_S = to_string(fc);
                    for (int k=1; k<=nm; k++) {
                        x_S = to_string(minp[k]);
                        if (system == "c") {
                            if (lmin[k] == 1 and lmin0[k] == 0) {
                                q_S = "start of precipitation of " + mineral_S[k] + " at fc = " + y_S;
                            } else if (lmin[k] == 1 and lmin1[k] == 0) {
                                if (min[k] < min0[k]) {
                                    lmin1[k] = 1;
                                    q_S = "end of precipitation and start of dissolution of " + mineral_S[k] + " at fc = " + y_S;
                                }
                            } else if (lmin[k] == 1 and lmin1[k] == 1) {
                                if (min[k] > min0[k]) {
                                    lmin1[k] = 0;
                                    q_S = "end of dissolution and start of precipitation of " + mineral_S[k] + " at fc = " + y_S;
                                }
                            } else if (lmin[k] == 0 and lmin1[k] == 1 and lmin0[k] == 1) {
                                lmin1[k] = 0;
                                q_S = "end of dissolution and of saturation of " + mineral_S[k] + " at fc = " + y_S;
                            } else if (lmin[k] == 0 and lmin1[k] == 0 and lmin0[k] == 1) {
                                lmin1[k] = 0;
                                q_S = "end of saturation of " + mineral_S[k] + " at fc = " + y_S;
                            }
                        } else if (system == "o") {
                            if (lmin[k] == 1 and lmin0[k] == 0) {
                                q_S = "start of precipitation of " + mineral_S[k] + " at fc = " + y_S;
                            } else if (lmin[k] == 0 and lmin0[k] == 1) {
                                q_S = "end of precipitation of " + mineral_S[k] + " at fc = " + y_S + ": moles = " + x_S;
                            }
                        }
                        if (q_S != "") {
                            outfile.open(event_file, ios::app);
                            outfile << q_S << endl;
                            outfile.close();
                            q_S = "";
                        }
                    }

                    if (ncpt == 1 or my_S != my0_S or npasf == output_step) {
                        if (units == "molal") ee = 1000;
                        if (units == "molar") ee = 1000000 * dens / (1000 + std);
                        outfile.open(chem_file, ios::app);
                            if (my_S == "") {
                                outfile << "No_minerals,";
                            } else {
                                outfile << my_S + ",";
                            }
                            outfile << fc << ",";
                            outfile << act[11] << "," << (std/(std+1000)) << ",";
                            outfile << mwev << "," << dens << "," << -log10(act[13]) << ",";
                            outfile << alk * ee << ",";
                            for (int i=1; i<=10; i++) {
                                if (tot[i] > 0) {
                                    s = 0;
                                    for (int j=1; j<=n; j++) {
                                        s += molal[j] * kmat[i][j];
                                    }
                                    outfile << s * ee << ",";
                                }
                            }
                            outfile << std * ee << endl;
                        outfile.close();

                        outfile.open(min_file, ios::app);
                            outfile << fc << ",";
                            if (system == "c") {
                                for (int i=1; i<=nm-1; i++) {
                                    if (min[i] >= 0) {
                                        outfile << min[i] << ",";
                                    } else if (min[i] < 0) {
                                        outfile << 0.0 << ",";
                                    }
                                }
                                if (min[nm] >= 0) {
                                    outfile << min[nm];
                                } else if (min[nm] < 0) {
                                    outfile << 0.0;
                                }
                            } else if (system == "o") {
                                for (int i=1; i<=nm-1; i++) {
                                    if (minp[i] >= 0) {
                                        outfile << minp[i] << ",";
                                    } else if (minp[i] < 0) {
                                        outfile << 0.0 << ",";
                                    }
                                }
                                if (minp[nm] >= 0) {
                                    outfile << minp[nm];
                                } else if (minp[nm] < 0) {
                                    outfile << 0.0;
                                }
                            }
                            outfile << endl;
                        outfile.close();
                    }
                }
                if (stdmax > 0 and std * ee >= stdmax) {
                    stop_simulation();
                }
                if (mwev > 0 and kinvariant == 0 and abs(act[11] - act0[11]) / (act[11] + act0[11]) < 0.0000000001) {
                    if (system == "c") {
                        nu = 0;
                        for (int k=1; k<=nm; k++) {
                            if (min[k] > 0 and min[k] < min0[k]) nu = 1;
                        }
                        if (nu == 0) {
                            if (nminer == nbmin) q_S = "invariant system / eutectic point / end";
                            if (nminer < nbmin) q_S = "invariant system /pseudo-eutectic point / end";
                        } else if (nu == 1) {
                            if (nperitec == 0) peritec();
                            if (nperitec == 0) {
                                if (nminer == nbmin) q_S = "invariant system / peritectic point / end";
                                if (nminer < nbmin) q_S = "invariant system / pseudo-peritectic point / end";
                            } else if (nperitec == 1) {
                                if (nminer == nbmin) q_S = "invariant system / peritectic / passing over";
                                if (nminer < nbmin) q_S = "invariant system / pseudo-peritectic / passing over";
                            }
                        }
                    } else if (system == "o") {
                        if (nminer == nbmin) q_S = "invariant system / eutectic / end";
                        if (nminer < nbmin) q_S = "invariant system / pseudo-eutectic / end";
                    }

                    if (verbose == 1) {
                        cout << endl;
                        cout << q_S << endl;
                        if (q_S.find("pseudo") != string::npos) {
                            q1_S = "maximum number of munerals allowed by the phase rule = ";
                            q2_S = "number of minerals in equilibrium with the invariant system = ";
                            cout << q1_S << nbmin << endl;
                            cout << q2_S << nbmin << endl;
                        }
                    }
                    if (output == 1 and q_S != q0_S) {
                        outfile.open(log_file, ios::app);
                        outfile << endl;
                        outfile << q_S << endl;
                        if (q_S.find("pseudo") != string::npos) {
                            q1_S = "maximum number of munerals allowed by the phase rule = ";
                            q2_S = "number of minerals in equilibrium with the invariant system = ";
                            outfile << q1_S << nbmin << endl;
                            outfile << q2_S << nbmin << endl;
                        }
                        outfile.close();

                        outfile.open(event_file, ios::app);
                            outfile << q_S << " at fc = " << fc << endl;
                            if (q_S.find("pseudo") != string::npos) {
                                q1_S = "maximum number of munerals allowed by the phase rule = ";
                                q2_S = "number of minerals in equilibrium with the invariant system = ";
                                outfile << q1_S << nbmin << endl;
                                outfile << q2_S << nbmin << endl;
                            }
                            q0_S = q_S;
                        outfile.close();
                    }
                    if (q_S.find("end") != string::npos) {
                        stop_simulation();
                    }
                } else if (kinvariant == 0 and abs(act[11] - act0[11]) / (act[11] + act0[11]) > 0.0000000001) {
                    nperitec = 0;
                    q_S = ""; q0_S = ""; q1_S = ""; q2_S = "";
                }

                if (npasi == print_step) npasi = 0;
                if (npasf == output_step) npasf = 0;
                if (my_S != my0_S) my0_S = my_S;
                for (int i=0; i<=n; i++) {
                    mol0[i] = mol[i];
                    molal0[i] = molal[i];
                    gact0[i] = gact[i];
                    act0[i] = act[i];
                }
                for (int i=1; i<=ntot; i++) {
                    tot0[i] = tot[i];
                }
                for (int k=1; k<=nm; k++) {
                    lmin0[k] = lmin[k];
                    pai0[k] = pai[k];
                    min0[k] = min[k];
                    minp0[k] = minp[k];
                }
                fi0 = fi;
                mwev0 = mwev;
                nminer0 = nminer;
                if (kinvariant == 0 and initdeseq == 0) {
                    if (inc_S == "auto") {
                        xi = (51 - 8 * log10(std * 1000)) / 700;
                        xi0 = xi;
                    } else if (inc_S == "manu") {
                        xi = xi0;
                    }
                    mwev = mwev + mol[11] * xi;
                    tot[11] = tot[11] - 2 * mol[11] * xi;
                }
            }

            LOOP2000:
            if (verbose == 1) {
                cout << "Free energy minimization" << endl;
            }
            if (kinvariant == 2) {
                if (wmin[kinvar[1]][ncm] > wmin[kinvar[2]][ncm]) {
                    ksupprim = kinvar[1];
                } else {
                    ksupprim = kinvar[2];
                }
            } else if (kinvariant > 2) {
                gmin = 0;
                for (int ii=1; ii<=kinvariant; ii++) {
                    for (int i=0; i<=n; i++) {
                        mol[i] = mol0[i];
                    }
                    for (int i=1; i<=ntot; i++) {
                        tot[i] = tot0[i];
                    }
                    mwev = mwev + mol[11] * xi;
                    tot[11] = tot[11] - 2 * mol[11] * xi;
                    for (int k=1; k<=nm; k++) {
                        lmin[k] = lmin0[k];
                        min[k] = min0[k];
                        minp[k] = minp0[k];
                    }
                    if (kinvar[ii] > 0 and kinvar[ii] != kneuf and min[kinvar[ii]] >= 0) {
                        l = lmin[kinvar[ii]]; m = min[kinvar[ii]]; p = psol[kinvar[ii]];
                        lmin[kinvar[ii]] = 0; min[kinvar[ii]] = 0;
                        psol[kinvar[ii]] = 1.e+50;
                        if (verbose == 1) {
                            cout << "mineral removed: " << mineral_S[kinvar[ii]] << endl;
                        }
                        ncmptinv = 0;
                        nw = 1;
                        while (nw != 0) {
                            ncmptinv += 1;
                            for (int i=0; i<=n; i++) gact0[i] = gact[i];
                            evp_actp();
                            for (int i=0; i<=n; i++) gact[i] = (gact0[i] + gact[i]) / 2;
                            for (int i=0; i<=n; i++) {
                                molal[i] = mol[i] * mh2o / mol[11];
                                act[i] = molal[i] * gact[i];
                            }
                            for (int k=1; k<=nm; k++) {
                                pai[k] = 1;
                                for (int i=1; i<=ncm; i++) {
                                    pai[k] = pai[k] * pow(act[i], wmin[k][i]);
                                }
                                if (pai[k] >= psol[k]) {
                                    if (min[k] >= 0) {
                                        lmin[k] = 1;
                                    } else if (min[k] < 0) {
                                        lmin[k] = 0;
                                        min[k] = 0;
                                    }
                                } else if (pai[k] < psol[k]) {
                                    if (min[k] <= 0) {
                                        lmin[k] = 0;
                                        min[k] = 0;
                                    } else if (min[k] > 0) {
                                        lmin[k] = 1;
                                    }
                                }
                            }
                            mineraux_S = "";
                            nminer = 0;
                            for (int k=1; k<=nm; k++) {
                                if (lmin[k] == 1) {
                                    nminer += 1;
                                    mineraux_S += "_" + mineral_S[k];
                                }
                            }
                            if (verbose == 1) {
                                cout << ncmptinv << mineraux_S << endl;
                            }
                            for (int i=0; i<=n; i++) mol1[i] = mol[i];
                            reseq();
                            for (int i=0; i<=n; i++) molal[i] = mol[i] * mh2o / mol[11];
                            nw = 0;
                            for (int i=0; i<=n; i++) {
                                if (abs(mol1[i] - mol[i]) > pkmol) nw = 1;
                            }
                        }

                        g = 0;
                        for (int i=0; i<=ncm; i++) {
                            if (mol[i] > 0) g = g + mol[i] * (mu[i] + log(act[i]));
                        }
                        for (int k=1; k<=nm; k++) {
                            g = g + min[k] * mum[k];
                        }
                        if (verbose == 1) {
                            cout << "g = " << g << endl;
                            for (int i=1; i<=kinvariant; i++) {
                                if (i != ii) {
                                    cout << mineral_S[kinvar[i]] << endl;
                                }
                            }
                            cout << endl;
                        }
                        if (g < gmin) {
                            gmin = g;
                            ksupprim = kinvar[ii];
                        }
                        lmin[kinvar[ii]] = l; min[kinvar[ii]] = m; psol[kinvar[ii]] = p;
                    }
                }
            }
            for (int i=0; i<=n; i++) {
                mol[i] = mol0[i];
            }
            for (int i=1; i<=ntot; i++) {
                tot[i] = tot0[i];
            }
            for (int k=1; k<=nm; k++) {
                lmin[k] = lmin0[k];
                min[k] = min0[k];
                minp[k] = minp0[k];
            }
            mwev = mwev0;
            nminer = nminer0;
            mwev = mwev + mol[11] * xi;
            tot[11] = tot[11] - 2 * mol[11] * xi;
            kinvariant = 0;
            lmin[ksupprim] = 0; min[ksupprim] = 0; psol[ksupprim] = 1.e+50;
            min0[ksupprim] = 0;
            if (verbose == 1) {
                cout << "mineral definitely removed: " << mineral_S[ksupprim] << endl;
            }
            goto LOOP500;

        STOP:
            stop_simulation();
        }

        void read_transfer() {
            // read stockage
            infile.open("resources/stockage");
            infile >> tra_file;
            infile.close(); 

            // read transfer file
            infile.open(tra_file);
                infile >> temp; infile >> tinit;
                infile >> ph; infile >> phinit; 
                infile >> po; infile >> poinit;
                infile >> diltot; infile >> constituant_S;
                nbmin = 0;
                for (int i=1; i<=10; i++) {
                    infile >> totinit[i];
                    if (totinit[i] != 0) nbmin += 1;
                }
                // allocate ica and kinvar based on nbmin
                for (int i=0; i<=n+nbmin; i++) ica.push_back(0);
                for (int i=0; i<=nbmin+3; i++) kinvar.push_back(0);
                for (int i=0; i<=n; i++) {
                    infile >> mol[i];
                    mol0[i] = mol[i];
                }
                infile >> system; infile >> xi;
                infile >> output_step; infile >> output;
                infile >> print_step; infile >> verbose;
                infile >> units;
                infile >> log_file; infile >> chem_file;
                infile >> event_file; infile >> min_file;
                infile >> miner_S; infile >> stdmax;
                infile >> pkmol; infile >> pkeq;
            infile.close();

        }

        void evp_actp() {
            double cat[10] = {0}; double ani[12] = {0}; double h[4] = {0};
            double at, bt, ct, dt, et;
            double u, z, w, v, f, fj, co, s;
            int nc, na, nn;

            cat[1] = mol[1]; cat[2] = mol[2]; cat[3] = mol[3];
            cat[4] = mol[4]; cat[5] = mol[22]; cat[6] = mol[5];
            cat[7] = mol[18]; cat[8] = mol[23]; cat[9] = mol[13];
            ani[1] = mol[6]; ani[2] = mol[7]; ani[3] = mol[25];
            ani[4] = mol[15]; ani[5] = mol[14]; ani[6] = mol[12];
            ani[7] = mol[24]; ani[8] = mol[19]; ani[9] = mol[20];
            ani[10] = mol[21]; ani[11] = mol[8];
            h[1] = mol[10]; h[2] = mol[9]; h[3] = mol[0];

            for (int i=1; i<=9; i++) {
                cat[i] = cat[i] * mh2o / mol[11];
            }
            for (int i=1; i<=11; i++) {
                ani[i] = ani[i] * mh2o / mol[11];
            }
            for (int i=1; i<=3; i++) {
                h[i] = h[i] * mh2o / mol[11];
            }

            {
                infile.open("resources/coefft4");
                getline(infile, line);
                stringstream linestream(line);
                string value;
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
                        getline(infile, line); stringstream linestream(line); string value;

                        getline(linestream, value, ','); string x_S = value;
                        getline(linestream, value, ','); nzc[i] = stoi(value);
                    }
                    for (int i=1; i<=na; i++) {
                        getline(infile, line); stringstream linestream(line); string value;

                        getline(linestream, value, ','); string x_S = value;
                        getline(linestream, value, ','); nza[i] = stoi(value);
                    }
                    {
                        getline(infile, line); stringstream linestream(line); string value;

                        getline(linestream, value, ','); string x_S = value;
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
                                getline(infile, line); stringstream linestream(line); string value;

                                getline(linestream, value, ','); string x_S = value;
                                getline(linestream, value, ','); at = stod(value);
                                getline(linestream, value, ','); bt = stod(value);
                                getline(linestream, value, ','); ct = stod(value);
                                getline(linestream, value, ','); dt = stod(value);
                                getline(linestream, value, ','); et = stod(value);

                                b0[i][j] = temperature(at, bt, ct, dt, et, temp);
                            }
                            {
                                getline(infile, line); stringstream linestream(line); string value;

                                getline(linestream, value, ','); string x_S = value;
                                getline(linestream, value, ','); at = stod(value);
                                getline(linestream, value, ','); bt = stod(value);
                                getline(linestream, value, ','); ct = stod(value);
                                getline(linestream, value, ','); dt = stod(value);
                                getline(linestream, value, ','); et = stod(value);

                                b1[i][j] = temperature(at, bt, ct, dt, et, temp);
                            }
                            {
                                getline(infile, line); stringstream linestream(line); string value;

                                getline(linestream, value, ','); string x_S = value;
                                getline(linestream, value, ','); at = stod(value);
                                getline(linestream, value, ','); bt = stod(value);
                                getline(linestream, value, ','); ct = stod(value);
                                getline(linestream, value, ','); dt = stod(value);
                                getline(linestream, value, ','); et = stod(value);

                                b2[i][j] = temperature(at, bt, ct, dt, et, temp);
                            }
                            {
                                getline(infile, line); stringstream linestream(line); string value;

                                getline(linestream, value, ','); string x_S = value;
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
                            getline(infile, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = value;
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
                            getline(infile, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = value;
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
                                getline(infile, line); stringstream linestream(line); string value;

                                getline(linestream, value, ','); string x_S = value;
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
                                getline(infile, line); stringstream linestream(line); string value;

                                getline(linestream, value, ','); string x_S = value;
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
                            getline(infile, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = value;
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
                            getline(infile, line); stringstream linestream(line); string value;

                            getline(linestream, value, ','); string x_S = value;
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
                    xip[2][9][1] = -0.0102;
                    xip[2][1][2] = 0.046;
                }
                
                infile.close();
            }

            bp0 = 1.2e0;

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

            u = 0; z = 0;
            for (int i=1; i<=nc; i++) {
                u += cat[i] * pow(nzc[i], 2); z += cat[i] * nzc[i];
            }
            for (int j=1; j<=na; j++) {
                u += ani[j] * pow(nza[j], 2); z += ani[j] * nza[j];
            }
            fi = u / 2; fj = sqrt(fi);
            u = 6 * ap0 * fj;
            for (int i=1; i<=nc-1; i++) {
                for (int j=i+1; j<=nc; j++) {
                    if (nzc[i] == nzc[j]) {
                        ec[i][j] = 0; fc[i][j] = 0;
                    } else {
                        xc[i][j] = 2 * u; xc[i][i] = pow(nzc[i], 2) * u; xc[j][j] = pow(nzc[j], 2) * u;
                        ec[i][j] = (j0(xc[i][j]) - j0(xc[i][i]) / 2 - j0(xc[j][j]) / 2) / fi / 2;
                        fc[i][j] = ((xc[i][j] * j1(xc[i][j]) - xc[i][i] * j1(xc[i][i]) / 2 - 
                                    xc[j][j] * j1(xc[j][j]) / 2) / pow(fi, 2) / 4 - ec[i][j] / fi);
                        ec[j][i] = ec[i][j]; fc[j][i] = fc[i][j];
                    }
                }
            }
            for (int i=1; i<=na-1; i++) {
                for (int j=i+1; j<=na; j++) {
                    if (nza[i] == nza[j]) {
                        ea[i][j] = 0; fa[i][j] = 0;
                    } else {
                        xa[i][j] = 2 * u; xa[i][i] = pow(nza[i], 2) * u; xa[j][j] = pow(nza[j], 2) * u;
                        ea[i][j] = (j0(xa[i][j]) - j0(xa[i][i]) / 2 - j0(xa[j][j]) / 2) / fi / 2;
                        fa[i][j] = ((xa[i][j] * j1(xa[i][j]) - xa[i][i] * j1(xa[i][i]) / 
                                        2 - xa[j][j] * j1(xa[j][j]) / 2) / pow(fi, 2) / 4 - ea[i][j] / fi);
                        ea[j][i] = ea[i][j]; fa[j][i] = fa[i][j];
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
                    f += cat[i] * ani[j] * bp[i][j];
                }
            }
            for (int i=1; i<=nc-1; i++) {
                for (int j=i+1; j<=nc; j++) {
                    f += cat[i] * cat[j] * pp[i][j];
                }
            }
            for (int i=1; i<=na-1; i++) {
                for (int j=i+1; j<=na; j++) {
                    f += ani[i] * ani[j] * qp[i][j];
                }
            }
            for (int ii=1; ii<=nc; ii++) {
                u = pow(nzc[ii], 2) * f;
                for (int j=1; j<=na; j++) {
                    u += ani[j] * (b[ii][j] * 2 + z * cc[ii][j]);
                } 
                for (int i=1; i<=nc; i++) {
                    if (i != ii) {
                        v = 0;
                        for (int j=1; j<=na; j++) {
                            v += ani[j] * scp[ii][i][j];
                        }
                        u += cat[i] * (p[ii][i] * 2 + v);
                    }
                }
                for (int i=1; i<=na-1; i++) {
                    for (int j=i+1; j<=na; j++) {
                        u += ani[i] * ani[j] * sap[i][j][ii];
                    }
                }
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += cat[i] * ani[j] * cc[i][j] * nzc[ii];
                    }
                }
                for (int i=1; i<=nn; i++) {
                    u += h[i] * lcp[i][ii] * 2;
                }
                for (int k=1; k<=nn; k++) {
                    for (int j=1; j<=na; j++) {
                        u += h[k] * ani[j] * xip[k][ii][j];
                    }
                }
                gc[ii] = exp(u);
            }
            for (int jj=1; jj<=na; jj++) {
                u = pow(nza[jj], 2) * f;
                for (int i=1; i<=nc; i++) {
                    u += cat[i] * (b[i][jj] * 2 + z * cc[i][jj]);
                }
                for (int i=1; i<=na; i++) {
                    if (i != jj) {
                        v = 0;
                        for (int j=1; j<=nc; j++) {
                            v += cat[j] * sap[jj][i][j];
                        }
                        u += ani[i] * (q[jj][i] * 2 + v);
                    }
                }
                for (int i=1; i<=nc-1; i++) {
                    for (int j=i+1; j<=nc; j++) {
                        u += cat[i] * cat[j] * scp[i][j][jj];
                    }
                }
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += cat[i] * ani[j] * cc[i][j] * nza[jj];
                    }
                }
                for (int j=1; j<=nn; j++) {
                    u += h[j] * lap[j][jj];
                }
                for (int k=1; k<=nn; k++) {
                    for (int i=1; i<=nc; i++) {
                        u += h[k] * cat[i] * xip[k][i][jj];
                    }
                }
                ga[jj] = exp(u);
            }
            for (int k=1; k<=nn; k++) {
                u = 0;
                for (int i=1; i<=nc; i++) {
                    u += cat[i] * lcp[k][i] * 2;
                }
                for (int j=1; j<=na; j++) {
                    u += ani[j] * lap[k][j] * 2;
                }
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += cat[i] * ani[j] * xip[k][i][j];
                    }
                }
                gn[k] = exp(u);
            }
            u = -ap0 * pow(fi, 1.5e0) / (1 + bp0 * fj);
            for (int i=1; i<=nc; i++) {
                for (int j=1; j<=na; j++) {
                    u += cat[i] * ani[j] * (bf[i][j] + z * cc[i][j]);
                }
            }
            for (int i=1; i<=nc-1; i++) {
                for (int j=i+1; j<=nc; j++) {
                    v = 0;
                    for (int k=1; k<=na; k++) {
                        v += ani[k] * scp[i][j][k];
                    }
                    u += cat[i] * cat[j] * (pf[i][j] + v);
                }
            }
            for (int i=1; i<=na-1; i++) {
                for (int j=i+1; j<=na; j++) {
                    v = 0;
                    for (int k=1; k<=nc; k++) {
                        v += cat[k] * sap[i][j][k];
                    }
                    u += ani[i] * ani[j] * (qf[i][j] + v);
                }
            }
            for (int k=1; k<=nn; k++) {
                for (int i=1; i<=nc; i++) {
                    u += h[k] * cat[i] * lcp[k][i];
                }
            }
            for (int k=1; k<=nn; k++) {
                for (int j=1; j<=na; j++) {
                    u += h[k] * ani[j] * lap[k][j];
                }
            }
            for (int k=1; k<=nn; k++) {
                for (int i=1; i<=nc; i++) {
                    for (int j=1; j<=na; j++) {
                        u += h[k] * cat[i] * ani[j] * xip[k][i][j];
                    }
                }
            }

            s = 0;
            for (int i=1; i<=nc; i++) s += cat[i];
            for (int j=1; j<=na; j++) s += ani[j];

            co = 1 + 2 * u / s;
            aw = exp(-s * co / mh2o);
            gact[0] = gn[3]; gact[1] = gc[1]; gact[2] = gc[2]; gact[3] = gc[3];
            gact[4] = gc[4]; gact[22] = gc[5]; gact[5] = gc[6];
            gact[18] = gc[7]; gact[23] = gc[8]; gact[13] = gc[9];
            gact[6] = ga[1]; gact[7] = ga[2]; gact[25] = ga[3];
            gact[15] = ga[4]; gact[14] = ga[5]; gact[12] = ga[6];
            gact[24] = ga[7]; gact[19] = ga[8]; gact[20] = ga[9];
            gact[21] = ga[10]; gact[8] = ga[11];
            gact[10] = aw * aw * pow(gn[1], log(10)); gact[9] = gn[2]; 
            gact[16] = 1e0; gact[17] = 1e0; gact[11] = aw / mh2o;
            ndepact = 1;
        }

        void evp_density() {
            int ncdens, nadens;
            double u, v;
            string x_S;

            ncdens = 5; nadens = 5;
            double s[ncdens+1][nadens+1] = {0};
            double cat[ncdens+1] = {0}; double ani[nadens+1] = {0};
            int ic[ncdens+1] = {0}; int ia[nadens+1] = {0};

            for (int i=1; i<=8; i++) {
                if (i <= ncdens) cat[i] = mol[i] / mol[11] * mh2o;
                if (i > ncdens) ani[i-ncdens] = mol[i] / mol[11] * mh2o;
            }
            ani[4] = mol[15] / mol[11] * mh2o;
            ani[5] = mol[14] / mol[11] * mh2o;

            for (int i=1; i<=ncdens; i++) ic[i] = nch[i];
            for (int i=1; i<=3; i++) ia[i] = -nch[i+5];
            ia[4] = -nch[15]; ia[5] = -nch[14];
            if (ncpt == 1) {
                for (int i=0; i<=ncdens; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=nadens; k++) temp1d.push_back(0); au.push_back(temp1d);
                }
                for (int i=0; i<=ncdens; i++) {
                    vector<double> temp1d;
                    for (int k=0; k<=nadens; k++) temp1d.push_back(0); bu.push_back(temp1d);
                }
                infile.open("resources/densite");
                for (int i=1; i<=5; i++) {
                    for (int j=1; j<=5; j++) {
                        getline(infile, line); stringstream linestream(line);
                        getline(linestream, value, ','); x_S = value;
                        getline(linestream, value, ','); u = stod(value);
                        getline(linestream, value, ','); v = stod(value);
                        getline(linestream, value, ','); au[i][j] = stod(value);
                        getline(linestream, value, ','); bu[i][j] = stod(value);
                    }
                }
                infile.close();
            }
            dens = 1e0; u = 0;
            for (int j=1; j<=nadens; j++) u += ani[j] * ia[j];
            for (int i=1; i<=ncdens; i++) {
                for (int j=1; j<=nadens; j++) {
                    s[i][j] = int((ic[i] + ia[j]) / 2) * cat[i] * ani[j] / u;
                    dens += au[i][j] * s[i][j] + bu[i][j] * pow(s[i][j], 2);
                }
            }
        }

        void evp_invar() {
            double swap, ah2o;
            int ninvar, l, n1, n2, det, nz, ii, jj, kk, det1, n3, n4, ncond;

            kinvariant = 0; ncm = 15;
            ninvar = nbmin + 3;
            vector<string> minv_S, minvar_S;

            for (int i=0; i<=ninvar; i++) {
                minv_S.push_back("");
                minvar_S.push_back("");
            }

            int kinv[ninvar+1] = {0}; 
            double psminv[ninvar+1] = {0}; double winv[ninvar+1][ncm+1] = {0};
            double psminvar[ninvar+1] = {0}; double t0[ninvar+1][ninvar+1] = {0};
            double t1[ninvar+1][ninvar+1] = {0}; double t2[ninvar+1][ninvar+1] = {0};
            double t3[ninvar+1][ncm+1] = {0}; double t4[ncm+1][ncm+1] = {0};
            double tt4[ncm+1] = {0};

            for (int k=1; k<=3; k++) {
                psminv[k] = log10(psc[k]);
            }
            winv[1][13] = 1; winv[1][12] = 1; winv[1][15] = -1;
            winv[2][13] = 1; winv[2][11] = -1; winv[2][14] = 1;
            winv[3][13] = 1; winv[3][11] = 1; winv[3][15] = -1;

            n1 = 3;
            for (int k=1; k<=nm; k++) {
                if (lmin[k] == 1) {
                    n1 += 1;
                    kinv[n1] = k;
                    minv_S[n1] = mineral_S[k];
                    psminv[n1] = log10(psol[k]);
                    for (int j=1; j<=ncm; j++) {
                        winv[n1][j] = wmin[k][j];
                    }
                    swap = winv[n1][11];
                    winv[n1][11] = winv[n1][15];
                    winv[n1][15] = swap;
                }
            }

            for (int i=1; i<=n1; i++) {
                for (int j=i; j<=n1; j++) {
                    t1[i][j] = 0;
                    for (int k=1; k<=ncm-1; k++) {
                        t1[i][j] = t1[i][j] + winv[i][k] * winv[j][k];
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
                n3 = 0; n2 = n1 - 1;
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
                            if (t2[i][k-1] != 0) {
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
                        for (int j=1; j<=ncm; j++) {
                            t3[n3][j] = winv[kk][j];
                        }
                    }
                }

                n4 = ncm;
                for (int j=ncm; j>=1; --j) {
                    u = 0;
                    for (int i=1; i<=n3; i++) {
                        u += pow(t3[i][j], 2);
                    }
                    if (u == 0) {
                        for (int k=j+1; k<=n4; k++) {
                            for (int i=1; i<=n3; i++) {
                                t3[i][k-1] = t3[i][k];
                            }
                        }
                        n4 -= 1;
                    }
                }

                for (int i=1; i<=n4; i++) {
                    for (int j=i; j<=n4; j++) {
                        t4[i][j] = 0;
                        for (int k=1; k<=n4; k++) {
                            t4[i][j] = t4[i][j] + t3[k][i] * t3[k][j];
                            t4[j][i] = t4[i][j];
                        }
                    }
                }
                for (int i=1; i<=n4; i++) {
                    tt4[i] = 0;
                    for (int k=1; k<=n3; k++) {
                        tt4[i] = tt4[i] + t3[k][i] * psminvar[k];
                    }
                }

                for (int k=2; k<=n4; k++) {
                    for (int i=k; i<=n4; i++) {
                        if (abs(t4[i][k-1]) > epsilon) {
                            u = t4[k-1][k-1] / t4[i][k-1];
                            for (int j=k; j<=n4; j++) {
                                t4[i][j] = t4[k-1][j] - t4[i][j] * u;
                                if (abs(t4[i][j]) < epsilon) t4[i][j] = 0;
                            }
                            tt4[i] = tt4[k-1] - tt4[i] * u;
                        }
                    }
                }

                ah2o = pow(10, (tt4[n4] / t4[n4][n4]));
                kinvariant = n3;
                for (int i=kinvariant; i>=1; --i) {
                    if (kinvar[i] == 0) {
                        for (int k=1; k<=kinvariant-1; k++) {
                            kinvar[k] = kinvar[k+1];
                        }
                        kinvariant -= 1;
                    }
                }
                if (verbose == 1) {
                    cout << "invariant system constrained by: " << endl;
                    for (int k=1; k<=kinvariant; k++) {
                        cout << " " << mineral_S[kinvar[k]] << endl;
                    }
                    cout << endl;
                    cout << "invariant aH2O = " << ah2o << endl;
                    cout << "simulation aH2O = " << aw << endl;
                }
                if (output == 1) {
                    outfile.open(event_file, ios::app);
                        outfile << "invariant system constrained by: " << endl;
                        for (int k=1; k<=kinvariant; k++) {
                            outfile << " " << mineral_S[kinvar[k]] << endl;
                        }
                        outfile << endl;
                        outfile << "invariant aH2O = " << ah2o << "      ";
                        outfile << "simulation aH2O = " << aw << endl;
                    outfile.close();
                }
            }
        }

        void reseq() {
            double s, sh, p, u;
            int nconv, nt, nu, kk, k, j, l, i, ni, nj;
            
            double z[n+nminer+1][n+nminer+1] = {0};
            double zz[n+nminer+1] = {0}; double xx[n+nminer+1] = {0};

            nconv = 1;
            ncm = 15; nt = n;
            nu = 1;
            while (nu != 0) {
                for (int i=1; i<=nt; i++) {
                    for (int j=1; j<=nt; j++) {
                        z[i][j] = 0;
                    }
                    zz[i] = 0;
                }
                for (int i=1; i<=12; i++) {
                    for (int j=1; j<=n; j++) {
                        if (mol[j] != 0) z[i][j] = kmat[i][j];
                    }
                }
                for (int i=13; i<=n; i++) {
                    kmat[i][0] = 0;
                    for (int j=1; j<=n; j++) {
                        if (j != 11) kmat[i][0] = kmat[i][0] + kmat[i][j];
                    }
                }
                for (int i=13; i<=n; i++) {
                    p = 1; u = 0;
                    for (int j=1; j<=n; j++) {
                        if (mol[j] != 0 and j != 11) {
                            z[i][j] = kmat[i][j] / mol[j];
                            p = p * pow(gact[j], kmat[i][j]);
                            u += kmat[i][j] * log(mol[j]);
                        } else if (j == 11) {
                            z[i][j] = -kmat[i][0] / mol[j];
                            p = p * pow(aw, kmat[i][j]);
                            u -= kmat[i][0] * log(mol[j]);
                        } else if (mol[j] == 0) {
                            z[i][j] = 0;
                        }
                    }
                    p = p * pow(mh2o, kmat[i][0]);
                    zz[i] = log(psc[i-12]) - log(p) - u;
                }

                l = 0;
                for (int k=1; k<=nm; k++) {
                    if (lmin[k] == 1) {
                        l += 1;
                        ica[n+l] = 1;
                        wmin[k][0] = 0;
                        p = 1; u = 0;
                        for (int j=1; j<=ncm; j++) {
                            if (j != 11) wmin[k][0] = wmin[k][0] + wmin[k][j];
                        }
                        for (int j=1; j<=ncm; j++) {
                            if (j != 11 and mol[j] > 0) {
                                z[n+l][j] = wmin[k][j] / mol[j];
                                p = p * pow(gact[j], wmin[k][j]);
                                u += wmin[k][j] * log(mol[j]);
                            } else if (j == 11) {
                                z[n+l][j] = -wmin[k][0] / mol[j];
                                p = p * pow(aw, wmin[k][j]);
                                u -= wmin[k][0] * log(mol[j]);
                            }
                        }
                        p = p * pow(mh2o, wmin[k][0]);
                        zz[n+l] = log(psol[k]) - log(p) - u;

                        sh = 0;
                        for (int j=1; j<=ncm; j++) {
                            sh += wmin[k][j] * kmat[11][j];
                        }
                        for (int i=1; i<=10; i++) {
                            z[i][n+l] = wmin[k][i];
                        }
                        z[11][n+l] = sh;
                    }
                }
                nt = n + l;

                for (int i=1; i<=10; i++) {
                    u = 0;
                    for (int j=1; j<=n; j++) {
                        u += kmat[i][j] * mol[j];
                    }
                    for (int k=1; k<=nm; k++) {
                        if (lmin[k] == 1) {
                            u += min[k] * wmin[k][i];
                        }
                    }
                    zz[i] = tot[i] - u;
                }
                u = 0;
                for (int j=1; j<=n; j++) {
                    u += kmat[11][j] * mol[j];
                }
                for (int j=1; j<=ncm; j++) {
                    for (int k=1; k<=nm; k++) {
                        u += wmin[k][j] * kmat[11][j] * min[k];
                    }
                }

                zz[11] = tot[11] - u;
                u = 0;
                for (int j=1; j<=n; j++) {
                    u += z[12][j] * mol[j];
                }
                zz[12] = tot[12] - u;

                for (int k=1; k<=10; k++) {
                    if (tot[k] == 0) {
                        ica[k] = 0;
                        for (int j=k+1; j<=n; j++) {
                            if (kmat[k][j] != 0) ica[j] = 0;
                        }
                    }
                }

                ni = nt; nj = nt;
                for (int k=nt; k>=1; --k) {
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
                        if (abs(z[i][k-1]) > epsilon) {
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

                for (int k=1; k<=nt; k++) {
                    if (ica[k] == 0) {
                        for (int i=ni; i>=k; --i) {
                            xx[i+1] = xx[i];
                        }
                        xx[k] = 0;
                        ni += 1;
                    }
                }
                nu = 1;
                while (nu != 0) {
                    for (int i=1; i<=n; i++) {
                        if (ica[i] == 1) {
                            nu = 0;
                            if (mol[i] + xx[i] / nconv < 0) {
                                nconv += 1;
                                nu = 1;
                                break;
                            }
                        }
                    }
                    if (nconv >= max_conv) {
                        if (verbose == 1) {
                            cout << endl;
                            cout << "The equation system diverges: end of simulation" << endl;
                        }
                        if (output == 1) {
                            outfile.open(event_file, ios::app);
                                outfile << endl;
                                outfile << "The equation system diverges: end of simulation" << endl;
                            outfile.close();
                        }
                        stop_simulation();
                    }
                }
                for (int i=1; i<=n; i++) {
                    mol[i] = mol[i] + xx[i] / nconv;
                }
                i = n;
                for (int k=1; k<=nm; k++) {
                    if (lmin[k] == 1) {
                        i += 1;
                        min[k] = min[k] + xx[i] / nconv;
                    }
                }
                nu = 0;
                for (int i=1; i<=nt; i++) {
                    if (ica[i] == 1) {
                        if (abs(xx[i]) > pkeq) nu = 1;
                    }
                }
            }
        }

        void reseqinv() {
            int nt, k, j, i, ii, kk, nmin;
            double hmin, s, swap, u;

            nt = ntot - 1;
            double t[nt+1][nt+1] = {0}; double tt[nt+1] = {0};
            double t0[nt+1][nt+1] = {0}; double tt0[nt+1] = {0};
            double xx[nt+1] = {0};

            if (ninv == 0) {
                hmin = 1000;
                for (int k=1; k<=nm; k++) {
                    for (int kk=1; kk<=kinvariant; kk++) {
                        cout << kk << endl;
                        if (k == kinvar[kk] and min[k] > 0) {
                            s = 0;
                            for (int i=1; i<=15; i++) {
                                s += wmin[k][i] * kmat[nt][i];
                            }
                            if (s > 0) {
                                if (s * min[k] < hmin) hmin = s * min[k];
                            }
                        }
                    }
                }
                xinv = hmin / 100;
                if (xinv <= 0) {
                    stop_simulation();
                }
            }

            ninv = 1;
            j = 0;

            for (int k=1; k<=nm; k++) {
                if (lmin[k] == 1) {
                    for (int kk=1; kk<=kinvariant; kk++) {
                        if (k == kinvar[kk]) {
                            j += 1;
                            for (int i=1; i<=nt-1; i++) {
                                t0[i][j] = wmin[k][i];
                            }
                            s = 0;
                            for (int i=1; i<=15; i++) {
                                s += wmin[k][i] * kmat[nt][i];
                            }
                            t0[nt][j] = s;
                        }
                    }
                }
            }
            nmin = j;

            for (int i=1; i<=nt-1; i++) {
                tt[i] = 0;
                for (int k=1; k<=nm; k++) {
                    for (int kk=1; kk<=kinvariant; kk++) {
                        if (k == kinvar[kk]) {
                            tt0[i] += min[k] * wmin[k][i];
                        }
                    }
                }
            }
            tt0[nt] = 0;
            for (int k=1; k<=nm; k++) {
                for (int kk=1; kk<=kinvariant; kk++) {
                    if (k == kinvar[kk]) {
                        for (int i=1; i<=15; i++) {
                            tt0[nt] += min[k] * wmin[k][i] * kmat[11][i];
                        }
                    }
                }
            }
            tt0[nt] -= xinv;

            for (int i=1; i<=nmin; i++) {
                for (int j=i; j<=nmin; j++) {
                    t[i][j] = 0;
                    for (int k=1; k<=nt; k++) {
                        t[i][j] += t0[k][i] * t0[k][j];
                        t[j][i] = t[i][j];
                    }
                }
            }
            for (int i=1; i<=nmin; i++) {
                tt[i] = 0;
                for (int k=1; k<=nt; k++) {
                    tt[i] += t0[k][i] * tt0[k];
                }
            }

            for (int k=2; k<=nmin; k++) {
                for (int i=k; i<=nmin; i++) {
                    if (abs(t[i][k-1]) > epsilon) {
                        u = t[k-1][k-1] / t[i][k-1];
                        for (int j=k; j<=nmin; j++) {
                            t[i][j] = t[k-1][j] - t[i][j] * u;
                        }
                        tt[i] = tt[k-1] - tt[i] * u;
                    }
                }
            }

            xx[nmin] = tt[nmin] / t[nmin][nmin];
            for (int i=nmin-1; i>=1; --i) {
                s = 0;
                for (int j=i+1; j<=nmin; j++) {
                    s += t[i][j] * xx[j];
                }
                xx[i] = (tt[i] - s) / t[i][i];
                if (abs(xx[i]) < epsilon) xx[i] = 0;
            }
            i = 0;
            for (int k=1; k<=kinvariant; k++) {
                i += 1;
                min[kinvar[k]] = xx[i];
                if (min[kinvar[k]] <= 0) ninv = 0;
            }
        }

        void peritec() {
            double u, s, swap;
            int nt, k, i, j, l, ni, nj, ii;
            string x_S;

            nt = ntot-1;
            double t[nt+1][nt+1] = {0}; double tt[nt+1] = {0};
            double t0[nt+1][nt+1] = {0}; double tt0[nt+1] = {0};
            double xx[nt+1] = {0};

            j = 0;
            for (int k=1; k<=nm; k++) {
                if (lmin[k] == 1) {
                    j += 1;
                    for (int i=1; i<=nt-1; i++) {
                        t0[i][j] = wmin[k][i];
                    }
                    s = 0;
                    for (int i=1; i<=ncm; i++) {
                        s += wmin[k][i] * kmat[11][i];
                    }
                    t0[11][j] = s;
                }
            }
            j += 1;
            t0[nt][j] = 2;
            nj = j;
            for (int i=1; i<=nt; i++) {
                tt0[i] = tot[i];
            }

            for (int i=1; i<=nj; i++) {
                for (int j=i; j<=nj; j++) {
                    t[i][j] = 0;
                    for (int k=1; k<=nt; k++) {
                        t[i][j] = t[i][j] + t0[k][i] * t0[k][j];
                        t[j][i] = t[i][j];
                    }
                }
            }
            for (int i=1; i<=nj; i++) {
                tt[i] = 0;
                for (int k=1; k<=nt; k++) {
                    tt[i] = tt[i] + t0[k][i] * tt0[k];
                }
            }

            for (int k=2; k<=nj; k++) {
                for (int i=k; i<=nj; i++) {
                    if (abs(t[i][k-1]) > epsilon) {
                        u = t[k-1][k-1] / t[i][k-1];
                        for (int j=k; j<=nj; j++) {
                            t[i][j] = t[k-1][j] - t[i][j] * u;
                        }
                        tt[i] = tt[k-1] - tt[i] * u;
                    }
                }
            }
            xx[nj] = tt[nj] / t[nj][nj];
            for (int i=nj-1; i>=1; --i) {
                s = 0;
                for (int j=i+1; j<=nj; j++) {
                    s += t[i][j] * xx[j];
                }
                xx[i] = (tt[i] - s) / t[i][i];
            }
            nperitec = 0;
            for (int i=1; i<=nj-1; i++) {
                if (xx[i] < 0) nperitec = 1;
            }
        }

        void compact() {
            int nrows = 0; int ncols = 0;
            double x;

            vector<string> col_names;
            vector<vector<string>> matrix;

            infile.open(min_file);
                getline(infile, line); 
                line = trim(line);
                stringstream linestream(line);

                while (linestream.good()) {
                    getline(linestream, value, ','); 
                    if (!value.empty()) {
                        col_names.push_back(value);
                        ncols++;
                    }
                }
                matrix.push_back(col_names);
                int include[ncols] = {0};

                while (infile.good()) {
                    vector<string> line_vect;
                    getline(infile, line); 
                    line = trim(line);

                    if (!line.empty()) {
                        stringstream linestream(line); 
                        for (int i=0; i<ncols; i++) {
                            getline(linestream, value, ',');
                            line_vect.push_back(value);

                            if (i > 0) {
                                x = stod(value);
                                if (x > 0) {
                                    include[i] = 1;
                                }
                            }

                        }
                        matrix.push_back(line_vect);
                        nrows++;
                    }
                }
            infile.close();

            outfile.open(min_file, ios::trunc);
                for (int i=0; i<=nrows; i++) {
                    for (int k=0; k<ncols; k++) {
                        if (k == 0) {
                            outfile << matrix[i][k];
                        } else {
                            if (include[k] == 1) {
                                outfile << "," << matrix[i][k];
                            }
                        }
                    }
                    outfile << endl;
                }
            outfile.close();
        }

        void stop_simulation() {
            if (output == 1) {
                compact();
            }
            exit(0);
        }

};

int main() {
    Simulation test;

    return 0;
}
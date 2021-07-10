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

double g0(double x) {
    return 2 * (1 - (1-x) * exp(-x)) / pow(x, 2);
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
    double mwev, mwev0, mh2o, psc3;
    int nc, na, ncm, nm, nminer, nbmin, kinvariant, ncpt, ncomp;
    int output_step, output, print_step, verbose;
    int ic, ia;
    string a_S, x_S, zone_S, q0_S, constituant_S, inc_S;
    string tra_file, log_file, chem_file, event_file, min_file;
    string system, units, miner_S;
    string line, value;
    ifstream infile; ofstream outfile; ofstream log;

    // declare arrays and vectors
    array<string,n+1> aq_S;
    double tot[ntot+1] = {0}; double totinit[ntot+1] = {0};
    double tot0[ntot+1] = {0}; double totest[ntot+1];
    double psc[ncomplex+1]; double gact[n+1] = {0}; 
    double gact0[n+1] = {0}; double gact1[n+1] = {0}; 
    double mol[n+1] = {0}; double mol0[n+1] = {0}; 
    double mol1[n+1] = {0}; double molal[n+1] = {0}; 
    double molal0[n+1] = {0}; double act[n+1] = {0}; 
    double act0[n+1] = {0};  double atom[n+1] = {0}; 
    int kmat[n+1][n+1] = {0}; int nch[n+1] = {0};


    // Declare vectors
    vector<int> ica, kinvar, linvar;
    vector<vector<double>> wmin;
    vector<double> mu, mum, psol, psol0, pai, pai0;
    vector<double> min, min0, minp, minp0, c_ion;
    vector<int> lmin, lmin0, lmin1;
    vector<string> nom_ion_S, mineral_S;




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
            mh2o = 55.51;

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
            infile.open("aquv.dat");
            for (int i=0; i<=n; i++) {
                getline(infile, line); stringstream linestream(line);
                getline(linestream, value, ','); aq_S[i] = value;
                getline(linestream, value, ','); atom[i] = stod(value);
                getline(linestream, value, ','); nch[i] = stoi(value);   
            }
            infile.close();

            // read coefficients of partial derivatives of EVP equation set
            infile.open("matrice2");
            for (int i=1; i<=n; i++) {
                getline(infile, line); stringstream linestream(line);
                for (int j=1; j<=n; j++) {
                    getline(linestream, value, ','); kmat[i][j] = stoi(value);
                }
            }
            infile.close();

            // open the log file in append mode
            if (print_step > 0 and output != 0) {
                log.open(log_file, ios::app);
            }

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

            infile.open("complex3");
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
            psc3 = psc[3];
            psc[3] = psc[3] * psc[14] * pow(10, po);
            infile.close();
            
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
                        vector<double> temp;
                        for (int k=0; k<=ncm; k++) temp.push_back(0);
                        wmin.push_back(temp);
                }
                for (int i=0; i<=ncm; i++) {
                    mu.push_back(0); nom_ion_S.push_back("");
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
                    getline(linestream, value, ','); a_S = upper_to_lower(value);
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
                    getline(linestream, value, ','); mineral_S[k] = value;
                    getline(linestream, value, ','); ncomp = stoi(value);

                    for (int i=1; i<=ncomp; i++) {
                        getline(linestream, value, ','); c_ion[i] = stod(value);
                        getline(linestream, value, ','); nom_ion_S[i] = value;
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
                    mum[k] = at + bt / 300 * temp + ct / 3000 * pow(temp, 2) + dt / 3000000 * pow(temp, 3) + et / 300000000 * pow(temp, 4);
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






        }

        void read_transfer() {
            // read stockage
            infile.open("stockage");
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

        // void eql_actp() {
        //     double c[10]; double a[12]; double h[4];
        //     double at, bt, ct, dt, et;
        //     double u, z, w, v, f, fj, co, s;
        //     int nc, na, nn;

        //     c[1] = molal[1]; c[2] = molal[2]; c[3] = molal[3];
        //     c[4] = molal[4]; c[5] = molal[22]; c[6] = molal[5];
        //     c[7] = molal[18]; c[8] = molal[23]; c[9] = molal[11];
        //     a[1] = molal[6]; a[2] = molal[7]; a[3] = molal[25];
        //     a[4] = molal[12]; a[5] = molal[14]; a[6] = molal[13];
        //     a[7] = molal[24]; a[8] = molal[19]; a[9] = molal[20];
        //     a[10] = molal[21]; a[11] = molal[8];
        //     h[1] = molal[10]; h[2] = molal[9]; h[3] = molal[15];

        //     {
        //     file.open("coefft4");
        //     getline(file, line);
        //     stringstream linestream(line);
        //     string value;
        //     getline(linestream, value, ','); nc = stoi(value);
        //     getline(linestream, value, ','); na = stoi(value);
        //     getline(linestream, value, ','); nn = stoi(value);
            

        //     if (ndepact == 0) {
        //         // initialize 1d, 2d, and 3d vectors to 0
        //         for (int i=0; i<=nc; i++) nzc.push_back(0);
        //         for (int i=0; i<=na; i++) nza.push_back(0);
        //         for (int i=0; i<=nc; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=na; k++) temp.push_back(0); b0.push_back(temp);
        //         }
        //         for (int i=0; i<=nc; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=na; k++) temp.push_back(0); b1.push_back(temp);
        //         }
        //         for (int i=0; i<=nc; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=na; k++) temp.push_back(0); b2.push_back(temp);
        //         }
        //         for (int i=0; i<=nc; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=na; k++) temp.push_back(0); c0.push_back(temp);
        //         }
        //         for (int i=0; i<=nc; i++) {
        //             vector<vector<double>> temp2d;
        //             for (int k=0; k<=nc; k++) {
        //                 vector<double> temp1d;
        //                 for (int j=0; j<=na; j++) temp1d.push_back(0); temp2d.push_back(temp1d);
        //             }
        //             scp.push_back(temp2d);
        //         }
        //         for (int i=0; i<=na; i++) {
        //             vector<vector<double>> temp2d;
        //             for (int k=0; k<=na; k++) {
        //                 vector<double> temp1d;
        //                 for (int j=0; j<=nc; j++) temp1d.push_back(0); temp2d.push_back(temp1d);
        //             }
        //             sap.push_back(temp2d);
        //         }
        //         for (int i=0; i<=nc; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=nc; k++) temp.push_back(0); tcp.push_back(temp);
        //         }
        //         for (int i=0; i<=na; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=na; k++) temp.push_back(0); tap.push_back(temp);
        //         }
        //         for (int i=0; i<=nn; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=nc; k++) temp.push_back(0); lcp.push_back(temp);
        //         }
        //         for (int i=0; i<=nn; i++) {
        //             vector<double> temp;
        //             for (int k=0; k<=na; k++) temp.push_back(0); lap.push_back(temp);
        //         }
        //         for (int i=0; i<=nn; i++) {
        //             vector<vector<double>> temp2d;
        //             for (int k=0; k<=nc; k++) {
        //                 vector<double> temp1d;
        //                 for (int j=0; j<=na; j++) temp1d.push_back(0); temp2d.push_back(temp1d);
        //             }
        //             xip.push_back(temp2d);
        //         }

        //         // Read data into matrices/vectors
        //         for (int i=1; i<=nc; i++) {
        //             getline(file, line); stringstream linestream(line); string value;

        //             getline(linestream, value, ','); string x_S = value;
        //             getline(linestream, value, ','); nzc[i] = stoi(value);
        //         }
        //         for (int i=1; i<=na; i++) {
        //             getline(file, line); stringstream linestream(line); string value;

        //             getline(linestream, value, ','); string x_S = value;
        //             getline(linestream, value, ','); nza[i] = stoi(value);
        //         }
        //         {
        //             getline(file, line); stringstream linestream(line); string value;

        //             getline(linestream, value, ','); string x_S = value;
        //             getline(linestream, value, ','); at = stod(value);
        //             getline(linestream, value, ','); bt = stod(value);
        //             getline(linestream, value, ','); ct = stod(value);
        //             getline(linestream, value, ','); dt = stod(value);
        //             getline(linestream, value, ','); et = stod(value);

        //             ap0 = temperature(at, bt, ct, dt, et, temp);
        //         }

        //         for (int i=1; i<=nc; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 {
        //                     getline(file, line); stringstream linestream(line); string value;

        //                     getline(linestream, value, ','); string x_S = value;
        //                     getline(linestream, value, ','); at = stod(value);
        //                     getline(linestream, value, ','); bt = stod(value);
        //                     getline(linestream, value, ','); ct = stod(value);
        //                     getline(linestream, value, ','); dt = stod(value);
        //                     getline(linestream, value, ','); et = stod(value);

        //                     b0[i][j] = temperature(at, bt, ct, dt, et, temp);
        //                 }
        //                 {
        //                     getline(file, line); stringstream linestream(line); string value;

        //                     getline(linestream, value, ','); string x_S = value;
        //                     getline(linestream, value, ','); at = stod(value);
        //                     getline(linestream, value, ','); bt = stod(value);
        //                     getline(linestream, value, ','); ct = stod(value);
        //                     getline(linestream, value, ','); dt = stod(value);
        //                     getline(linestream, value, ','); et = stod(value);

        //                     b1[i][j] = temperature(at, bt, ct, dt, et, temp);
        //                 }
        //                 {
        //                     getline(file, line); stringstream linestream(line); string value;

        //                     getline(linestream, value, ','); string x_S = value;
        //                     getline(linestream, value, ','); at = stod(value);
        //                     getline(linestream, value, ','); bt = stod(value);
        //                     getline(linestream, value, ','); ct = stod(value);
        //                     getline(linestream, value, ','); dt = stod(value);
        //                     getline(linestream, value, ','); et = stod(value);

        //                     b2[i][j] = temperature(at, bt, ct, dt, et, temp);
        //                 }
        //                 {
        //                     getline(file, line); stringstream linestream(line); string value;

        //                     getline(linestream, value, ','); string x_S = value;
        //                     getline(linestream, value, ','); at = stod(value);
        //                     getline(linestream, value, ','); bt = stod(value);
        //                     getline(linestream, value, ','); ct = stod(value);
        //                     getline(linestream, value, ','); dt = stod(value);
        //                     getline(linestream, value, ','); et = stod(value);

        //                     c0[i][j] = temperature(at, bt, ct, dt, et, temp);
        //                 }
        //             }
        //         }

        //         for (int i=1; i<=nc-1; i++) {
        //             for (int j=i+1; j<=nc; j++) {
        //                 getline(file, line); stringstream linestream(line); string value;

        //                 getline(linestream, value, ','); string x_S = value;
        //                 getline(linestream, value, ','); at = stod(value);
        //                 getline(linestream, value, ','); bt = stod(value);
        //                 getline(linestream, value, ','); ct = stod(value);
        //                 getline(linestream, value, ','); dt = stod(value);
        //                 getline(linestream, value, ','); et = stod(value);

        //                 tcp[i][j] = temperature(at, bt, ct, dt, et, temp);
        //                 tcp[j][i] = tcp[i][j];
        //             }
        //         }

        //         for (int i=1; i<=na-1; i++) {
        //             for (int j=i+1; j<=na; j++) {
        //                 getline(file, line); stringstream linestream(line); string value;

        //                 getline(linestream, value, ','); string x_S = value;
        //                 getline(linestream, value, ','); at = stod(value);
        //                 getline(linestream, value, ','); bt = stod(value);
        //                 getline(linestream, value, ','); ct = stod(value);
        //                 getline(linestream, value, ','); dt = stod(value);
        //                 getline(linestream, value, ','); et = stod(value);

        //                 tap[i][j] = temperature(at, bt, ct, dt, et, temp);
        //                 tap[j][i] = tap[i][j];
        //             }
        //         }

        //         for (int k=1; k<=nc-1; k++) {
        //             for (int i=k+1; i<=nc; i++) {
        //                 for (int j=1; j<=na; j++) {
        //                     getline(file, line); stringstream linestream(line); string value;

        //                     getline(linestream, value, ','); string x_S = value;
        //                     getline(linestream, value, ','); at = stod(value);
        //                     getline(linestream, value, ','); bt = stod(value);
        //                     getline(linestream, value, ','); ct = stod(value);
        //                     getline(linestream, value, ','); dt = stod(value);
        //                     getline(linestream, value, ','); et = stod(value);

        //                     scp[k][i][j] = temperature(at, bt, ct, dt, et, temp);
        //                     scp[i][k][j] = scp[k][i][j];
        //                 }
        //             }
        //         }

        //         for (int k=1; k<=na-1; k++) {
        //             for (int i=k+1; i<=na; i++) {
        //                 for (int j=1; j<=nc; j++) {
        //                     getline(file, line); stringstream linestream(line); string value;

        //                     getline(linestream, value, ','); string x_S = value;
        //                     getline(linestream, value, ','); at = stod(value);
        //                     getline(linestream, value, ','); bt = stod(value);
        //                     getline(linestream, value, ','); ct = stod(value);
        //                     getline(linestream, value, ','); dt = stod(value);
        //                     getline(linestream, value, ','); et = stod(value);

        //                     sap[k][i][j] = temperature(at, bt, ct, dt, et, temp);
        //                     sap[i][k][j] = sap[k][i][j];
        //                 }
        //             }
        //         }

        //         for (int i=1; i<=nn; i++) {
        //             for (int j=1; j<=nc; j++) {
        //                 getline(file, line); stringstream linestream(line); string value;

        //                 getline(linestream, value, ','); string x_S = value;
        //                 getline(linestream, value, ','); at = stod(value);
        //                 getline(linestream, value, ','); bt = stod(value);
        //                 getline(linestream, value, ','); ct = stod(value);
        //                 getline(linestream, value, ','); dt = stod(value);
        //                 getline(linestream, value, ','); et = stod(value);

        //                 lcp[i][j] = temperature(at, bt, ct, dt, et, temp);
        //             }
        //         }

        //         for (int i=1; i<=nn; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 getline(file, line); stringstream linestream(line); string value;

        //                 getline(linestream, value, ','); string x_S = value;
        //                 getline(linestream, value, ','); at = stod(value);
        //                 getline(linestream, value, ','); bt = stod(value);
        //                 getline(linestream, value, ','); ct = stod(value);
        //                 getline(linestream, value, ','); dt = stod(value);
        //                 getline(linestream, value, ','); et = stod(value);

        //                 lap[i][j] = temperature(at, bt, ct, dt, et, temp);
        //             }
        //         }

        //         for (int k=1; k<=nn; k++) {
        //             for (int i=1; i<=nc; i++) {
        //                 for (int j=1; j<=na; j++) {
        //                     xip[k][i][j] = 0;
        //                 }
        //             }
        //         }
        //         xip[2][9][1] = -0.0102;
        //         xip[2][1][2] = 0.046;
        //     }
            
        //     file.close();
        //     }

        //     // declare and initialize temporary matrices to 0
        //     double ec[(nc+1)][(nc+1)] = { 0 }; double ea[(na+1)][(na+1)] = { 0 };
        //     double fc[(nc+1)][(nc+1)] = { 0 }; double fa[(na+1)][(na+1)] = { 0 };
        //     double xc[(nc+1)][(nc+1)] = { 0 }; double xa[(na+1)][(na+1)] = { 0 };
        //     double pp[(nc+1)][(nc+1)] = { 0 }; double qp[(na+1)][(na+1)] = { 0 };
        //     double p[(nc+1)][(nc+1)] = { 0 }; double q[(na+1)][(na+1)] = { 0 };
        //     double pf[(nc+1)][(nc+1)] = { 0 }; double qf[(na+1)][(na+1)] = { 0 };
        //     double cc[(nc+1)][(na+1)] = { 0 }; double bf[(nc+1)][(na+1)] = { 0 };
        //     double b[(nc+1)][(na+1)] = { 0 }; double bp[(nc+1)][(na+1)] = { 0 };
        //     double gc[(nc+1)] = { 0 }; double ga[(na+1)] = { 0 }; double gn[(nn+1)] = { 0 };

        //     bp0 = 1.2; mh2o = 55.51;

        //     u = 0; z = 0;
        //     for (int i=1; i<=nc; i++) {
        //         u += c[i] * pow(nzc[i], 2); z += c[i] * nzc[i];
        //     }
        //     for (int j=1; j<=na; j++) {
        //         u += a[j] * pow(nza[j], 2); z += a[j] * nza[j];
        //     }
        //     fi = u / 2; fj = sqrt(fi);
        //     u = 6 * ap0 * fj;
        //     for (int i=1; i<=nc-1; i++) {
        //         for (int j=i+1; j<=nc; j++) {
        //             if (nzc[i] == nzc[j]) {
        //                 ec[i][j] = 0; fc[i][j] = 0;
        //             } else {
        //                 xc[i][j] = 2 * u;
        //                 xc[i][i] = pow(nzc[i], 2) * u;
        //                 xc[j][j] = pow(nzc[j], 2) * u;
        //                 ec[i][j] = (j0(xc[i][j]) - j0(xc[i][i]) / 2 - j0(xc[j][j]) / 2) / fi / 2;
        //                 fc[i][j] = ((xc[i][j] * j1(xc[i][j]) - xc[i][i] * j1(xc[i][i]) / 2 - 
        //                             xc[j][j] * j1(xc[j][j]) / 2) / pow(fi, 2) / 4 - ec[i][j] / fi);
        //                 ec[j][i] = ec[i][j]; 
        //                 fc[j][i] = fc[i][j];
        //             }
        //         }
        //     }
        //     for (int i=1; i<=na-1; i++) {
        //         for (int j=i+1; j<=na; j++) {
        //             if (nza[i] == nza[j]) {
        //                 ea[i][j] = 0; fa[i][j] = 0;
        //             } else {
        //                 xa[i][j] = 2 * u; 
        //                 xa[i][i] = pow(nza[i], 2) * u;
        //                 xa[j][j] = pow(nza[j], 2) * u;
        //                 ea[i][j] = (j0(xa[i][j]) - j0(xa[i][i]) / 2 - j0(xa[j][j]) / 2) / fi / 2;
        //                 fa[i][j] = ((xa[i][j] * j1(xa[i][j]) - xa[i][i] * j1(xa[i][i]) / 
        //                                 2 - xa[j][j] * j1(xa[j][j]) / 2) / pow(fi, 2) / 4 - ea[i][j] / fi);
        //                 ea[j][i] = ea[i][j];
        //                 fa[j][i] = fa[i][j];
        //             }
        //         }
        //     }
        //     for (int i=1; i<=nc-1; i++) {
        //         for (int j=i+1; j<=nc; j++) {
        //             pp[i][j] = fc[i][j]; 
        //             p[i][j] = tcp[i][j] + ec[i][j];
        //             pf[i][j] = p[i][j] + pp[i][j] * fi;
        //             pp[j][i] = pp[i][j]; 
        //             p[j][i] = p[i][j]; 
        //             pf[j][i] = pf[i][j];
        //         }
        //     }
        //     for (int i=1; i<=na-1; i++) {
        //         for (int j=i+1; j<=na; j++) {
        //             qp[i][j] = fa[i][j]; 
        //             q[i][j] = tap[i][j] + ea[i][j];
        //             qf[i][j] = q[i][j] + qp[i][j] * fi;
        //             qp[j][i] = qp[i][j]; 
        //             q[j][i] = q[i][j]; 
        //             qf[j][i] = qf[i][j];
        //         }
        //     }
        //     w = fj * 12;
        //     for (int i=1; i<=nc; i++) {
        //         for (int j=1; j<=na; j++) {
        //             cc[i][j] = c0[i][j] / sqrt(nzc[i] * nza[j]) / 2;
        //             if (nzc[i] == 2 and nza[j] == 2) v = fj * 1.4e0;
        //             if (nzc[i] == 1 or nza[j] == 1) v = fj * 2;
        //             bf[i][j] = (b0[i][j] + b1[i][j] * exp(-v) + b2[i][j] * exp(-w));
        //             b[i][j] = (b0[i][j] + b1[i][j] * (g0(v)) + b2[i][j] * (g0(w)));
        //             bp[i][j] = (b1[i][j] * (g1(v)) / fi + b2[i][j] * (g1(w)) / fi);
        //         }
        //     }
        //     f = -ap0 * (fj / (1 + bp0 * fj) + 2 / bp0 * log(1 + bp0 * fj));
        //     for (int i=1; i<=nc; i++) {
        //         for (int j=1; j<=na; j++) {
        //             f += c[i] * a[j] * bp[i][j];
        //         }
        //     }
        //     for (int i=1; i<=(nc-1); i++) {
        //         for (int j=(i+1); j<=nc; j++) {
        //             f += c[i] * c[j] * pp[i][j];
        //         }
        //     }
        //     for (int i=1; i<=(na-1); i++) {
        //         for (int j=(i+1); j<=na; j++) {
        //             f += a[i] * a[j] * qp[i][j];
        //         }
        //     }
        //     for (int ii=1; ii<=nc; ii++) {
        //         u = pow(nzc[ii], 2) * f;
        //         for (int j=1; j<=na; j++) {
        //             u += a[j] * (b[ii][j] * 2 + z * cc[ii][j]);
        //         } 
        //         for (int i=1; i<=nc; i++) {
        //             if (i != ii) {
        //                 v = 0;
        //                 for (int j=1; j<=na; j++) {
        //                     v += a[j] * scp[ii][i][j];
        //                 }
        //                 u += c[i] * (p[ii][i] * 2 + v);
        //             }
        //         }
        //         for (int i=1; i<=(na-1); i++) {
        //             for (int j=(i+1); j<=na; j++) {
        //                 u += a[i] * a[j] * sap[i][j][ii];
        //             }
        //         }
        //         for (int i=1; i<=nc; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 u += c[i] * a[j] * cc[i][j] * nzc[ii];
        //             }
        //         }
        //         for (int i=1; i<=nn; i++) {
        //             u += h[i] * lcp[i][ii] * 2;
        //         }
        //         for (int k=1; k<=nn; k++) {
        //             for (int j=1; j<=na; j++) {
        //                 u += h[k] * a[j] * xip[k][ii][j];
        //             }
        //         }
        //         gc[ii] = exp(u);
        //     }
        //     for (int jj=1; jj<=na; jj++) {
        //         u = pow(nza[jj], 2) * f;
        //         for (int i=1; i<=nc; i++) {
        //             u += c[i] * (b[i][jj] * 2 + z * cc[i][jj]);
        //         }
        //         for (int i=1; i<=na; i++) {
        //             if (i != jj) {
        //                 v = 0;
        //                 for (int j=1; j<=nc; j++) {
        //                     v += c[j] * sap[jj][i][j];
        //                 }
        //                 u += a[i] * (q[jj][i] * 2 + v);
        //             }
        //         }
        //         for (int i=1; i<=nc-1; i++) {
        //             for (int j=i+1; j<=nc; j++) {
        //                 u += c[i] * c[j] * scp[i][j][jj];
        //             }
        //         }
        //         for (int i=1; i<=nc; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 u += c[i] * a[j] * cc[i][j] * nza[jj];
        //             }
        //         }
        //         for (int j=1; j<=nn; j++) {
        //             u += h[j] * lap[j][jj];
        //         }
        //         for (int k=1; k<=nn; k++) {
        //             for (int i=1; i<=nc; i++) {
        //                 u += h[k] * c[i] * xip[k][i][jj];
        //             }
        //         }
        //         ga[jj] = exp(u);
        //     }
        //     for (int k=1; k<=nn; k++) {
        //         u = 0;
        //         for (int i=1; i<=nc; i++) {
        //             u += c[i] * lcp[k][i] * 2;
        //         }
        //         for (int j=1; j<=na; j++) {
        //             u += a[j] * lap[k][j] * 2;
        //         }
        //         for (int i=1; i<=nc; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 u += c[i] * a[j] * xip[k][i][j];
        //             }
        //         }
        //         gn[k] = exp(u);
        //     }
        //     u = -ap0 * pow(fi, 1.5e0) / (1 + bp0 * fj);
        //     for (int i=1; i<=nc; i++) {
        //         for (int j=1; j<=na; j++) {
        //             u += c[i] * a[j] * (bf[i][j] + z * cc[i][j]);
        //         }
        //     }
        //     for (int i=1; i<=nc-1; i++) {
        //         for (int j=i+1; j<=nc; j++) {
        //             v = 0;
        //             for (int k=1; k<=na; k++) {
        //                 v += a[k] * scp[i][j][k];
        //             }
        //             u += c[i] * c[j] * (pf[i][j] + v);
        //         }
        //     }
        //     for (int i=1; i<=na-1; i++) {
        //         for (int j=i+1; j<=na; j++) {
        //             v = 0;
        //             for (int k=1; k<=nc; k++) {
        //                 v += c[k] * sap[i][j][k];
        //             }
        //             u += a[i] * a[j] * (qf[i][j] + v);
        //         }
        //     }
        //     for (int k=1; k<=nn; k++) {
        //         for (int i=1; i<=nc; i++) {
        //             u += h[k] * c[i] * lcp[k][i];
        //         }
        //     }
        //     for (int k=1; k<=nn; k++) {
        //         for (int j=1; j<=na; j++) {
        //             u += h[k] * a[j] * lap[k][j];
        //         }
        //     }
        //     for (int k=1; k<=nn; k++) {
        //         for (int i=1; i<=nc; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 u += h[k] * c[i] * a[j] * xip[k][i][j];
        //             }
        //         }
        //     }

        //     s = 0;
        //     for (int i=1; i<=nc; i++) s += c[i];
        //     for (int j=1; j<=na; j++) s += a[j];

        //     co = 1 + 2 * u / s;
        //     aw = exp(-s * co / mh2o);
        //     gact[1] = gc[1]; gact[2] = gc[2]; gact[3] = gc[3];
        //     gact[4] = gc[4]; gact[22] = gc[5]; gact[5] = gc[6];
        //     gact[18] = gc[7]; gact[23] = gc[8]; gact[11] = gc[9];
        //     gact[6] = ga[1]; gact[7] = ga[2]; gact[25] = ga[3];
        //     gact[12] = ga[4]; gact[14] = ga[5]; gact[13] = ga[6];
        //     gact[24] = ga[7]; gact[19] = ga[8]; gact[20] = ga[9];
        //     gact[21] = ga[10]; gact[8] = ga[11];
        //     gact[10] = aw * aw * pow(gn[1], log(10));
        //     gact[9] = gn[2]; gact[15] = gn[3];
        //     gact[16] = 1; gact[17] = 1;
        //     ndepact = 1;
        // }

        // void eql_density() {
        //     int nc, na;

        //     nc = 5; na = 5;

        //     double s[nc+1][na+1] = {0};
        //     double ao[nc+1][na+1] = {0}; double bo[nc+1][na+1] = {0};
        //     double au[nc+1][na+1] = {0}; double bu[nc+1][na+1] = {0};

        //     file.open("densite");
        //     for (int i=1; i<=5; i++) {
        //         for (int j=1; j<=5; j++) {
        //             getline(file, line);
        //             stringstream linestream(line);
        //             string value;

        //             getline(linestream, value, ','); x_S = value;
        //             getline(linestream, value, ','); ao[i][j] = stod(value);
        //             getline(linestream, value, ','); bo[i][j] = stod(value);
        //             getline(linestream, value, ','); au[i][j] = stod(value);
        //             getline(linestream, value, ','); bu[i][j] = stod(value);
        //         }
        //     }
        //     file.close();

        //     if (units == "molar") {
        //         dens = 1; u = 0;
        //         for (int j=1; j<=na; j++) u += ani[j];
        //         for (int i=1; i<=nc; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 s[i][j] = (nchcat[i] + nchani[j]) / 2 * cat[i] * ani[j] / nchcat[i] / nchani[j] / u;
        //                 dens += ao[i][j] * s[i][j] + bo[i][j] * pow(s[i][j], 2);
        //             }
        //         }
        //     } else if (units == "molal") {
        //         dens = 1; u = 0;
        //         for (int j=1; j<=na; j++) u += ani[j] * nchani[j];
        //         for (int i=1; i<=nc; i++) {
        //             for (int j=1; j<=na; j++) {
        //                 s[i][j] = (nchcat[i] + nchani[j]) / 2 * cat[i] * ani[j] / u;
        //                 dens += au[i][j] * s[i][j] + bu[i][j] * pow(s[i][j], 2);
        //             }
        //         }
        //     }
        // }

        // void eql_invar() {
        //     vector<string> minv_S, minvar_S;
        //     int n1, n2, n3, n4, ncond, nbmin, ninvar, det, det1, z, ncm, i1, i2;
        //     double swap, ah2o;
        //     int ii, jj, kk;

        //     kinvariant = 0; ncm = 14;
        //     nbmin = 10;
        //     ninvar = nbmin + 3;

        //     for (int i=0; i<=ninvar; i++) minv_S.push_back("");
        //     for (int i=0; i<=ninvar; i++) minvar_S.push_back("");
        //     int kinv[ninvar+1] = {0}; double psminv[ninvar+1] = {0};
        //     double winv[ninvar+1][ncm+1] = {0};
        //     double psminvar[ninvar+1] = {0}; 
        //     double t0[ninvar+1][ninvar+1] = {0}; double t1[ninvar+1][ninvar+1] = {0};
        //     double t2[ninvar+1][ninvar+1] = {0}; double t3[ninvar+1][ncm+1] = {0};
        //     double t4[ncm+1][ncm+1] = {0}; double tt4[ncm+1] = {0};

        //     for (int k=1; k<=3; k++) {
        //         psminv[k] = log10(psc[k]);
        //     }
        //     winv[1][11] = 1; winv[1][13] = 1; winv[1][0] = -1;
        //     winv[2][11] = 1; winv[2][12] = -1; winv[2][14] = 1;
        //     winv[3][11] = 1; winv[3][12] = 1; winv[3][0] = -1;

        //     n1 = 3;
        //     for (int l=1; k<=nm; k++) {
        //         if (lmin[k] == 1) {
        //             n1 += 1;
        //             kinv[n1] = k;
        //             minv_S[n1] = mineral_S[k];
        //             psminv[n1] = log10(psol[k]);
        //             for (int j=0; j<=ncm; j++) {
        //                 winv[n1][j] = wmin[k][j];
        //             }
        //         }
        //     }
        //     for (int i=1; i<=n1; i++) {
        //         swap = winv[i][0];
        //         winv[i][0] = winv[i][14];
        //         winv[i][14] = swap;
        //     }
        //     for (int i=1; i<=n1; i++) {
        //         for (int j=i; j<=n1; j++) {
        //             t1[i][j] = 0;
        //             for (int k=0; k<=ncm-1; k++) {
        //                 t1[i][j] += winv[i][k] * winv[j][k];
        //                 t1[j][i] = t1[i][j];
        //                 t0[i][j] = t1[i][j];
        //                 t0[j][i] = t0[i][j];
        //             }
        //         }
        //     }
        //     for (int k=2; k<=n1; k++) {
        //         for (int i=k; i<=n1; i++) {
        //             if (abs(t1[i][k-1]) > epsilon) {
        //                 u = t1[k-1][k-1] / t1[i][k-1];
        //                 for (int j=k; j<=n1; j++) {
        //                     t1[i][j] = t1[k-1][j] - t1[i][j] * u;
        //                     if (abs(t1[i][j]) < epsilon) t1[i][j] = 0;
        //                 }
        //             }
        //         }
        //     }
        //     det = 1;
        //     for (int i=1; i<=n1; i++) {
        //         if (abs(t1[i][i]) < epsilon) {
        //             det = 0;
        //             break;
        //         }
        //     }

        //     if (det == 0) {
        //         n3 = 0;
        //         n2 = n1-1;
        //         for (int kk=1; kk<=n1; kk++) {
        //             ii = 0;
        //             for (int i=1; i<=n1; i++) {
        //                 if (i != kk) {
        //                     ii += 1;
        //                     jj = 0;
        //                     for (int j=1; j<=n1; j++) {
        //                         if (j != kk) {
        //                             jj += 1;
        //                             t2[ii][jj] = t0[i][j];
        //                         }
        //                     }
        //                 }
        //             }

        //             for (int k=2; k<=n2; k++) {
        //                 for (int i=k; i<=n2; i++) {
        //                     if (abs(t2[i][k-1]) > epsilon) {
        //                         u = t2[k-1][k-1] / t2[i][k-1];
        //                         for (int j=k; j<=n2; j++) {
        //                             t2[i][j] = t2[k-1][j] - t2[i][j] * u;
        //                             if (abs(t2[i][j]) < epsilon) t2[i][j] = 0;
        //                         }
        //                     }
        //                 }
        //             }

        //             det1 = 1;
        //             for (int i=1; i<=n2; i++) {
        //                 if (abs(t2[i][i]) < epsilon) {
        //                     det1 = 0;
        //                     break;
        //                 }
        //             }
        //             if (det1 == 1) {
        //                 n3 += 1;
        //                 kinvar[n3] = kinv[kk];
        //                 minvar_S[n3] = minv_S[kk];
        //                 psminvar[n3] = psminv[kk];
        //                 for (int j=0; j<=ncm; j++) {
        //                     t3[n3][j] = winv[kk][j];
        //                 }
        //             }
        //         }

        //         if (n3 == 0) {
        //             kinvariant = -1;
        //         } else if (n3 > 0) {
        //             n4 = ncm;
        //             for (int j=ncm; j>=1; --j) {
        //                 u = 0;
        //                 for (int i=1; i<=n3; i++) {
        //                     u += pow(t3[i][j], 2);
        //                 }
        //                 if (u < epsilon) {
        //                     for (int k=j+1; k<=n4; k++) {
        //                         for (int i=1; i<=n3; i++) {
        //                             t3[i][k-1] = t3[i][k];
        //                         }
        //                     }
        //                     n4 -= 1;
        //                 }
        //             }

        //             for (int i=1; i<=n4; i++) {
        //                 for (int j=1; j<=n4; j++) {
        //                     t4[i][j] = 0;
        //                     for (int k=1; k<=n3; k++) {
        //                         t4[i][j] += t3[k][i] * t3[k][j];
        //                         t4[j][i] = t4[i][j];
        //                     }
        //                 }
        //             }
        //             for (int i=1; i<=n4; i++) {
        //                 tt4[i] = 0;
        //                 for (int k=1; k<=n3; k++) {
        //                     tt4[i] += t3[k][i] * psminvar[k];
        //                 }
        //             }

        //             for (int k=2; k<=n4; k++) {
        //                 for (int i=k; i<=n4; i++) {
        //                     if (abs(t4[i][k-1]) > epsilon) {
        //                         u = t4[k-1][k-1] / t4[i][k-1];
        //                         for (int j=k; j<=n4; j++) {
        //                             t4[i][j] = t4[k-1][j] - t4[i][j] * u;
        //                             if (abs(t4[i][j]) < epsilon) {
        //                                 t4[i][j] = 0;
        //                             }
        //                         }
        //                         tt4[i] = tt4[k-1] - tt4[i] * u;
        //                     }
        //                 }
        //             }

        //             if (abs(t4[n4][n4]) > epsilon) {
        //                 ah2o = pow(10, (tt4[n4] / t4[n4][n4]));
        //                 if (ah2o > 1 or ah2o <= 0.02) {
        //                     kinvariant = -2;
        //                 } else {
        //                     kinvariant = n3;
        //                     for (int i=kinvariant; i>=1; --i) {
        //                         if (kinvar[i] == 0) {
        //                             for (int k=1; k<=kinvariant-1; k++) {
        //                                 kinvar[k] = kinvar[k+1];
        //                             }
        //                             kinvariant -= 1;
        //                         }
        //                     }
        //                 }
        //             } else if(abs(t4[n4][n4]) <= epsilon) {
        //                 kinvariant = -2;
        //             }
        //         }
        //     }
        // }
};

int main() {
    Simulation test;

    return 0;
}
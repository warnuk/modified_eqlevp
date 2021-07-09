#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include <new>

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

void actp(double *molal, double *gact, double &aw, double &fi, double &temp, double &ap0, double &bp0, double &mh2o) {

    // Declare arrays and variables
    double c[10]; double a[12]; double h[4];
    double at, bt, ct, dt, et;
    double u, z, w, v, f, s, co, fj;
    int nc, na, nn;

    ifstream file; string line;

    c[1] = molal[1]; c[2] = molal[2]; c[3] = molal[3];
    c[4] = molal[4]; c[5] = molal[22]; c[6] = molal[5];
    c[7] = molal[18]; c[8] = molal[23]; c[9] = molal[11];
    a[1] = molal[6]; a[2] = molal[7]; a[3] = molal[25];
    a[4] = molal[12]; a[5] = molal[14]; a[6] = molal[13];
    a[7] = molal[24]; a[8] = molal[19]; a[9] = molal[20];
    a[10] = molal[21]; a[11] = molal[8];
    h[1] = molal[10]; h[2] = molal[9]; h[3] = molal[15];

    // Open coefft4 file
    file.open("coefft4");

    // Read first line for nc, na, & nn
    {   
        getline(file, line);
        stringstream linestream(line);
        string value;
        getline(linestream, value, ',');
        nc = stoi(value);
        getline(linestream, value, ',');
        na = stoi(value);
        getline(linestream, value, ',');
        nn = stoi(value);
    }

    // Initialize 0 arrays for the rest of the data
    int nzc[(nc+1)] = { 0 }; int nza[(na+1)] = { 0 };
    double b0[(nc+1)][(na+1)] = { 0 }; double b1[(nc+1)][(na+1)] = { 0 };
    double b2[(nc+1)][(na+1)] = { 0 }; double c0[(nc+1)][(na+1)] = { 0 };
    double sc[(nc+1)][(nc+1)][(na+1)] = { 0 }; 
    double sa[(na+1)][(na+1)][(nc+1)] = { 0 };
    double tc[(nc+1)][(nc+1)] = { 0 }; double ta[(na+1)][(na+1)] = { 0 };
    double lc[(nn+1)][(nc+1)] = { 0 }; double la[(nn+1)][(na+1)] = { 0 };
    double xi[(nn+1)][(nc+1)][(na+1)] = { 0 };
    
    for (int i=1; i<=nc; i++) {
        getline(file, line);
        stringstream linestream(line);
        string value;

        getline(linestream, value, ',');
        string x_S = value;
        getline(linestream, value, ',');
        nzc[i] = stoi(value);
    }
    for (int i=1; i<=na; i++) {
        getline(file, line);
        stringstream linestream(line);
        string value;

        getline(linestream, value, ',');
        string x_S = value;
        getline(linestream, value, ',');
        nza[i] = stoi(value);
    }
    {
        getline(file, line);
        stringstream linestream(line);
        string value;

        getline(linestream, value, ',');
        string x_S = value;
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

        ap0 = temperature(at, bt, ct, dt, et, temp);
    }

    for (int i=1; i<=nc; i++) {
        for (int j=1; j<=na; j++) {
            {
                getline(file, line);
                stringstream linestream(line);
                string value;

                getline(linestream, value, ',');
                string x_S = value;
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

                b0[i][j] = temperature(at, bt, ct, dt, et, temp);
            }
            {
                getline(file, line);
                stringstream linestream(line);
                string value;

                getline(linestream, value, ',');
                string x_S = value;
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

                b1[i][j] = temperature(at, bt, ct, dt, et, temp);
            }
            {
                getline(file, line);
                stringstream linestream(line);
                string value;

                getline(linestream, value, ',');
                string x_S = value;
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

                b2[i][j] = temperature(at, bt, ct, dt, et, temp);
            }
            {
                getline(file, line);
                stringstream linestream(line);
                string value;

                getline(linestream, value, ',');
                string x_S = value;
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

                c0[i][j] = temperature(at, bt, ct, dt, et, temp);
            }
        }
    }
    for (int i=1; i<=nc-1; i++) {
        for (int j=i+1; j<=nc; j++) {
            getline(file, line);
            stringstream linestream(line);
            string value;

            getline(linestream, value, ',');
            string x_S = value;
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

            tc[i][j] = temperature(at, bt, ct, dt, et, temp);
            tc[j][i] = tc[i][j];
        }
    }
    for (int i=1; i<=na-1; i++) {
        for (int j=i+1; j<=na; j++) {
            getline(file, line);
            stringstream linestream(line);
            string value;

            getline(linestream, value, ',');
            string x_S = value;
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

            ta[i][j] = temperature(at, bt, ct, dt, et, temp);
            ta[j][i] = tc[i][j];
        }
    }
    for (int k=1; k<=nc-1; k++) {
        for (int i=k+1; i<=nc; i++) {
            for (int j=1; j<=na; j++) {
            getline(file, line);
            stringstream linestream(line);
            string value;

            getline(linestream, value, ',');
            string x_S = value;
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

            sc[k][i][j] = temperature(at, bt, ct, dt, et, temp);
            sc[i][k][j] = sc[k][i][j];
            }
        }
    }
    for (int k=1; k<=na-1; k++) {
        for (int i=k+1; i<=na; i++) {
            for (int j=1; j<=nc; j++) {
            getline(file, line);
            stringstream linestream(line);
            string value;

            getline(linestream, value, ',');
            string x_S = value;
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

            sa[k][i][j] = temperature(at, bt, ct, dt, et, temp);
            sa[i][k][j] = sa[k][i][j];
            }
        }
    }
    for (int i=1; i<=nn; i++) {
        for (int j=1; j<=nc; j++) {
            getline(file, line);
            stringstream linestream(line);
            string value;

            getline(linestream, value, ',');
            string x_S = value;
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

            lc[i][j] = temperature(at, bt, ct, dt, et, temp);
        }
    }
    for (int i=1; i<=nn; i++) {
        for (int j=1; j<=na; j++) {
            getline(file, line);
            stringstream linestream(line);
            string value;

            getline(linestream, value, ',');
            string x_S = value;
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

            la[i][j] = temperature(at, bt, ct, dt, et, temp);
        }
    }
    for (int k=1; k<=nn; k++) {
        for (int i=1; i<=nc; i++) {
            for (int j=1; j<=na; j++) {
            xi[k][i][j] = 0;
            }
        }
    }
    xi[2][9][1] = -0.0102e0;
    xi[2][1][2] = 0.046e0;
    
    double ec[(nc+1)][(nc+1)] = { 0 }; double ea[(na+1)][(na+1)] = { 0 };
    double fc[(nc+1)][(nc+1)] = { 0 }; double fa[(na+1)][(na+1)] = { 0 };
    double xc[(nc+1)][(nc+1)] = { 0 }; double xa[(na+1)][(na+1)] = { 0 };
    double pp[(nc+1)][(nc+1)] = { 0 }; double qp[(na+1)][(na+1)] = { 0 };
    double p[(nc+1)][(nc+1)] = { 0 }; double q[(na+1)][(na+1)] = { 0 };
    double pf[(nc+1)][(nc+1)] = { 0 }; double qf[(na+1)][(na+1)] = { 0 };
    double cc[(nc+1)][(na+1)] = { 0 }; double bf[(nc+1)][(na+1)] = { 0 };
    double b[(nc+1)][(na+1)] = { 0 }; double bp[(nc+1)][(na+1)] = { 0 };
    double gc[(nc+1)] = { 0 }; double ga[(na+1)]; double gn[(nn+1)] = { 0 };

    bp0 = 1.2e0; mh2o = 55.51e0;

    u = 0; z = 0;
    for (int i=1; i<=nc; i++) {
        u += c[i] * pow(nzc[i], 2);
        z += c[i] * nzc[i];
    }
    for (int j=1; j<=na; j++) {
        u += a[j] * pow(nza[j], 2);
        z += a[j] * nza[j];
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
                xc[j][i] = pow(nzc[j], 2) * u;
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
                xa[i][j] = 2 * u; xa[i][i] = pow(nza[i], 2) * u;
                xa[j][j] = pow(nza[j], 2) * u;
                ea[i][j] = (j0(xa[i][j]) - j0(xa[i][i]) / 2 - j0(xa[j][j]) / 2) / fi / 2;
                fa[i][j] = ((xa[i][j] * j1(xa[i][j]) - xa[i][i] * j1(xa[i][i]) / 
                                2 - xa[j][j] * j1(xa[j][j]) / 2) / pow(fi, 2) / 4 - ea[i][j] / fi);
                ea[j][i] = ea[i][j]; fa[j][i] = fa[i][j];
            }
        }
    }
    for (int i=1; i<=nc-1; i++) {
        for (int j=i+1; j<=nc; j++) {
            pp[i][j] = fc[i][j]; p[i][j] = tc[i][j] + ec[i][j];
            pf[i][j] = p[i][j] + pp[i][j] * fi;
            pp[j][i] = pp[i][j]; p[j][i] = p[i][j]; pf[j][i] = pf[i][j];
        }
    }
    for (int i=1; i<=na-1; i++) {
        for (int j=i+1; j<=na; j++) {
            qp[i][j] = fa[i][j]; q[i][j] = ta[i][j] + ea[i][j];
            qf[i][j] = q[i][j] + qp[i][j] * fi;
            qp[j][i] = qp[i][j]; q[j][i] = q[i][j]; qf[j][i] = qf[i][j];
        }
    }
    w = fj * 12;
    for (int i=1; i<=nc; i++) {
        for (int j=1; j<=na; j++) {
            cc[i][j] = c0[i][j] / sqrt(nzc[i] * nza[j]) / 2;
            if (nzc[i] == 2 and nza[j] == 2) v = fj * 1.4e0;
            if (nzc[i] == 1 or nza[j] == 1) v = fj * 2;
            bf[i][j] = b0[i][j] + b1[i][j] * exp(-v) + b2[i][j] * exp(-w);
            b[i][j] = b0[i][j] + b1[i][j] * (g0(v)) + b2[i][j] * (g0(w));
            bp[i][j] = b1[i][j] * (g1(v)) / fi + b2[i][j] * (g1(w)) / fi;
        }
    }
    f = -ap0 * (fj / (1 + bp0 * fj) + 2 / bp0 * log(1 + bp0 * fj));
    for (int i=1; i<=nc; i++) {
        for (int j=1; j<=na; j++) {
            f += c[i] * a[j] * bp[i][j];
        }
    }
    for (int i=1; i<=(nc-1); i++) {
        for (int j=(i+1); j<= nc; j++) {
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
            u += a[j] * b[ii][j] * 2 + z * cc[ii][j];
        } 
        for (int i=1; i<=nc; i++) {
            if (i != 11) {
                v = 0;
                for (int j=1; j<=na; j++) {
                    v += a[j] * sc[ii][i][j];
                }
                u += c[i] * p[ii][i] * 2 + v;
            }
        }
        for (int i=1; i<=(na-1); i++) {
            for (int j=(i+1); j<=na; j++) {
                u += a[i] * a[j] * sa[i][j][ii];
            }
        }
        for (int i=1; i<=nc; i++) {
            for (int j=1; j<=na; j++) {
                u += c[i] * a[j] * cc[i][j] * nzc[ii];
            }
        }
        for (int i=1; i<=nn; i++) {
            u += h[i] * lc[i][ii] * 2;
        }
        for (int k=1; k<=nn; k++) {
            for (int j=1; j<=na; j++) {
                u += h[k] * a[j] * xi[k][ii][j];
            }
        }
        gc[ii] = exp(u);
    }
    // line 1610

}

int main() {
    // data common
    const int max_constit = 30; const int max_conv = 100;
    double pk; double pkstd;
    const double epsilon = 1.0e-8; 

    // declare and initialize constants
    const int n = 25; const int ntot = 12; const int n0 = 10; 
    const int max_cat = 5; const int max_an = 5;

    // declare dynamic arrays
    double* tot; double* tot0; double* totinit; double* molal;
    double* act; double* gact; double* atom; double* zz;
    double* xx; double* psc; double* cat; double* ani;
    int* nchcat; int* nchani; int* nch; int* ica;
    string* constit_S; string* aq_S; string* mineral_S; string* mineral0_S;
    string* nom_ion_S;
    double* c_ion;

    // declare variables
    int i, j, k, l, ii, nminer, nm, nt, nc, na, ncomp, icat, iani, nconstit;
    int nm0, nu, ni, nj, nconv, num, nwmp, nwm, kinvariant, nw, nmneuf, npasfich, npasimp, ncompt, npasecran;
    string num_S, unit_S, x_S, e_S, repinit_S, min_S, pco2_S, pco_S, syst_S, rep_S;
    string zone_S, m_S, p_S, y_S, fichmin_S, molemin_S, fich_S, prnt_S, constituant_S;
    string at_S, bt_S, ct_S, dt_S, et_S;
    double temp, dens, ph, at, bt, ct, dt, et, pk0, pkf, pkmol, pkeq, dph, dil, diltot, mh2o;
    double ap0, bp0;
    double a, b, c, u, sc, cmax, sa, amax, dca, delta, xu, eq, s, ee, stdi, aw, fi, std, ef, tinit, stdmax;
    double po, po0, poinit, phinit, ph0, pc, poa, alcar, albor, alsil, aloh, alh, altest, d, v, xi, eps;

    string label, log_file, event_file, chem_file, min_file;
    double pco2, max_salinity;
    string system; string output_units, units;
    double dilute;
    string add_minerals, rem_minerals;
    vector<string> add_min; vector<string> rem_min;
    int print_step, output_step, verbose, output;
    ifstream file; 
    string line;

    // Allocate dynamic arrays
    constit_S = new (nothrow) string [max_constit+1];
    psc = new (nothrow) double [15] {0};
    tot = new (nothrow) double [ntot+1] {0};
    tot0 = new (nothrow) double [ntot+1] {0};
    totinit = new (nothrow) double [ntot+1] {0};
    nch = new (nothrow) int [n+1] {0};
    molal = new (nothrow) double [n+1] {0};
    act = new (nothrow) double [n+1] {0};
    gact = new (nothrow) double [n+1] {0};
    aq_S = new (nothrow) string [n+1];
    atom = new (nothrow) double [n+1] {0};
    auto kmat = new (nothrow) int [n+1][n+1]();
    ica = new (nothrow) int [n+1] {0};
    auto z = new (nothrow) double [n+1][n+1]();
    zz = new (nothrow) double [n+1] {0};
    xx = new (nothrow) double [n+1] {0};
    cat = new (nothrow) double [max_cat+1] {0};
    ani = new (nothrow) double [max_an+1] {0};
    nchcat = new (nothrow) int [max_cat+1] {0};
    nchani = new (nothrow) int [max_an+1] {0};

    // read input
    {
        // Declare temporary input chemistry variables
        double na, k, li, ca, mg, cl, so4, no3, b, si, alk;

        // Open the input file into the file stream
        file.open("example.txt");

            // Read the input data into simulation attributes
            file >> label; file >> temp; file >> dens; 
            file >> ph; file >> na; file >> k;
            file >> li; file >> ca; file >> mg;
            file >> cl; file >> so4; file >> no3;
            file >> b; file >> si; file >> alk;
            file >> pco2; file >> system; 
            file >> output_units; file >> dilute; 
            
            // Parse minerals to be added to database
            {
                file >> x_S; 
                stringstream linestream(x_S);
                while (linestream.good()) {
                    string substr;
                    getline(linestream, substr, ',');
                    add_min.push_back(lower_to_upper(substr));
                }     
            }
            
            // Parse minerals to be removed from database
            {
                file >> x_S; 
                stringstream linestream(x_S);
                while (linestream.good()) {
                    string substr;
                    getline(linestream, substr, ',');
                    rem_min.push_back(lower_to_upper(substr));
                }     
            }
            
            file >> max_salinity;
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

    // Set outfile filenames
    log_file = label + ".log";
    event_file = label + ".j" + system + "@";
    chem_file = label + ".j" + system + "&";
    min_file = label + ".j" + system + "%";

    // EQL greeting (line 114 of f90 source)
    if (verbose == 1) {
        time_t now;
        time (&now);

        cout << "\nThis is EQL..............\n";
        cout << asctime(localtime(&now)) << "\n";
    };

    // Set more variables
    pk = 0.1e0; pk0 = pk; pkf = 0.0001e0;
    pkstd = 0.1e0; dph = 0.2e0; dil = 0.; 
    diltot = 1.; mh2o = 55.51e0; 
    nconv = 2; eps = 1.0e-12;
    

    // Set initial temperature to temperature
    tinit = temp;

    // Set nminer based on nonzero input chemistry
    nminer = 0;
    for (int i = 1; i <= ntot; i++) {
      if (tot[i] > 0) nminer ++;
    }

    // Interpret units based on reported density
    if (dens == 0) dens = 1;
    if (dens == 1) {
      units = "molal";
    } else {
      units = "molar";
    }

    // Open and read "aqu.dat"
    file.open("aqu.dat");
    for (int i = 0; i <= n; i++) {
      getline(file, line);
      stringstream linestream(line);
      string value;

      getline(linestream, value, ',');
      aq_S[i] = value;

      getline(linestream, value, ',');
      atom[i] = stod(value);

      getline(linestream, value, ',');
      nch[i] = stoi(value);
    }
    file.close();

    // Open and read "matrice1"
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

    {
        // Open and read murtf2
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
    
                if (x_S == aq_S[j]) {
    
                    mu[j] = (at + 
                             bt / 300 * temp + 
                             ct / 30000 * pow(temp, 2) +
                             dt / 3000000 * pow(temp, 3) + 
                             et / 300000000 * pow(temp, 4));
                }
            }
        }

        for (int k=1; k<=nm; k++) {
            getline(file, line);
            stringstream linestream(line);
            string value, x_S, nom_ion_S;
            double at, bt, ct, dt, et;
            int ncomp;
            double c_ion;
    
            getline(linestream, value, ',');
            mineral_S[k] = value;
            getline(linestream, value, ',');
            ncomp = stoi(value);
    
            for (int i=1; i<=ncomp; i++) {
                getline(linestream, value, ',');
                c_ion = stod(value);
                getline(linestream, value, ',');
                nom_ion_S = value;
                x_S = upper_to_lower(nom_ion_S);

                for (int j=0; j<=nt; j++) {
                    if (x_S == aq_S[j]) {
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

    // Set psol
    for (int k=1; k<=nm; k++) {
      double u = mum[k];
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

        for (int k=1; k<=nm; k++) {
            getline(file, line);
            stringstream linestream(line);
            string value;
            getline(linestream, value, ',');
            mineral0_S[k] = value;
        }
        for (int k=1; k<=nm; k++) {
            nwmin[k] = 0;
        }
        for (int k=1; k<=nm; k++) {
            for (int l=1; l<=nm0; l++) {
                if (mineral0_S[l] == mineral_S[k]) {
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
LOOP500:
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

    s = 0;
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

LOOP200: 
    do {
        int nu = 1;
        int ncompt = 0;
        while (nu != 0) {
        actp(molal, gact, aw, fi, temp, ap0, bp0, mh2o);
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
            for (int j=1; j<=n; j++) {
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

        int ni = n; int nj = n;
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
                double u = z[k-1][k-1] / z[i][k-1];
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
                goto STOP;
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
        exit;
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
    } while (true);

    if (units == "molal" and dil == 0) {
        //line 525
    }


STOP:
    // delete dynamically allocated arrays
    delete[] tot, tot0, totinit, molal, act, gact, atom;
    delete[] zz, xx, psc, cat, ani;
    delete[] z;
    delete[] kmat;
    delete[] nchcat, nchani, nch, ica;
    delete[] constit_S, aq_S, mineral_S, mineral0_S;
    delete[] nom_ion_S;
    delete[] c_ion;

    // Successful exit
    return 0;
}
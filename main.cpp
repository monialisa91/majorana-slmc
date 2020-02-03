#include <iostream>
#include <armadillo>
#include "math.h"
#include "hamiltonian.h"

using namespace std;
using namespace arma;

int main() {
// 1. CONSTANTS

    int n, lx, ly;
    double t, delta, alphay, alphaz, Bx, Bz, cp, temp, J_spin, gamma_N, phi, disorder;
    bool magnetic_disorder;
    temp = 1.0;
    lx = 3; // dlugosc lancucha
    ly = 2; // drugi wymiar
    delta = 0.27;
    alphay = 0;	// Rashba SOC w kierunku y
    alphaz = 0;  // Rashba SOC w kierunku z 0.15
    Bx = 0; // pole magnetyczne w kierunku x
    Bz = 0;
    J_spin=-1;
    gamma_N=0.01;
    temp=0.00001;
    cp=1.25;
    t = -1;
    disorder=0;
    magnetic_disorder = false;
    n = lx*ly;

    vec ang_xy (n);

    for(int j=0; j<n; j++) {
        ang_xy(j) = 2*M_PI*((double) rand() / (RAND_MAX));
    }

    cx_mat test;
    test = hamiltonian(lx, ly, cp, t, delta, J_spin, ang_xy);

    cout << IsHermitian(test, n) << endl;



}
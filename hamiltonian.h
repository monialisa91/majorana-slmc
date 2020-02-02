//
// Created by fraktal on 27.01.2020.
//

#ifndef RASHBA_HAMILTONIAN_H
#define RASHBA_HAMILTONIAN_H

#include <iostream>
#include <armadillo>
#include "math.h"

using namespace std;
using namespace arma;

cx_mat hamiltonian(int lx, int ly, double mu, double t, double delta, double J_spin, vec ang_xy){
    // tau matrices
    int n = lx*ly;
    int i, j;
    cx_mat ham = zeros<cx_mat>(4*n,4*n);
    for(int ix=0; ix<lx; ix++) {
        for(int iy=0; iy<ly; iy++) {
            i = ix + lx*iy;

            for(int idxy=0; idxy<2; idxy++){
                if(idxy == 0){
                    if(ix == lx-1) continue;
                    j = i+1;
                }
                else {
                    if (iy == ly-1) continue;
                    j = i + lx;
                }


                ham(i,j) = t;            //  C^+_up C_up
                ham(j,i) = t;
                ham (i, i) = -mu;

                ham(i+n,j+n) = t;         //  C^+_dn C_dn
                ham(j+n,i+n) = t;
                ham(i+n, i+n) = -mu;


                ham(i+2*n,j+2*n) = -t;   // -C^+_dn C_dn
                ham(j+2*n,i+2*n) = -t;
                ham(i+2*n, i+2*n) = mu;


                ham(i+3*n,j+3*n) = -t;   // -C^+_up C_up
                ham(j+3*n,i+3*n) = -t;
                ham(i+3*n, i+3*n) = mu;

                ham(i, i+2*n) = delta;
                ham(i+2*n , i) = delta;
                ham (i+n,i+3*n) = delta;
                ham(i+3*n, i+n) = delta;

                cx_double angle = cx_double(cos(ang_xy(i)), -sin(ang_xy(i)));
                cx_double angle_conj = cx_double(cos(ang_xy(i)), sin(ang_xy(i)));


                ham(i,i+n) = 0.5*J_spin*angle;
                ham(i+n, i) = 0.5*J_spin*angle_conj;
                ham(i+2*n,i+3*n) = 0.5*J_spin*angle;
                ham(i+3*n, i+2*n) = 0.5*J_spin*angle_conj;

            }


        }

    }

    ham(n-1, n-1) = -mu;
    ham(2*n-1, 2*n-1) = -mu;
    ham(3*n-1, 3*n-1) = mu;
    ham(4*n-1, 4*n-1) = mu;
    ham(n-1, 3*n-1) = delta;
    ham(3*n-1, n-1) = delta;
    ham(2*n-1,4*n-1) = delta;
    ham(4*n-1, 2*n-1) = delta;

    return ham;
}


double energy(cx_mat ham, double beta, double cp) {
    int n = ham.n_rows;
    vec eigval;
    eig_sym(eigval, ham);
    double T = 1.0/beta;
    double E = 0;
    for(int i=0; i<4*n; i++) {
        E += log(1+exp(-beta*(eigval(i)-cp)));
    }
    return -T*E;

}




#endif //RASHBA_HAMILTONIAN_H

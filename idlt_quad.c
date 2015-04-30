/*******************************************************************************
 * Quad precision implementation of discrete Legendre transform (DLT). 
 * Copyright (c) 2015 Nick Hale and Alex Townsend
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <quadmath.h>

void main( void);
void drive(int N, long double *v, long double *x, long double *w, long double *c);
void gaussq(__float128 **a, __float128 *b, __float128 *x, int n);

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void drive(int N, long double *v, long double *x, long double *w, long double *c)
{
    long int i, j, l;
    __float128 xi, vi, P, Pm1, Pm2, n, Pp, Ppm1, Ppm2;
    __float128 **a, *cq, *vq;
    FILE * fid;
    
    // Dynamic storage:
    a = calloc(N, sizeof(__float128 *)); --a;
    for ( i = 1; i <= N; ++i ) {
        a[i] =  calloc(N, sizeof(__float128)); --a[i];
    }
    vq = calloc(N, sizeof(__float128)); --vq;
    cq = calloc(N, sizeof(__float128)); --cq;

    for ( i = 0 ; i < N ; i++ ) {
        xi = ((__float128)x[i]);
        
        for ( l = 0 ; l < 2 ; l++ ) {
            // Improve accuracy of root via Newton.
            Pm2 = 1.;
            Pm1 = xi;
            Ppm2 = 0.;
            Ppm1 = 1.;
            for ( j = 2 ; j <= N ; j++ ) {
                n = ((__float128)j);
                P = (2. - 1./n) * Pm1 * xi - (1. - 1./n) * Pm2;
                Pp = (2. - 1./n) * (Pm1 + Ppm1 * xi) - (1. - 1./n) * Ppm2;
                
                Pm2 = Pm1; Pm1 = P;
                Ppm2 = Ppm1; Ppm1 = Pp;
            }
            xi -= P / Pp;
        }
        
        Pm2 = 1.;
        Pm1 = xi;
        a[i+1][1] = Pm2;
        a[i+1][2] = Pm1;
        for ( j = 2 ; j < N ; j++ ) {
            n = ((__float128)j);
            P = (2. - 1./n) * Pm1 * xi - (1. - 1./n) * Pm2;
            a[i+1][j+1] = P;
            Pm2 = Pm1;
            Pm1 = P;
        }
    }
    
    for ( i = 0 ; i < N ; i++ ) {
         vq[i+1] = (__float128)v[i];
    }
    
    gaussq(a, vq, cq, N);
    
    for ( i = 0 ; i < N ; i++ ) {
        c[i] = (double)cq[i+1];
    }
    
    for ( i = 1; i <= N; ++i ) {
        free(++a[i]);
    }
    free(++a);
    free(++vq);
    free(++cq);
    
    return;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* Gaussian elimination code. */
void gaussq(__float128 **a, __float128 *b, __float128 *x, int n)
{
    int i,j,k,m;
    __float128 xfac, temp, temp1, amax;
    
    /* Forward reduction */
    for ( k = 1; k <= n - 1; ++k ) {
        amax = (__float128) fabsq(a[k][k]) ;
        m = k;
        
        /* Determine pivot */
        for ( i = k + 1; i <= n; i++ ){   
            xfac = (__float128) fabsq(a[i][k]);
            if ( xfac > amax ) {
                amax = xfac; 
                m = i;
            }
        }
        
        /* Row interchanges */
        if ( m != k ) {
            temp1 = b[k];
            b[k]  = b[m];
            b[m]  = temp1;
            for ( j = k; j <= n; j++ ) {
                temp = a[k][j];
                a[k][j] = a[m][j];
                a[m][j] = temp;
            }
        }
        
        /* Row ellimination */
        for ( i = k + 1; i <= n; ++i ) {
            xfac = a[i][k]/a[k][k];
            for ( j = k + 1; j <= n; ++j ) {
                a[i][j] = a[i][j] - xfac*a[k][j];
            }
            b[i] = b[i] - xfac*b[k];
        }   
    }
    
    /* Back substitution */
    for ( j = 1; j <= n; ++j ) {
        k = n - j + 1;
        x[k] = b[k];
        for ( i = k + 1; i <= n; ++i ) {
            x[k] = x[k] - a[k][i]*x[i];
        }
        x[k] = x[k]/a[k][k];
    }
    
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void main( void ) {
    
    int k, N = 0;
    long double v[100002];
    long double x[100002];
    long double w[100002];
    long double c[100002];
    FILE * fid;
    ;
    fid = fopen("data.bin", "r");
    while (!feof(fid))
    {
        if ( fscanf(fid, "%Lf", &c[N]) == 1 )
            ++N;
    }
    fclose(fid);
    
    fid = fopen("data_x.bin", "r");
    for ( k = 0; k < N ; k++ ) {
        fscanf(fid, "%Lf", &x[k]);
    }
    fclose(fid);
    
    fid = fopen("data_w.bin", "r");
    for ( k = 0; k < N ; k++ ) {
        fscanf(fid, "%Lf", &w[k]);
    }
    fclose(fid);
    
    drive( N, c, x, w, v );
    
    fid = fopen("results.bin", "w+");
    for ( k = 0; k < N ; k++ )
        fprintf(fid, "%33.33Lg\n", v[k]);
    fclose(fid);
}
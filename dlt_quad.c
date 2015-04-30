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
#include <quadmath.h> 

void main( void);
void drive(int N, long double *c, long double *t, long double *v);

void drive(int N, long double *c, long double *t, long double *v)
{   
    long int i, j, l;
    __float128 ti, x, vi, P, Pm1, Pm2, n, one, two, Pp, Ppm1, Ppm2;
    
    one = ((__float128)1);
    two = ((__float128)2);
    
    for ( i = 0 ; i < N ; i++ ) {
        ti = ((__float128)t[i]);
        x = cosq(ti);
        
        for ( l = 0 ; l < 2 ; l++ ) {
            /* Improve accuracy of root via Newton. */
            Pm2 = one;
            Pm1 = x;
            Ppm2 = 0;
            Ppm1 = one;
            for ( j = 2 ; j <= N ; j++ ) {
                n = ((__float128)j);
                P = (2. - 1./n) * Pm1 * x - (1. - 1./n) * Pm2; 
                Pp = (2. - 1./n) * (Pm1 + Ppm1 * x) - (1. - 1./n) * Ppm2; 
                
                Pm2 = Pm1; Pm1 = P;
                Ppm2 = Ppm1; Ppm1 = Pp;
            }
            x -= P / Pp;
        }
        
        Pm2 = one;
        Pm1 = x;
        vi = ((__float128)c[0]) + x*((__float128)c[1]);
        for ( j = 2 ; j < N ; j++ ) {
            n = ((__float128)j);
            P = (2. - 1./n) * Pm1 * x - (1. - 1./n) * Pm2; 
            vi += P * ((__float128)c[j]);
            Pm2 = Pm1;
            Pm1 = P;
        }
        v[i] = vi;

    }
   
    return;
}

void main( void ) {
    
    int k, N = 0;
    long double v[100002];
    long double t[100002];
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
    
    fid = fopen("data_t.bin", "r");
    for ( k = 0; k < N ; k++ ) {
        fscanf(fid, "%Lf", &t[k]);
    }  
    fclose(fid);
    
    drive( N, c, t, v );
    
    fid = fopen("results.bin", "w+");
    for ( k = 0; k < N ; k++ )
        fprintf(fid, "%33.33Lg\n", v[k]);
    fclose(fid);
}
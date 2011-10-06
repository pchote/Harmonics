/*
* Copyright 2011 Paul Chote
* This file is part of harmonics, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include <gsl_sf_legendre.h>

void set_color_table()
{
    float l[9] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float r[9] = {0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
    float g[9] = {0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
    float b[9] = {0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
    cpgctab(l, r, g, b, 9, -1.0, 0.5);
    
    // Set the lowest index to black
    int icmin,icmax;
    cpgqcir(&icmin,&icmax);
    cpgscr(icmin, 0, 0, 0);
}

int plot_harmonic(int l, int m)
{
    if (cpgopen("6/xs") <= 0)
    {
        fprintf(stderr, "Unable to open PGPLOT window");
        return 1;
    }

    // 960 x 540
    //cpgpap(11.29, 0.562);
    cpgask(0);
    cpgslw(3);
    cpgscf(2);
    cpgsfs(2);
    set_color_table();

    int ang_res = 300;
    int proj_res = 300;

    /*
     * Plot the harmonic value as a function of theta,phi
     */
    
    // Initial sample
    cpgsvp(0.1, 0.48, 0.1, 0.9);
    cpgwnad(0, ang_res, 0, ang_res);

    float min = 1e9;
    float max = -1e9;
    float *angle_data = (float *)malloc(ang_res*ang_res*sizeof(float));
    for (int t = 0; t < ang_res; t++)
        for (int p = 0; p < ang_res; p++)
        {
            double theta = t*M_PI/ang_res;
            double phi = p*2*M_PI/ang_res;
            
            int tp = ang_res*t + p;
            float val = (float)cos(m*phi)*gsl_sf_legendre_sphPlm(l, m, cos(theta));
            angle_data[tp] = val;
            min = fmin(min, val);
            max = fmax(max, val);
        }
    
    // Compress color range slightly so we can use the bottom-most index for a null color
    min *=1.1;
    max *=1.1;
    
    cpgsci(1);
    cpgswin(0, 2*M_PI, 0, M_PI);
    cpgbox("bcstn", M_PI, 4, "bcstn", M_PI, 4);
    cpgmtxt("l", 2.5, 0.5, 0.5, "theta");
    cpgmtxt("b", 2.5, 0.5, 0.5, "phi");
    
    
    /*
     * Plot the harmonic value on the surface of a sphere
     */

    float *proj_data = (float *)malloc(proj_res*proj_res*sizeof(float));
    float tr[] = {-0.5, 1, 0, -0.5, 0, 1};
    
    double rt = 0;
    double pt = -M_PI/4;
    double t = 0;
    while (1)
    {
        // Update theta/phi plot
        cpgsvp(0.1, 0.48, 0.1, 0.9);
        cpgwnad(0, ang_res, 0, ang_res);
        
        for (int j = 0; j < ang_res; j++)
            for (int i = 0; i < ang_res; i++)
            {
                double theta = j*M_PI/ang_res;
                double phi = i*2*M_PI/ang_res;
                
                int tp = ang_res*j + i;
                float val = (float)cos(m*phi)*gsl_sf_legendre_sphPlm(l, m, cos(theta))*cos(t);
                angle_data[tp] = val;
                min = fmin(min, val);
                max = fmax(max, val);
            }
        
        cpgimag(angle_data, ang_res, ang_res, 1, ang_res, 1, ang_res, min, max, tr);

        // Draw 3d projection
        cpgsvp(0.52, 0.9, 0.1, 0.9);
        cpgwnad(0, proj_res, 0, proj_res);

        // Project from x,y screen coords to theta,phi surface coords
        for (int j = 0; j < proj_res; j++)
            for (int i = 0; i < proj_res; i++)
            {
                int ij = proj_res*j + i;
                proj_data[ij] = -1e9;
                
                double x = i*4.0f/proj_res - 2;
                double y = j*4.0f/proj_res - 2;
                double z = sqrt(1 - x*x - y*y);
                
                // Rotate away from the plane so we are looking slightly down on the rotation axis
                double xxx = x;
                double yyy = y*cos(pt) - z*sin(pt);
                double zzz = y*sin(pt) + z*cos(pt);
                
                
                // Rotate around the axis with time
                double xx = xxx*cos(rt) + zzz*sin(rt);
                double yy = yyy;
                
                double theta = acos(yy);
                double phi = acos(xx/sin(theta));
                
                if (x*x + y*y >= 1)
                    continue;
                
                float val = (float)cos(m*phi)*gsl_sf_legendre_sphPlm(l, m, cos(theta))*cos(t);
                proj_data[ij] = val;
            }

        cpgimag(proj_data, proj_res, proj_res, 1, proj_res, 1, proj_res, min, max, tr);
        
        // Define the coordinate axis lines
        double b = 75;
        double e = 140;
        double h = 145;
        double points[36] = 
        {
            // Heads
            b, 0, 0,
            h, 0, 0,
            0, 0, b,
            0, 0, h,
            0, b, 0,
            0, h, 0,
            
            // Tails
            -b, 0, 0,
            -e, 0, 0,
            0, 0, -b,
            0, 0, -e,
            0, -b, 0,
            0, -e, 0,
        };
        double tp[24];
        
        for (int i = 0; i < 18; i++)
        {
            double x = points[3*i + 0];
            double y = points[3*i + 1];
            double z = points[3*i + 2];
                        
            // Rotate around the axis with time
            double xx = x*cos(rt) + z*sin(rt);
            double yy = y;
            double zz = -x*sin(rt) + z*cos(rt);
            
            // Rotate away from the plane so we are looking slightly down on the rotation axis
            double xxx = xx;
            double yyy = yy*cos(pt) - zz*sin(pt);
            double zzz = yy*sin(pt) + zz*cos(pt);
            
            tp[2*i+0] = xxx;
            tp[2*i+1] = yyy;
        }
        
        float offset[] = {150,150};
        int axis_colors[3] = {2,3,5};
        
        cpgsci(1);
        for (int i = 0; i < 6; i++)
        {
            // Correct for occlusion of the axis by the sphere
            if ((i != 2 && tp[4*i+1] > 0) || i == 5)
            {
                // Line from p1 to p2
                double dpx = tp[4*i+2] - tp[4*i];
                double dpy = tp[4*i+3] - tp[4*i+1];
                
                // Line from c to p1
                double dcx = tp[4*i];
                double dcy = tp[4*i+1];
                
                // Polynomial coefficients
                double a = dpx*dpx + dpy*dpy;
                double b = 2*(dcx*dpx + dcy*dpy);
                double c = dcx*dcx + dcy*dcy - 75*75;
                
                // Solve for line parameter x.
                double d = sqrt(b*b - 4*a*c);
                double x1 = (-b + d)/(2*a);
                double x2 = (-b - d)/(2*a);
                
                // The solution we want will be 0<=x<=1
                double sol = (x1 >= 0 && x1 <= 1) ? x1 : x2;
                if (sol > 0)
                {
                    tp[4*i] += sol*dpx;
                    tp[4*i+1] += sol*dpy;
                }
            }
            
            cpgsci(axis_colors[i%3]);
            if (i < 3)
                cpgarro(offset[0]+tp[4*i],
                        offset[1]+tp[4*i+1],
                        offset[0]+tp[4*i+2],
                        offset[1]+tp[4*i+3]);
            else
            {
                cpgmove(offset[0]+tp[4*i+2],
                        offset[1]+tp[4*i+3]);
                cpgdraw(offset[0]+tp[4*i],
                        offset[1]+tp[4*i+1]);
            }
        }
        rt += 0.02;
        t += 0.15;
    }
    free(angle_data);
    free(proj_data);
    cpgend();
    return 0;
}

int main( int argc, char *argv[] )
{
    if (argc == 3)
        plot_harmonic(atoi(argv[1]), atoi(argv[2]));
    else
    {
        fprintf(stderr, "Invalid args");
        return 1;
    }
    return 0;
}

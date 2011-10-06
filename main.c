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

#pragma mark Generic Helpers

typedef struct
{
    int l;
    int m; // Requires -l <= m <= +l
    float a; // Amplitude
    float w; // Angular frequency
} mode;

/*
 * Calculate the linear combination of spherical harmonics at a given point
 */
float calculate_surface_value(mode *modes, int num_modes, double theta, double phi, double t)
{
    float val = 0;
    for (int i = 0; i < num_modes; i++)
        val += modes[i].a*cos(modes[i].m*phi)*gsl_sf_legendre_sphPlm(modes[i].l, modes[i].m, cos(theta))*cos(modes[i].w*t);
    return val;
}

/*
 * Set the color table for plotting
 * Based on the rainbow color table in the PGPLOT example pgdemo4.f
 */
void set_color_table()
{
    float l[9] = {0.0, 0.005, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float r[9] = {0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
    float g[9] = {0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
    float b[9] = {0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
    cpgctab(l, r, g, b, 9, 1.0, 0.5);
}

/*
 * Determine min and max values for the color palette
 * to best use the available dynamic range
 */
void set_color_scale(mode *modes, int num_modes, int ang_res, float *omin, float *omax)
{
    double min = 1e9;
    double max = -1e9;
    for (int j = 0; j < ang_res; j++)
        for (int i = 0; i < ang_res; i++)
        {
            double val = calculate_surface_value(modes, num_modes, j*M_PI/ang_res, i*2*M_PI/ang_res, 0);
            min = fmin(min, val);
            max = fmax(max, val);
        }

    // Compress color range slightly so we can use the bottom-most index for the background
    min *= 1.2f;
    max *= 1.2f;

    *omin = (float)min;
    *omax = (float)max;
}

/*
 * Rotate a point about the x axis by an angle a
 */
void rotate_x(double x, double y, double z, double a, double *ox, double *oy, double *oz)
{
    *ox = x;
    *oy = y*cos(a) - z*sin(a);
    *oz = y*sin(a) + z*cos(a);
}

/*
 * Rotate a point about the y axis by an angle a
 */
void rotate_y(double x, double y, double z, double a, double *ox, double *oy, double *oz)
{
    *ox = x*cos(a) + z*sin(a);
    *oy = y;
    *oz = -x*sin(a) + z*cos(a);
}

/*
 * Find the 2 values of the line parameter u that describe the
 * intersection of the line between p1 and p2 with the sphere
 * centered at p3 with radius r
 * p1,p2,p3 are double[3]
 * u is a double[2], for the results
 * returns FALSE if there is no intersection
 */
int line_sphere_intersection(double *p1, double *p2, double *p3, double r, double *u)
{
    // Line from p1 to p2
    double dpx = p2[0] - p1[0];
    double dpy = p2[1] - p1[1];
    double dpz = p2[2] - p1[2];

    // Line from p3 to p1
    double dcx = p1[0] - p3[0];
    double dcy = p1[1] - p3[1];
    double dcz = p1[2] - p3[2];

    // Polynomial coefficients
    double a = dpx*dpx + dpy*dpy + dpz*dpz;
    double b = 2*(dcx*dpx + dcy*dpy + dcz*dpz);
    double c = dcx*dcx + dcy*dcy + dcz*dcz - r*r;
    double disc = b*b - 4*a*c;

    // No real solutions
    if (disc < 0)
        return 0;

    // Solve for line parameter u.
    double d = sqrt(disc);
    u[0] = (-b + d)/(2*a);
    u[1] = (-b - d)/(2*a);

    return 1;
}

#pragma mark Drawing Functions

/*
 * Draw the axis lines on top of the projection plot
 */
void plot_projection_axes(double width, double scale, double length, double theta, double phi)
{
    // Sphere radius
    double r = width*scale/2;

    // Axis offset
    float o = width/2;

    // Axis lengths
    double e = r+length;
    double h = e+5;
    double points[36] =
    {
        // Heads
        r, 0, 0,
        h, 0, 0,
        0, 0, r,
        0, 0, h,
        0, r, 0,
        0, h, 0,

        // Tails
        -r, 0, 0,
        -e, 0, 0,
        0, 0, -r,
        0, 0, -e,
        0, -r, 0,
        0, -e, 0,
    };
    double tp[24];
    int axis_colors[3] = {2,3,5};

    for (int i = 0; i < 18; i++)
    {
        double x = points[3*i];
        double y = points[3*i + 1];
        double z = points[3*i + 2];

        // Rotate around the axis with time
        rotate_y(x, y, z, theta, &x, &y, &z);

        // Rotate away from the plane so we are looking slightly down on the rotation axis
        rotate_x(x, y, z, phi, &x, &y, &z);

        tp[2*i+0] = x;
        tp[2*i+1] = y;
    }

    for (int i = 0; i < 6; i++)
    {
        // Correct for occlusion of the axis by the sphere
        if ((i != 2 && tp[4*i+1] > 0) || i == 5)
        {
            double p1[3] = {tp[4*i], tp[4*i+1], 0};
            double p2[3] = {tp[4*i+2], tp[4*i+3], 0};
            double p3[3] = {0, 0, 0};
            double u[2] = {0,0};

            if (line_sphere_intersection(p1, p2, p3, r, u))
            {
                // The solution we want will be positive
                double sol = u[0] >= 0 ? u[0] : u[1];

                // The entire line is hidden if sol > 1
                if (sol > 1)
                    continue;

                tp[4*i] += sol*(p2[0]-p1[0]);
                tp[4*i+1] += sol*(p2[1]-p1[1]);
            }
        }

        cpgsci(axis_colors[i%3]);
        if (i < 3)
        {
            // Head
            cpgarro(o+tp[4*i],
                    o+tp[4*i+1],
                    o+tp[4*i+2],
                    o+tp[4*i+3]);
        }
        else
        {
            // Tail
            cpgmove(o+tp[4*i+2],
                    o+tp[4*i+3]);
            cpgdraw(o+tp[4*i],
                    o+tp[4*i+1]);
        }
    }
}

int plot_harmonic(int l, int m)
{
    /*
     * Simulation setup
     */

    // Pulsation modes
    int num_modes = 1;
    mode modes[] =
    {
        {l, m, 1.0f, 1.0f}
    };

    /*
     * Plot setup
     */

    // Calculation resolution for each subplot
    int ang_res = 300;
    int proj_res = 300;

    // Sphere radius as a fraction of the plot half-width
    double scale = 0.5;

    if (cpgopen("6/xs") <= 0)
    {
        fprintf(stderr, "Unable to open PGPLOT window");
        return 1;
    }

    cpgask(0);
    cpgslw(3);
    cpgscf(2);
    cpgsfs(2);
    cpgsci(1);
    set_color_table();

    float min, max;
    set_color_scale(modes, num_modes, ang_res, &min, &max);

    float *proj_data = (float *)malloc(proj_res*proj_res*sizeof(float));
    float *angle_data = (float *)malloc(ang_res*ang_res*sizeof(float));
    float tr[] = {-0.5, 1, 0, -0.5, 0, 1};

    // Plot border and decorations for angle plot
    cpgsvp(0.1, 0.48, 0.1, 0.9);
    cpgwnad(0, ang_res, 0, ang_res);
    cpgswin(0, 2*M_PI, 0, M_PI);
    cpgbox("bcstn", M_PI, 4, "bcstn", M_PI, 4);
    cpgmtxt("l", 2.5, 0.5, 0.5, "theta");
    cpgmtxt("b", 2.5, 0.5, 0.5, "phi");

    double rt = 0;
    double pt = -M_PI/4;

    // Current simulation time
    double t = 0;
    while (1)
    {
        /*
         * Plot the harmonic value as a function of theta,phi
         */
        cpgsvp(0.1, 0.48, 0.1, 0.9);
        cpgwnad(0, ang_res, 0, ang_res);

        for (int j = 0; j < ang_res; j++)
            for (int i = 0; i < ang_res; i++)
                angle_data[ang_res*j + i] = calculate_surface_value(modes, num_modes, j*M_PI/ang_res, i*2*M_PI/ang_res, t);

        cpgimag(angle_data, ang_res, ang_res, 1, ang_res, 1, ang_res, min, max, tr);

        /*
         * Plot the harmonic value on the surface of a sphere
         */
        cpgsvp(0.52, 0.9, 0.1, 0.9);
        cpgwnad(0, proj_res, 0, proj_res);

        // Project from x,y screen coords to theta,phi surface coords
        for (int j = 0; j < proj_res; j++)
            for (int i = 0; i < proj_res; i++)
            {
                int ij = proj_res*j + i;
                proj_data[ij] = -1e9;

                double x = (2*i - proj_res)/(scale*proj_res);
                double y = (2*j - proj_res)/(scale*proj_res);

                if (x*x + y*y >= 1)
                    continue;

                double z = sqrt(1 - x*x - y*y);

                // Rotate away from the plane so we are looking slightly down on the rotation axis
                rotate_x(x, y, z, pt, &x, &y, &z);

                // Rotate away from the plane so we are looking slightly down on the rotation axis
                rotate_y(x, y, z, rt, &x, &y, &z);

                double theta = acos(y);
                double phi = acos(x/sin(theta));
                proj_data[ij] = calculate_surface_value(modes, num_modes, theta, phi, t);
            }

        cpgimag(proj_data, proj_res, proj_res, 1, proj_res, 1, proj_res, min, max, tr);
        plot_projection_axes(proj_res, scale, 60, rt, pt);
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

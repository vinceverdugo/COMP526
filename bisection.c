#include <stdio.h>  // needed to be able to use printf
#include <stdlib.h> // needed to be able to define some variables in the heap
#include <math.h>   // needed to be able to use some math primitive functions, such as fabs() for absolute value of floating point numbers

// Declare bisection method
int bisect_hist(double *hist, double a, double b, double tol); 

//Declare endpoints as static variables in memory pool and tolerance variable
    static double a1 = -4;
    static double b1 = -2;
    static double a2 = 2;
    static double b2 = 4;

// Define function e^(-x)(3.2sin(x)-0.5cos(x)) = 3
double f(double x0){
    double x;
    x = exp(-x0)*(3.2*sin(x0)-0.5*cos(x0)) - 3;
    return x;
}

// Define bisection method
int bisect_hist(double *hist, double a, double b, double tol) {
    // Create element to access add elements to heap
    int index = 0;

    // Loop that iterates if f(x) crosses the x axis between points x=a and x=b
    while (fabs(b-a) > tol) {
        //Create midpoint of range
        double mid = (a+b)/2;
        // Store midpoint value
        hist[index] = mid;
        // Check to see if [a, mid] has a root
        if (f(a)*f(mid) < 0) {
            // If yes, b=mid and while loop will then iterate under range [a,mid]
            b = mid;
        }
        else {
            // a=mid and while loop will then iterate under range [mid,b]
            a = mid;
        }
        index++;
    }
    return index;
}

int main(void) {
    double tol = pow(10,-4);
    // Allocated two heaps with 100 entries that are of double 
    double *roots1 = malloc(100*sizeof(double));
    double *roots2 = malloc(100*sizeof(double));

    // Total number of tests performed on first interval
    int numTests1 = bisect_hist(roots1, a1, b1, tol);
    // Root approximation (should be)
    double x_star1 = roots1[numTests1 - 1];
    // Check if root approximation is valid 
    if (round(f(x_star1)) == 0) {
        printf("After %i iterations:\n", numTests1);
        printf("x_* = %f, f(x_*) = %f\n", x_star1, f(x_star1));
    }
    // Root approximation is false
    else {
        printf("Error, root could not be found on the interval [%f, %f]", a1, b1);
    }
    // Retest for second interval
    int numTests2 = bisect_hist(roots2, a2, b2, tol);
    double x_star2 = roots2[numTests2 - 1];
    if (round(f(x_star2)) == 0) {
        printf("After %i iterations:\n", numTests2);
        printf("x_* = %f, f(x_*) = %f\n", x_star2, f(x_star2));
    }
    else {
        printf("Error, root could not be found on the interval [%f, %f]", a2, b2);
    }
    // Free memory
    free(roots1);
    free(roots2);
    return 0;
}
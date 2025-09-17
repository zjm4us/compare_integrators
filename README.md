# Compare different quadrature rules for integration

There are two examples provided for calculating the weights and abscissas for gaussian quadrature rules, try:

```
make
./gqconstants
```

or

```
python gqconstants.py
```

You can also use the C++ example as a guide to build your own executable

There is no need to look at rules >~25 for Gaussian quadrature.  And you can also stop at ~ 1000 divisions for the trapezoidal and Simpson's rules.  If you run much longer you'll see the numerical errors bevome visible for the trapezoidal, but hyou'll need to think about how to code efficiently or the running time may be very long.

# Numerical Integration: Error Assessment

This project uses three numerical integration methods to see how accurate they are for different types of functions. The goal was to explore how trapezoid, Simpson, and Gaussian quadrature behave for smooth versus hard-to-integrate functions, following Section 6.2.5 in Landau's *Computational Physics*.

## Files

- `integration_errors.py` – Python script that computes integrals and relative errors.  
- `Errors.png` – Log-log plot of relative errors for a smooth function (\( e^{-t} \)).  
- `BadErrors.png` – Log-log plot of relative errors for a highly oscillatory function (\( \sin(200 t) \)).  

## Methods Used

1. **Trapezoid Rule** – Approximates the integral using linear interpolation. This is the simplest method I tried first to see how well it works.  
2. **Simpson’s Rule** – Uses quadratic interpolation and ensures an even number of intervals. I expected this to be more accurate than trapezoid for smooth functions.  
3. **Gaussian Quadrature** – Uses Gauss-Legendre nodes and weights. This method is supposed to be very precise even with fewer points.

## Error Assessment

- **Smooth Function** (\( f(t) = e^{-t} \)):
  - Trapezoid converges as expected (\( O(N^{-2}) \)).  
  - Simpson converges faster (\( O(N^{-4}) \)).  
  - Gaussian quadrature performed extremely well, reaching near machine precision even with few nodes. This was the most accurate method for the smooth function.

- **Hard-to-Integrate Function** (\( f(t) = \sin(200 t) \)):
  - Trapezoid and Simpson only give reasonable results for very large `N`.  
  - Gaussian quadrature struggles because the standard nodes do not capture the rapid oscillations.  
  - This shows that even a high-precision method can perform poorly if the function varies too quickly between nodes.

## Discussion

- I noticed that the oscillatory function jumps so quickly that the methods can’t catch all the peaks and valleys, which explains why the errors stay large for small `N`.  
- Possible ways to improve accuracy include:  
  - Adaptive quadrature that places more points where the function changes rapidly.  
  - Special methods designed for oscillatory functions, such as Filon-type quadrature or using variable transformations.

## Terminal Output

Below is the output printed by the script:

Errors for Errors for ∫ e^{-t} dt on [0,1]
N= 2 trap=2.075e-02 simp=3.372e-04 gauss=2.243e-04
N= 10 trap=8.332e-04 simp=5.549e-07 gauss=1.000e-16
N= 20 trap=2.083e-04 simp=3.471e-08 gauss=1.756e-16
N= 40 trap=5.208e-05 simp=2.170e-09 gauss=1.756e-16
N= 80 trap=1.302e-05 simp=1.356e-10 gauss=1.756e-16
N= 160 trap=3.255e-06 simp=8.477e-12 gauss=1.756e-16
N= 320 trap=8.138e-07 simp=5.299e-13 gauss=1.756e-16
N= 640 trap=2.035e-07 simp=3.302e-14 gauss=1.756e-16

Errors for Errors for ∫ sin(200x) dx on [0,1]
N= 2 trap=1.849e+02 simp=1.894e+02 gauss=7.506e+01
N= 10 trap=1.442e+01 simp=1.658e+01 gauss=1.064e+01
N= 20 trap=2.479e+00 simp=8.113e+00 gauss=3.619e+01
N= 40 trap=4.347e+00 simp=4.969e+00 gauss=3.619e+01
N= 80 trap=5.847e-01 simp=6.693e-01 gauss=3.619e+01
N= 160 trap=1.337e-01 simp=1.658e-02 gauss=3.619e+01
N= 320 trap=3.277e-02 simp=8.888e-04 gauss=3.619e+01
N= 640 trap=8.151e-03 simp=5.360e-05 gauss=3.619e+01

Plots saved as Errors.png and BadErrors.png


## Usage

To reproduce the results, just run:

```bash
python integration_errors.py


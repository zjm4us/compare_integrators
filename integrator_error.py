import numpy as np
from numpy.polynomial.legendre import leggauss
import matplotlib.pyplot as plt

# Integration methods
def trapezoid(f, a, b, N):
    x = np.linspace(a, b, N+1)
    y = f(x)
    h = (b - a)/N
    return h*(0.5*y[0] + y[1:-1].sum() + 0.5*y[-1])

def simpson(f, a, b, N):
    if N % 2 == 1:
        N += 1
    x = np.linspace(a, b, N+1)
    y = f(x)
    h = (b - a)/N
    return h/3 * (y[0] + y[-1] + 4*y[1:-1:2].sum() + 2*y[2:-1:2].sum())

def gauss_legendre(f, a, b, n):
    nodes, weights = leggauss(n)
    t = 0.5*(nodes + 1)*(b - a) + a
    return 0.5*(b - a) * np.dot(weights, f(t))

# Error analysis
def error_analysis(f, a, b, I_exact, Ns, fname, title):
    trap_errs, simp_errs, gauss_errs = [], [], []
    eps_min = 1e-16
    max_gauss = 20

    for N in Ns:
        I_t = trapezoid(f, a, b, N)
        I_s = simpson(f, a, b, N)
        n_g = min(N, max_gauss)
        I_g = gauss_legendre(f, a, b, n_g)

        trap_errs.append(abs(I_t - I_exact)/abs(I_exact))
        simp_errs.append(abs(I_s - I_exact)/abs(I_exact))
        gauss_errs.append(max(abs(I_g - I_exact)/abs(I_exact), eps_min))

    # print errors
    print(f"\nErrors for {title}")
    for i, N in enumerate(Ns):
        print(f"N={N:4d}  trap={trap_errs[i]:.3e}  simp={simp_errs[i]:.3e}  gauss={gauss_errs[i]:.3e}")

    # plot errors
    plt.figure(figsize=(7,5))
    plt.loglog(Ns, trap_errs, 'o-', label="Trapezoid")
    plt.loglog(Ns, simp_errs, 's-', label="Simpson")
    plt.loglog(Ns, gauss_errs, '^-', label="Gauss-Legendre")
    plt.xlabel("N")
    plt.ylabel("Relative Error")
    plt.title(title)
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig(fname, dpi=200)
    plt.show()

# Main program
if __name__ == "__main__":
    Ns = [2, 10, 20, 40, 80, 160, 320, 640]

    # normal function
    f1 = lambda t: np.exp(-t)
    I_exact1 = 1 - np.exp(-1)
    error_analysis(f1, 0.0, 1.0, I_exact1, Ns, "Errors.png",
                   "Errors for ∫ e^{-t} dt on [0,1]")

    # hard function
    f2 = lambda t: np.sin(200*t)
    I_exact2 = (1 - np.cos(200)) / 200
    error_analysis(f2, 0.0, 1.0, I_exact2, Ns, "BadErrors.png",
                   "Errors for ∫ sin(200x) dx on [0,1]")

    print("\nPlots saved as Errors.png and BadErrors.png")


#include <iostream>
#include <cmath>
#include <iomanip>

double f(double t, double y) {
    double f1 = std::exp(-t * t) * std::cos(0.8 * t);
    double f2 = std::exp(2 * t - 3);
    double g = 2 * (t - 2);
    return f1 * f2 - g * y;
}

double rk4_step(double t, double y, double h) {
    double k1 = f(t, y);
    double k2 = f(t + h / 2.0, y + h * k1 / 2.0);
    double k3 = f(t + h / 2.0, y + h * k2 / 2.0);
    double k4 = f(t + h, y + h * k3);
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

double rk4_adaptive(double& t, double& y, double& h,
    double eps = 1e-5, double epsM = 1e-6) {
    while (true) {
        double y_big = rk4_step(t, y, h);

        double y_half = rk4_step(t, y, h / 2.0);
        y_half = rk4_step(t + h / 2.0, y_half, h / 2.0);

        double error = std::fabs(y_half - y_big) / 15.0;

        if (error < eps) {
            t += h;
            y = y_half;

            double factor = pow(eps / (error + 1e-16), 0.25); 
            if (factor > 2.0) factor = 2.0;  
            if (factor < 0.5) factor = 0.5;  
            h *= factor;

            return error;
        }
        else {
            h *= 0.5;
        }
    }
}


int main() {
    double t0 = 1.0;
    double T = 6.0;
    double y = 10.0;
    double h = 0.1;   
    double eps = 1e-5;
    double epsM = 1e-6;

    std::cout << std::fixed << std::setprecision(10);

    double t = t0;
    while (t < T) {
        double error = rk4_adaptive(t, y, h, eps, epsM);
        std::cout << "t=" << t << " y=" << y << " error=" << error << " h=" << h << "\n";
        if (t + h > T) h = T - t;
    }

    std::cout << "\nFinal result: y(" << T << ") = " << y << "\n";
    return 0;
}

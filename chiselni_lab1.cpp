//#include <iostream>
//#include <cmath>
//#include <iomanip>
//
//double f(double t, double y) {
//    double f1 = std::exp(-t * t) * std::cos(0.8 * t);
//    double f2 = std::exp(2 * t - 3);
//    double g = 2 * (t - 2);
//    return f1 * f2 - g * y;
//}
//
//double rk4_step(double t, double y, double h) {
//    double k1 = f(t, y);
//    double k2 = f(t + h / 2.0, y + h * k1 / 2.0);
//    double k3 = f(t + h / 2.0, y + h * k2 / 2.0);
//    double k4 = f(t + h, y + h * k3);
//    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
//}
//
//double rk4_adaptive(double& t, double& y, double& h,
//    double eps = 1e-5, double epsM = 1e-6) {
//    while (true) {
//        double y_big = rk4_step(t, y, h);
//
//        double y_half = rk4_step(t, y, h / 2.0);
//        y_half = rk4_step(t + h / 2.0, y_half, h / 2.0);
//
//        double error = std::fabs(y_half - y_big) / 15.0;
//
//        if (error < eps) {
//            t += h;
//            y = y_half;
//
//            double factor = pow(eps / (error + 1e-16), 0.25); 
//            if (factor > 2.0) factor = 2.0;  
//            if (factor < 0.5) factor = 0.5;  
//            h *= factor;
//
//            return error;
//        }
//        else {
//            h *= 0.5;
//        }
//    }
//}
//
//
//int main() {
//    double t0 = 1.0;
//    double T = 6.0;
//    double y = 10.0;
//    double h = 0.1;   
//    double eps = 1e-5;
//    double epsM = 1e-6;
//    double solution = 0.000003066334874;
//
//    std::cout << std::fixed << std::setprecision(10);
//
//    int steps = 0;
//    double t = t0;
//
//    std::cout << "t=" << t << " y=" << y << " u(t)=" << rk4_adaptive(t, y, h, eps, epsM) << " |y-u(t)|=" << y - rk4_adaptive(t, y, h, eps, epsM) << "\n";
//    while (t < T) {
//        double error = rk4_adaptive(t, y, h, eps, epsM);
//        std::cout << "t=" << t << " y=" << y << " u(t)=" << rk4_adaptive(t, y, h, eps, epsM) << " |y-u(t)|=" << y - error << "\n";
//        if (t + h > T) h = T - t;
//        steps++;
//    }
//
//    std::cout << "\nFinal result: y(" << T << ") = " << y << "\n";
//    std::cout << "Maple result: " << solution << std::endl;
//    std::cout << "Ammount of steps: " << steps << std::endl;
//    double error = std::fabs(solution - y);
//    std::cout << "\nError: " << std::scientific << error << std::fixed << " (" << error << ")\n";
//    return 0;
//}

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

double f(double t, double y) {
    double f1 = std::exp(-t * t) * std::cos(0.8 * t);
    double f2 = std::exp(2 * t - 3);
    double g = 2 * (t - 2);
    return f1 * f2 - g * y;
}

double solution(double t) {
    // t1 — частина з cos і sin
    double t1 = std::exp(-t * t + 2.0 * t - 3.0) *
        (-25.0 / 58.0 * std::cos(0.8 * t)
            + 5.0 / 29.0 * std::sin(0.8 * t));

    // t2 — друга частина, пов'язана з початковими умовами
    double A = -10.0
        - (25.0 / 28.0) * std::cos(0.8)
        + (5.0 / 29.0) * std::sin(0.8);

    double t2 = (std::exp(-t * t + 4.0 * t) * A) / std::exp(5);

    return t1 - t2;
}

double rk4_step(double t, double y, double h) {
    double k1 = f(t, y);
    double k2 = f(t + h / 2.0, y + h * k1 / 2.0);
    double k3 = f(t + h / 2.0, y + h * k2 / 2.0);
    double k4 = f(t + h, y + h * k3);
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

int main() {
    double t0 = 1.0;
    double T = 6.0;
    double y = 10.0;
    double tau = 0.1;        // початковий крок
    double eps = 1e-5;
    double epsM = 1e-6;
    double e_max = 0.0;

    double t = t0;
    const int p = 4;

    std::cout << std::fixed << std::setprecision(10);

    double u0 = solution(t0);
    std::cout << "t=" << t
        << " y=" << y
        << " u(t)=" << u0
        << " |y-u(t)|=" << std::fabs(y - u0)
        << "\n";

    int steps = 0;

    while (true) {
        if (std::fabs(T - t) < epsM) break;
        if (t + tau > T) tau = T - t;

        double w = rk4_step(t, y, tau);
        double y_half = rk4_step(t, y, tau / 2.0);
        y_half = rk4_step(t + tau / 2.0, y_half, tau / 2.0);

        double E = std::fabs(y_half - w) / std::max(1.0, std::fabs(y_half));
        double tauH = std::min(5.0, std::max(0.1, 0.9 * std::pow(eps / (E + 1e-16), 1.0 / (p + 1))));

        if (E <= eps) {
            t += tau;
            y = y_half;

            double u_t = solution(t); //analitichniy
            double local_error = std::fabs(y - u_t);
            if (local_error > e_max) e_max = local_error;

            std::cout << "t=" << t
                << " y=" << y
                << " u(t)=" << u_t
                << " |y-u(t)|=" << local_error << "\n";

            tau *= tauH;
            steps++;
        }
        else {
            tau *= tauH;
            if (tau < epsM) {
                std::cerr << "Step too small, aborting.\n";
                break;
            }
        }

        if (t >= T) break;
    }

    std::cout << "||e||_max = " << e_max << "\n";

    double uT = solution(T);
    double error = std::fabs(uT - y);
    std::cout << "Global error: " << std::scientific << error << std::fixed << "\n";
    std::cout << "Number of steps: " << steps << "\n";

    return 0;
}

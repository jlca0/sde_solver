#pragma once
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include <stdexcept>

using namespace std;

class BrownianPath {
    public:
        vector<double> camino;

        BrownianPath(double T, int pasos) {
            double dt = T / pasos;
            double sqrtdt = sqrt(dt);
            random_device rd;
            mt19937 gen(rd());
            normal_distribution<double> normal(0.0, 1.0);

            camino.resize(pasos + 1);
            camino[0] = 0.0;
            for (int i = 0; i < pasos; ++i) {
                double W = normal(gen) * sqrtdt;
                camino[i+1] = camino[i] + W;
            }
        }
};

class BrownianBridge {
    public:
        vector<double> camino;

        BrownianBridge(double T, const BrownianPath& W) {
            int pasos = W.camino.size() - 1;
            double dt = T / pasos;

            camino.resize(pasos + 1);
            for (int i = 0; i <= pasos; ++i) {
                double t = i * dt;
                camino[i] = W.camino[i] - (t / T) * W.camino[pasos];
            }
        }
};

class ExactGeomPath {
    private:
        vector<double> camino;

    public:

        ExactGeomPath(double X0, double T, double mu, double sigma, const BrownianPath& W) {
            int pasos = W.camino.size() - 1;
            double dt = T / pasos;
            
            camino.resize(pasos + 1);
            camino[0] = X0;

            for (int i = 0; i < pasos; ++i) {
                double dW = W.camino[i+1] - W.camino[i];
                camino[i+1] = camino[i] * exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);
            }
        }

        vector<double> get_trayectoria() const {
            return camino;
        }
};

class EulerMaruyama {
    private:
        vector<double> trayectoria;

    public:
        EulerMaruyama(int pasos) {
            trayectoria.resize(pasos + 1);
        }

        vector<double> get_trayectoria() const {
            return trayectoria;
        }

        void solve(double X0, double T, 
                function<double(double,double)> mu, 
                function<double(double,double)> sigma,
                const BrownianPath& W) {

            int pasos_camino = W.camino.size() - 1;
            int pasos_trayectoria = trayectoria.size() - 1;

            if (pasos_camino % pasos_trayectoria != 0) {
                throw runtime_error("El tamaño del camino browniano no es múltiplo del número de pasos.");
            }

            int inc = pasos_camino / pasos_trayectoria;
            double dt = T / pasos_trayectoria;
            
            trayectoria[0] = X0;

            for (int i = 0; i < pasos_trayectoria; ++i) {
                double t_i = i * dt;
                double dW = W.camino[(i + 1) * inc] - W.camino[i * inc];
                
                trayectoria[i+1] = trayectoria[i] + mu(t_i, trayectoria[i]) * dt + sigma(t_i, trayectoria[i]) * dW;
            }
        }
};

class Milstein {
    private:
        vector<double> trayectoria;

    public:
        Milstein(int pasos) {
            trayectoria.resize(pasos + 1);
        }

        vector<double> get_trayectoria() const {
            return trayectoria;
        }

        void solve(double X0, 
                   double T, 
                   function<double(double,double)> mu, 
                   function<double(double,double)> sigma,
                   function<double(double,double)> sigma_x,
                   const BrownianPath& W) {

            int pasos_camino = W.camino.size() - 1;
            int pasos_trayectoria = trayectoria.size() - 1;

            if (pasos_camino % pasos_trayectoria != 0) {
                throw runtime_error("El tamaño del camino browniano no es múltiplo del número de pasos.");
            }

            int inc = pasos_camino / pasos_trayectoria;
            double dt = T / pasos_trayectoria;
            
            trayectoria[0] = X0;

            for (int i = 0; i < pasos_trayectoria; ++i) {
                double t_i = i * dt;
                double dW = W.camino[(i + 1) * inc] - W.camino[i * inc];
                
                trayectoria[i+1] = trayectoria[i] + mu(t_i, trayectoria[i]) * dt + sigma(t_i, trayectoria[i]) * dW + 0.5 * sigma(t_i, trayectoria[i]) * sigma_x(t_i, trayectoria[i]) * (dW*dW - dt);
            }
        }
};
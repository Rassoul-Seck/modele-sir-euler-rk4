"""Modèle épidémique SIR résolu par Euler explicite et Runge-Kutta 4.

Système :
    dS/dt = -beta * S * I / N
    dI/dt =  beta * S * I / N - gamma * I
    dR/dt =  gamma * I
"""
import numpy as np
from scipy.integrate import solve_ivp


def sir(t, y, beta, gamma):
    """Second membre du système SIR."""
    S, I, R = y
    N = S + I + R
    infection = beta * S * I / N
    return np.array([-infection, infection - gamma * I, gamma * I])


def euler(f, y0, t0, tf, h, args=()):
    """Méthode d'Euler explicite (ordre 1)."""
    n = int(round((tf - t0) / h))
    t = t0 + h * np.arange(n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for k in range(n):
        y[k + 1] = y[k] + h * f(t[k], y[k], *args)
    return t, y


def rk4(f, y0, t0, tf, h, args=()):
    """Méthode de Runge-Kutta classique (ordre 4)."""
    n = int(round((tf - t0) / h))
    t = t0 + h * np.arange(n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for k in range(n):
        k1 = f(t[k], y[k], *args)
        k2 = f(t[k] + h / 2, y[k] + h / 2 * k1, *args)
        k3 = f(t[k] + h / 2, y[k] + h / 2 * k2, *args)
        k4 = f(t[k] + h, y[k] + h * k3, *args)
        y[k + 1] = y[k] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return t, y


def reference(y0, t0, tf, args):
    """Solution de référence très précise (SIR n'a pas de formule exacte simple)."""
    return solve_ivp(sir, (t0, tf), y0, args=args, method="DOP853",
                     rtol=1e-12, atol=1e-12, dense_output=True)

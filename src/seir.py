"""Modèle épidémique SEIR (avec période d'incubation).

Système :
    dS/dt = -beta * S * I / N
    dE/dt =  beta * S * I / N - sigma * E
    dI/dt =  sigma * E - gamma * I
    dR/dt =  gamma * I

E : exposés (infectés mais pas encore contagieux), 1/sigma = durée d'incubation.
"""
import numpy as np


def seir(t, y, beta, sigma, gamma):
    """Second membre du système SEIR."""
    S, E, I, R = y
    N = S + E + I + R
    infection = beta * S * I / N
    return np.array([-infection,
                     infection - sigma * E,
                     sigma * E - gamma * I,
                     gamma * I])

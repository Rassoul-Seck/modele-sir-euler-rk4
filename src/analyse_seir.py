"""Comparaison SIR / SEIR, effet de l'incubation et de R0 sur le pic épidémique."""
import os
import numpy as np
import matplotlib.pyplot as plt
from sir import sir, rk4
from seir import seir

N = 1000.0
gamma = 0.1
t0, tf, h = 0.0, 250.0, 0.1
os.makedirs("results", exist_ok=True)


def simuler_sir(beta, I0=1.0):
    y0 = np.array([N - I0, I0, 0.0])
    return rk4(sir, y0, t0, tf, h, (beta, gamma))


def simuler_seir(beta, sigma, E0=1.0):
    y0 = np.array([N - E0, E0, 0.0, 0.0])
    return rk4(seir, y0, t0, tf, h, (beta, sigma, gamma))


def pic(t, infectes):
    k = np.argmax(infectes)
    return t[k], infectes[k]


# --- Vérification : la population totale reste constante ---
t, y = simuler_seir(0.3, 0.2)
print(f"Écart max sur N (SEIR) : {np.max(np.abs(y.sum(axis=1) - N)):.2e}")

# --- Figure 1 : SIR vs SEIR (mêmes beta et gamma) ---
beta, sigma = 0.3, 0.2                     # incubation moyenne = 5 jours
t1, y1 = simuler_sir(beta)
t2, y2 = simuler_seir(beta, sigma)
tp1, ip1 = pic(t1, y1[:, 1])
tp2, ip2 = pic(t2, y2[:, 2])
print(f"SIR  : pic de {ip1:.0f} infectés au jour {tp1:.0f}")
print(f"SEIR : pic de {ip2:.0f} infectés au jour {tp2:.0f}")

plt.figure(figsize=(8, 5))
plt.plot(t1, y1[:, 1], color="tab:red", label="Infectés (SIR)")
plt.plot(t2, y2[:, 2], color="tab:purple", label="Infectés (SEIR)")
plt.plot(t2, y2[:, 1], "--", color="tab:orange", label="Exposés (SEIR)")
plt.xlabel("Temps (jours)")
plt.ylabel("Nombre d'individus")
plt.title(f"SIR vs SEIR (β = {beta}, γ = {gamma}, σ = {sigma})")
plt.xlim(0, 160)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("results/sir_vs_seir.png", dpi=150)
plt.close()

# --- Figure 2 : effet de la durée d'incubation ---
plt.figure(figsize=(8, 5))
for sigma in [0.1, 0.2, 0.5, 2.0]:
    t, y = simuler_seir(beta, sigma)
    tp, ip = pic(t, y[:, 2])
    print(f"sigma = {sigma:<4} (incubation {1 / sigma:4.1f} j) : "
          f"pic {ip:.0f} au jour {tp:.0f}")
    plt.plot(t, y[:, 2], label=f"σ = {sigma} (incubation {1 / sigma:.1f} j)")
plt.xlabel("Temps (jours)")
plt.ylabel("Infectés")
plt.title("Effet de l'incubation sur l'épidémie (SEIR)")
plt.xlim(0, 200)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("results/effet_incubation.png", dpi=150)
plt.close()

# --- Figure 3 : pic épidémique en fonction de R0 (SIR) ---
R0s = np.linspace(1.2, 6, 25)
pics_num, pics_th = [], []
for R0 in R0s:
    t, y = simuler_sir(R0 * gamma)
    pics_num.append(y[:, 1].max())
    # Formule exacte (intégrale première du SIR) : I_max = I0 + S0 - (N/R0)(1 + ln(R0 S0 / N))
    S0, I0 = N - 1.0, 1.0
    pics_th.append(I0 + S0 - (N / R0) * (1 + np.log(R0 * S0 / N)))

ecart = np.max(np.abs(np.array(pics_num) - np.array(pics_th)))
print(f"Écart max entre pic numérique et formule : {ecart:.3f} individus")

plt.figure(figsize=(7, 5))
plt.plot(R0s, np.array(pics_num) / N * 100, "o", label="Simulation (RK4)")
plt.plot(R0s, np.array(pics_th) / N * 100, "-", label="Formule théorique")
plt.xlabel("R0 = β/γ")
plt.ylabel("Pic d'infectés (% de la population)")
plt.title("Pic épidémique en fonction de R0 (SIR)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("results/pic_vs_R0.png", dpi=150)
plt.close()

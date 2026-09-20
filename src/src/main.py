"""Simulation SIR : trajectoires et étude de convergence."""
import os
import numpy as np
import matplotlib.pyplot as plt
from sir import sir, euler, rk4, reference

# --- Paramètres du modèle ---
N = 1000.0
I0 = 1.0
y0 = np.array([N - I0, I0, 0.0])
beta, gamma = 0.3, 0.1          # R0 = beta / gamma = 3
t0, tf = 0.0, 160.0
args = (beta, gamma)

os.makedirs("results", exist_ok=True)
ref = reference(y0, t0, tf, args)

# --- Figure 1 : trajectoires ---
t, y = rk4(sir, y0, t0, tf, 0.1, args)
plt.figure(figsize=(8, 5))
for i, (nom, c) in enumerate(zip(["Susceptibles", "Infectés", "Guéris"],
                                 ["tab:blue", "tab:red", "tab:green"])):
    plt.plot(t, y[:, i], label=nom, color=c)
plt.xlabel("Temps (jours)")
plt.ylabel("Nombre d'individus")
plt.title(f"Modèle SIR (RK4), R0 = {beta / gamma:.1f}")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("results/trajectoires.png", dpi=150)
plt.close()

# --- Figure 2 : erreur en fonction du pas (convergence) ---
pas = [4, 2, 1, 0.5, 0.25, 0.125]
err_e, err_r = [], []
for h in pas:
    for methode, liste in [(euler, err_e), (rk4, err_r)]:
        tt, yy = methode(sir, y0, t0, tf, h, args)
        exact = ref.sol(tt).T
        liste.append(np.max(np.abs(yy - exact)))

ordre_e = np.polyfit(np.log(pas), np.log(err_e), 1)[0]
ordre_r = np.polyfit(np.log(pas), np.log(err_r), 1)[0]
print(f"Ordre observé Euler : {ordre_e:.2f} (théorie : 1)")
print(f"Ordre observé RK4   : {ordre_r:.2f} (théorie : 4)")

plt.figure(figsize=(7, 5))
plt.loglog(pas, err_e, "o-", label=f"Euler (pente ≈ {ordre_e:.1f})")
plt.loglog(pas, err_r, "s-", label=f"RK4 (pente ≈ {ordre_r:.1f})")
plt.xlabel("Pas h")
plt.ylabel("Erreur maximale")
plt.title("Convergence des méthodes")
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.savefig("results/convergence.png", dpi=150)
plt.close()

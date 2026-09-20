# Modèle épidémique SIR : Euler vs Runge-Kutta 4

Simulation numérique du modèle SIR (Susceptibles, Infectés, Guéris) et comparaison
de deux méthodes de résolution d'EDO. Projet réalisé dans le cadre de la Licence
Modélisation Mathématique, Analyse et Simulation Numériques (UNCHK).

## Le modèle

Population totale constante N = S + I + R :

    dS/dt = -β·S·I/N
    dI/dt =  β·S·I/N - γ·I
    dR/dt =  γ·I

- β : taux de transmission
- γ : taux de guérison
- R0 = β/γ : nombre de reproduction de base (l'épidémie croît si R0 > 1)

Paramètres utilisés : N = 1000, 1 infecté au départ, β = 0.3, γ = 0.1 (R0 = 3).

## Méthodes

- **Euler explicite** (ordre 1)
- **Runge-Kutta 4** (ordre 4)
- **Référence** : `scipy.integrate.solve_ivp` (DOP853, tolérance 1e-12), car le
  système n'a pas de solution exacte simple.

## Installation et lancement

    pip install -r requirements.txt
    python src/main.py

## Résultats

Trajectoires (RK4) :

![Trajectoires](results/trajectoires.png)

Erreur en fonction du pas h (échelle log-log) :

![Convergence](convergence.png)

Ordres de convergence observés : ≈ 0.96 pour Euler et ≈ 3.84 pour RK4, ce qui
est cohérent avec la théorie (1 et 4).

## Limites

- Modèle simple : population homogène et constante, pas de naissances, de décès
  ni de vaccination.
- Les paramètres sont fictifs, sans calage sur des données réelles.

## Pistes d'amélioration

- Ajouter la vaccination ou le modèle SEIR (période d'incubation).
- Estimer β et γ à partir de données réelles (moindres carrés).
- Tester la sensibilité du pic épidémique à R0.

## Références

- Kermack & McKendrick (1927), modèle SIR.
- Cours d'analyse numérique, UNCHK.

## Licence

MIT

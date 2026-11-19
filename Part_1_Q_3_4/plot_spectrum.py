import numpy as np
import matplotlib.pyplot as plt

Hs = 6.0      # significant wave height Hm0 [m]
Tp = 8.8      # peak period [s], deep-water steepness s = 2*pi*Hs/(g*Tp^2) = 0.05
gamma = 3.3   # JONSWAP peak enhancement (1 = Pierson-Moskowitz)

fp = 1.0 / Tp
f = np.linspace(0.01, 0.5, 2000)   # frequency axis [Hz]

# Pierson-Moskowitz, written with Hs and Tp
E_pm = (5 / 16) * Hs**2 * fp**4 * f**-5 * np.exp(-1.25 * (fp / f)**4)

# JONSWAP = PM x normalisation x peak enhancement
sigma = np.where(f <= fp, 0.07, 0.09)
A_gamma = 1 - 0.287 * np.log(gamma)
E_js = A_gamma * E_pm * gamma**np.exp(-(f - fp)**2 / (2 * sigma**2 * fp**2))



plt.plot(f, E_pm, label="Pierson-Moskowitz")
plt.plot(f, E_js, label=f"JONSWAP, γ = {gamma}")
plt.xlabel("f [Hz]")
plt.ylabel("E(f) [m²/Hz]")
plt.title(f"Variance density spectrum, Hs = {Hs} m, Tp = {Tp} s")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("spectrum.png", dpi=150)
plt.show()

# spectral periods and steepness
g, h = 9.81, 20.0                       # water depth [m]
df = f[1] - f[0]
m = {n: np.sum(f**n * E_js) * df for n in (-1, 0, 1, 2)}
Tm01 = m[0] / m[1]
Tm02 = np.sqrt(m[0] / m[2])
Tm10 = m[-1] / m[0]
print(f"Tm01 = {Tm01:.2f} s, Tm02 = {Tm02:.2f} s, Tm-1,0 = {Tm10:.2f} s_o_p, Tp = {Tp} s")


def wavelength(T):
    "Calculate wavelength from period using linear wave theory"
    return g * T**2 / (2 * np.pi) 

L_0 = wavelength(Tp)
L_0_tm01 = wavelength(Tm01)
L_0_tm02 = wavelength(Tm02)
L_0_tm10 = wavelength(Tm10)
s_L_0 = Hs / L_0
s_L_0_tm01 = Hs / L_0_tm01
s_L_0_tm02 = Hs / L_0_tm02
s_L_0_tm10 = Hs / L_0_tm10
print(f"Deep-water wavelength L0 = {L_0:.2f} m, steepness s = {s_L_0:.3f}")
print(f"Tm01 wavelength L0 = {L_0_tm01:.2f} m, steepness s = {s_L_0_tm01:.3f}")
print(f"Tm02 wavelength L0 = {L_0_tm02:.2f} m, steepness s = {s_L_0_tm02:.3f}")
print(f"Tm10 wavelength L0 = {L_0_tm10:.2f} m, steepness s = {s_L_0_tm10:.3f}")

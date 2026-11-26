import numpy as np

def jonswap(Tp, Hs, f):
    """
    Calculate the JONSWAP spectrum for a given period Tp, significant wave height Hs 
and frequency array f.
    
    Parameters:
    Tp : float
        Peak period [s].
    Hs : float
        Significant wave height [m].
    f : numpy.ndarray
        Frequency array [Hz].
        
    Returns:
    E :JONSWAP spectrum values [m²/Hz].
    """
    g = 9.81  # acceleration due to gravity [m/s²]
    gamma = 3.3  # JONSWAP peak enhancement factor

    fp = 1.0 / Tp

    # Pierson-Moskowitz spectrum
    E_pm = (5 / 16) * Hs**2 * fp**4 * f**-5 * np.exp(-1.25 * (fp / f)**4)

    # JONSWAP spectrum
    sigma = np.where(f <= fp, 0.07, 0.09)
    A_gamma = 1 - 0.287 * np.log(gamma)
    E_js = A_gamma * E_pm * gamma**np.exp(-(f - fp)**2 / (2 * sigma**2 * fp**2))

    return E_js

D = 3 * 3600 #record length [s] (3 hours)
df = 1 / D  # frequency resolution [Hz]

h = 20.0 # water depth [m]
g = 9.81 # gravity [m/s^2]
Hs = 6.0      # significant wave height Hm0 [m]
Tp = 8.8      # peak period [s], deep-water steepness s = 2*pi*Hs/(g*Tp^2) = 0.05 (5% of the steepness limit)

f_min = 0.5 * Tp**(-1)  # minimum frequency [Hz]
f_max = 4.5 * Tp**(-1)   # maximum frequency [Hz]

f = np.arange(f_min, f_max, df) + df / 2   # new frequency array [Hz]
E = jonswap(Tp, Hs, f=f) #  spectrum from Q3, evaluated on f
a = np.sqrt(2 * E * df) # amplitude spectrum [m]
phi = np.random.default_rng(2026).uniform(0, 2 * np.pi, f.size) # random phase spectrum [rad]

dt = 0.1
t = np.arange(0, D, dt)
eta = np.zeros_like(t)

# chunk over time so the (components x times) matrix stays small
chunk = 5000
for s in range(0, t.size, chunk):
    tt = t[s:s + chunk]
    eta[s:s + chunk] = (a[:, None] * np.cos(2*np.pi*f[:, None]*tt[None, :] + phi[:, None])).sum(axis=0)

omega = 2 * np.pi * f # angular frequency [rad/s]

#using dispersion relation omega^2 = g * k * np.tanh(k * h) to get wavenumber k from omega  

k = omega**2 / g                       # deep-water first guess
for i in range(100):                   # fixed-point iteration, all components together
    k = omega**2 / (g * np.tanh(k * h))

u = np.zeros_like(t) # horizontal velocity at z = 0 [m/s], u = a*omega*cosh(k(z+h))/sinh(kh)*cos(...)
for s in range(0, t.size, chunk):
    tt = t[s:s + chunk]
    u[s:s + chunk] = (a[:, None] * omega[:, None] * np.cosh(k[:, None]*h) / np.sinh(k[:, None]*h)
                      * np.cos(2*np.pi*f[:, None]*tt[None, :] + phi[:, None])).sum(axis=0)

print("Hm0 =", 4 * eta.std()) 
print("mean =", eta.mean())
print(np.max(np.abs(omega**2 - g*k*np.tanh(k*h))))

# Everything below (Q4b and all figures) only runs when this file is run directly,
# not when another script (e.g. Question 6) imports the Q4a time series from it.
if __name__ == "__main__":

    #---------4b----------------

    #using a coarser frequency resolution to get a shorter time series

    N2 = 300 #coarser frequency resolution
    df2 = (f_max - f_min) / N2 #coarser frequency resolution [Hz]
    f2 = np.arange(f_min, f_max, df2) + df2 / 2   # new frequency array [Hz]
    D2 = 1 / df2 # new duration [s]

    E2 = jonswap(Tp, Hs, f=f2) #  spectrum from Q3, evaluated on f
    a2 = np.sqrt(2 * E2 * df2) # amplitude spectrum [m]
    phi2 = np.random.default_rng(2026).uniform(0, 2 * np.pi, f2.size) # random phase spectrum [rad]

    omega2 = 2 * np.pi * f2 # angular frequency [rad/s]
    k2 = omega2**2 / g                       # deep-water first guess
    for i in range(100):                   # fixed-point iteration, all components together 
        k2 = omega2**2 / (g * np.tanh(k2 * h))


    #code from Forristall (2000) to calculate the second-order sum-frequency kernel Bplus

    def Bplus(ki, wi, kj, wj, h):
        """Sum-frequency kernel, Sharma & Dean (1981), long-crested,
        notation of Forristall (2000).
        eta2 = sum_i sum_j a_i a_j Bplus_ij cos(theta_i + theta_j)"""

        Ri, Rj = wi**2 / g, wj**2 / g
        si, sj = np.sqrt(Ri), np.sqrt(Rj)
        kp = ki + kj
        Dp = ((si + sj) * (si * (kj**2 - Rj**2) + sj * (ki**2 - Ri**2))
              + 2 * (si + sj)**2 * (ki * kj - Ri * Rj)) \
             / ((si + sj)**2 - kp * np.tanh(kp * h))
        return 0.25 * ((Dp - (ki * kj - Ri * Rj)) / np.sqrt(Ri * Rj) + (Ri + Rj))

    # test 1: i = j, deep water -> Bplus = k/2
    kt = 0.06; ht = 1e4; wt = np.sqrt(g * kt)
    print(Bplus(kt, wt, kt, wt, ht), 0.5 * kt)

    # test 2: i = j, h = 20 m -> finite-depth Stokes second order
    ht = 20.0; wt = np.sqrt(g * kt * np.tanh(kt * ht))
    print(Bplus(kt, wt, kt, wt, ht),
          kt/4 * np.cosh(kt*ht) * (2 + np.cosh(2*kt*ht)) / np.sinh(kt*ht)**3)

    #velocity kernel: x-derivative of the Sharma and Dean (1981) second-order potential, D+ in the notation of Forristall (2000)

    def Dplus(ki, wi, kj, wj, h):
        Ri, Rj = wi**2 / g, wj**2 / g
        si, sj = np.sqrt(Ri), np.sqrt(Rj)
        kp = ki + kj
        return ((si + sj) * (si * (kj**2 - Rj**2) + sj * (ki**2 - Ri**2))
                + 2 * (si + sj)**2 * (ki * kj - Ri * Rj)) \
               / ((si + sj)**2 - kp * np.tanh(kp * h))

    def Gplus(ki, wi, kj, wj, h, z=0.0):
        kp, wp = ki + kj, wi + wj
        cosh_ratio = np.exp(kp*z) * (1 + np.exp(-2*kp*(z + h))) / (1 + np.exp(-2*kp*h))
        return 0.25 * g**2 * Dplus(ki, wi, kj, wj, h) * kp / (wi * wj * wp) * cosh_ratio

    # test 3: i = j, h = 20 m -> Stokes second order u2 = (3/4) w k a^2 cosh 2k(z+h) / sinh^4 kh (at z = 0)
    print(Gplus(kt, wt, kt, wt, ht), 0.75*wt*kt*np.cosh(2*kt*ht)/np.sinh(kt*ht)**4)

    #eta1, u1 (first order) and eta2, u2 (second-order sum-frequency) on the same a2, f2, phi2 grid
    z = 0.0
    t = np.arange(0, D2, dt)
    eta1, u1, eta2, u2 = (np.zeros_like(t) for _ in range(4))
    for s in range(0, t.size, 2000):
        tt = t[s:s + 2000][None, :]
        c = np.cos(2*np.pi*f2[:, None]*tt + phi2[:, None])
        eta1[s:s + 2000] = (a2[:, None] * c).sum(axis=0)
        u1[s:s + 2000] = (a2[:, None] * omega2[:, None] * np.cosh(k2[:, None]*(z + h))
                          / np.sinh(k2[:, None]*h) * c).sum(axis=0)

    for i in range(f2.size):
        Bp = Bplus(k2[i], omega2[i], k2, omega2, h)
        Gp = Gplus(k2[i], omega2[i], k2, omega2, h, z)
        ph = (phi2[i] + phi2)[:, None]
        fr = (f2[i] + f2)[:, None]
        for s in range(0, t.size, 2000):
            c = np.cos(2*np.pi*fr*t[s:s + 2000][None, :] + ph)
            eta2[s:s + 2000] += ((a2[i]*a2*Bp)[:, None] * c).sum(axis=0)
            u2[s:s + 2000] += ((a2[i]*a2*Gp)[:, None] * c).sum(axis=0)

    eta_total = eta1 + eta2
    u_total = u1 + u2

    #---------plots for the report (Q3 and Q4a, 4b), saved next to this script----------------

    import os
    import matplotlib.pyplot as plt

    out_dir = os.path.dirname(os.path.abspath(__file__))

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, name), dpi=200)
        print("saved", name)

    def welch_psd(x, dt, nseg):
        """One-sided spectral estimate, Hann window, 50 % overlap (Lecture 5 slide 8)."""
        win = np.hanning(nseg)
        segs = [x[i:i + nseg] for i in range(0, x.size - nseg + 1, nseg // 2)]
        P = np.mean([np.abs(np.fft.rfft((s_ - s_.mean()) * win))**2 for s_ in segs], axis=0)
        return np.fft.rfftfreq(nseg, dt), 2 * dt * P / np.sum(win**2)

    def skewness(x):
        return np.mean((x - x.mean())**3) / x.std()**3

    def crests_troughs(x):
        """Zero up-crossing analysis: crest height and trough depth of every wave."""
        up = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
        c = np.array([x[i0:i1 + 1].max() for i0, i1 in zip(up[:-1], up[1:])])
        tr = np.array([-x[i0:i1 + 1].min() for i0, i1 in zip(up[:-1], up[1:])])
        return c, tr

    # --- Q3: spectrum, spectral periods and steepness ---
    fq = np.linspace(0.01, 0.5, 2000)                      # same frequency axis as plot_spectrum.py
    fp = 1 / Tp
    E_q = jonswap(Tp, Hs, fq)
    E_pm = (5 / 16) * Hs**2 * fp**4 * fq**-5 * np.exp(-1.25 * (fp / fq)**4)
    dfq = fq[1] - fq[0]
    m = {n: np.sum(fq**n * E_q) * dfq for n in (-1, 0, 1, 2)}
    Tm01, Tm02, Tm10 = m[0] / m[1], np.sqrt(m[0] / m[2]), m[-1] / m[0]
    print(f"Q3: m0 = {m[0]:.3f} m2, Hm0 = {4*np.sqrt(m[0]):.2f} m, Tm01 = {Tm01:.2f} s, "
          f"Tm02 = {Tm02:.2f} s, Tm-1,0 = {Tm10:.2f} s")
    for name, T in (("Tp", Tp), ("Tm-1,0", Tm10), ("Tm01", Tm01), ("Tm02", Tm02)):
        print(f"   s0 with {name}: {2*np.pi*Hs/(g*T**2):.3f}")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(fq, E_pm, color="0.6", label="Pierson-Moskowitz")
    ax.plot(fq, E_q, "C0", lw=2, label=f"JONSWAP, γ = 3.3")
    for T, lab, col in ((Tp, "$T_p$", "k"), (Tm10, "$T_{m-1,0}$", "C1"), (Tm01, "$T_{m01}$", "C2"), (Tm02, "$T_{m02}$", "C3")):
        ax.axvline(1 / T, color=col, ls="--", lw=1, label=f"1/{lab} = {1/T:.3f} Hz")
    ax.set_xlabel("f [Hz]"); ax.set_ylabel("E(f) [m²/Hz]"); ax.set_xlim(0, 0.5)
    ax.set_title(f"Variance density spectrum, $H_{{m0}}$ = {Hs:g} m, $T_p$ = {Tp:g} s")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    save(fig, "fig_q3_spectrum.png")

    # --- Q4a: linear series on the fine grid (3 h record) ---
    tA = np.arange(eta.size) * dt                          # time axis of eta and u (t is reused for 4b)
    fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    w4a = (tA >= 1000) & (tA <= 1120)
    ax[0].plot(tA[w4a], eta[w4a], "C0"); ax[0].set_ylabel("η [m]")
    ax[1].plot(tA[w4a], u[w4a], "C1"); ax[1].set_ylabel("u(z = 0) [m/s]"); ax[1].set_xlabel("t [s]")
    for a_ in ax:
        a_.grid(alpha=0.3); a_.axhline(0, color="0.5", lw=0.5)
    ax[0].set_title("Q4a: linear random-phase realisation (120 s of the 3 h record)")
    save(fig, "fig_q4a_timeseries.png")

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    fw, P = welch_psd(eta, dt, 4096)
    ax[0].semilogy(fw, P, "C0", lw=1, label="time series (Welch)")
    ax[0].semilogy(fq, E_q, "k--", lw=1, label="JONSWAP target")
    ax[0].set_xlim(0, 0.6); ax[0].set_ylim(1e-4, 50); ax[0].set_xlabel("f [Hz]"); ax[0].set_ylabel("E(f) [m²/Hz]")
    ax[0].set_title("Surface elevation"); ax[0].legend(fontsize=8)
    fw, P = welch_psd(u, dt, 4096)
    Hu = omega * np.cosh(k * h) / np.sinh(k * h)
    ax[1].semilogy(fw, P, "C1", lw=1, label="time series (Welch)")
    ax[1].semilogy(f, Hu**2 * E, "k--", lw=1, label="$|H_u|^2 E$")
    ax[1].set_xlim(0, 0.6); ax[1].set_ylim(1e-4, 50); ax[1].set_xlabel("f [Hz]"); ax[1].set_ylabel("$S_u$(f) [(m/s)²/Hz]")
    ax[1].set_title("Velocity at z = 0"); ax[1].legend(fontsize=8)
    sig = eta.std()
    xs = np.linspace(-4, 4, 200)
    ax[2].hist(eta / sig, bins=80, density=True, color="C0", alpha=0.6, label="time series")
    ax[2].plot(xs, np.exp(-xs**2 / 2) / np.sqrt(2 * np.pi), "k--", lw=1, label="Gaussian")
    ax[2].set_xlabel("η/σ"); ax[2].set_ylabel("pdf"); ax[2].set_title(f"Distribution, skewness = {skewness(eta):.3f}")
    ax[2].legend(fontsize=8)
    for a_ in ax:
        a_.grid(alpha=0.3)
    save(fig, "fig_q4a_verification.png")

    # --- Q4b: first order vs first + second order (coarse grid, 660 s record) ---
    i_max = np.argmax(eta_total)
    t0 = min(max(t[i_max] - 50, 0), t[-1] - 100)        # 100 s window around the highest crest, inside the record
    w4b = (t >= t0) & (t <= t0 + 100)
    fig, ax = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True)
    ax[0].plot(t[w4b], eta1[w4b], "C0", label="η₁ (linear)")
    ax[0].plot(t[w4b], eta_total[w4b], "k", lw=1.2, label="η₁ + η₂⁺")
    ax[0].plot(t[w4b], eta2[w4b], "C3", lw=0.8, label="η₂⁺ (super-harmonic)")
    ax[0].set_ylabel("η [m]"); ax[0].legend(fontsize=8, loc="upper right")
    ax[1].plot(t[w4b], u1[w4b], "C1", label="u₁ (linear)")
    ax[1].plot(t[w4b], u_total[w4b], "k", lw=1.2, label="u₁ + u₂⁺")
    ax[1].plot(t[w4b], u2[w4b], "C3", lw=0.8, label="u₂⁺ (super-harmonic)")
    ax[1].set_ylabel("u(z = 0) [m/s]"); ax[1].set_xlabel("t [s]"); ax[1].legend(fontsize=8, loc="upper right")
    for a_ in ax:
        a_.grid(alpha=0.3); a_.axhline(0, color="0.5", lw=0.5)
    ax[0].set_title("Q4b: effect of the second-order super-harmonic bound waves (around the highest crest)")
    save(fig, "fig_q4b_timeseries.png")

    c1, tr1 = crests_troughs(eta1)
    c2, tr2 = crests_troughs(eta_total)
    print(f"Q4b: skewness eta1 = {skewness(eta1):.3f}, eta1+eta2 = {skewness(eta_total):.3f}; "
          f"skewness u1 = {skewness(u1):.3f}, u1+u2 = {skewness(u_total):.3f}")
    print(f"Q4b: waves = {c1.size}, max crest {c1.max():.2f} -> {c2.max():.2f} m, "
          f"max trough {tr1.max():.2f} -> {tr2.max():.2f} m, "
          f"mean crest/trough {c1.mean()/tr1.mean():.3f} -> {c2.mean()/tr2.mean():.3f}")

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    fw, P1 = welch_psd(eta1, dt, 2048)
    fw, P2 = welch_psd(eta_total, dt, 2048)
    ax[0].semilogy(fw, P1, "C0", lw=1, label="η₁")
    ax[0].semilogy(fw, P2, "k", lw=1, label="η₁ + η₂⁺")
    ax[0].semilogy(fq, E_q, "k--", lw=0.8, label="JONSWAP target")
    ax[0].set_xlim(0, 1.0); ax[0].set_ylim(1e-4, 50); ax[0].set_xlabel("f [Hz]"); ax[0].set_ylabel("E(f) [m²/Hz]")
    ax[0].set_title("Spectrum: energy at sum frequencies"); ax[0].legend(fontsize=8)
    sig1 = eta1.std()
    ax[1].hist(eta1 / sig1, bins=50, density=True, histtype="step", color="C0", label=f"η₁, skew {skewness(eta1):.2f}")
    ax[1].hist(eta_total / sig1, bins=50, density=True, histtype="step", color="k", label=f"η₁+η₂⁺, skew {skewness(eta_total):.2f}")
    ax[1].plot(xs, np.exp(-xs**2 / 2) / np.sqrt(2 * np.pi), "k--", lw=0.8, label="Gaussian")
    ax[1].set_xlabel("η/σ₁"); ax[1].set_ylabel("pdf"); ax[1].set_title("Distribution of η"); ax[1].legend(fontsize=8)
    xr = np.linspace(0, 1.2, 200)
    for v, col, lab in ((c1, "C0", "crests η₁"), (c2, "k", "crests η₁+η₂⁺"), (tr2, "C3", "troughs η₁+η₂⁺")):
        vs = np.sort(v)[::-1] / Hs
        ax[2].semilogy(vs, np.arange(1, vs.size + 1) / vs.size, ".", color=col, ms=4, label=lab)
    ax[2].semilogy(xr, np.exp(-8 * xr**2), "k--", lw=0.8, label="Rayleigh")
    ax[2].set_xlim(0, 1.2); ax[2].set_ylim(5e-3, 1.1); ax[2].set_xlabel("η_c/H_s, η_t/H_s"); ax[2].set_ylabel("exceedance probability")
    ax[2].set_title("Crests and troughs"); ax[2].legend(fontsize=8)
    for a_ in ax:
        a_.grid(alpha=0.3)
    save(fig, "fig_q4b_verification.png")

    plt.show()

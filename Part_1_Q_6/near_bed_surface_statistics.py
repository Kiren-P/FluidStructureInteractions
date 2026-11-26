"""
Part I, Question 6: statistics of
  (a) the horizontal orbital velocity near the bed, and the maximum horizontal
      pressure gradient near the bed, from the time series;
  (b) the horizontal orbital velocity near the surface.

Linear wave theory. The wave field is NOT rebuilt here: the spectrum, the
frequency grid, the amplitudes, random phases, wavenumbers and the time series
eta(t) and u(z = 0, t) are imported from the Question 4a code
(Part_1_Q_3_4/random_phase_timeseries.py), so Question 6 uses exactly the same
realisation (JONSWAP, Hm0 = 6 m, Tp = 8.8 s, h = 20 m, 3 h, seed 2026).
Only the extra quantities are computed here, by passing every component
through the linear transfer function of the quantity (Lecture 5, slide 15):

  u(z, t)       = sum a_i w_i   cosh k_i(z+h)/sinh k_i h  cos(w_i t + phi_i)
  du/dt(z, t)   = sum a_i w_i^2 cosh k_i(z+h)/sinh k_i h  sin(w_i t + phi_i)
  dp/dx(z, t)   = -rho du/dt            (linearised Euler, Lecture 4 slide 14)

near the bed z = -h (top of the wave boundary layer; potential flow, the
boundary layer itself is not modelled) and near the surface z = 0 (still
water level, as in Question 4).

Statistics: rms value from the series and from the spectrum
(r_rms = sqrt(m0,r) = sqrt(sum |H_r|^2 E df)), mean zero-crossing period,
spectral bandwidth, the maximum and minimum of every zero-crossing cycle
(the "amplitude" of each individual wave in u, du/dt), their exceedance
distribution compared with the Rayleigh distribution
P(r_hat > x) = exp(-x^2 / (2 r_rms^2)), significant values (mean of the
highest third), r_1/10 and the largest value in the record compared with the
expected maximum of N Rayleigh-distributed peaks.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

out_dir = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# Wave field of Question 4a (imported, not recomputed)
# Only the Q4a part of random_phase_timeseries.py runs on import; Q4b and the
# Q3/Q4 figures sit under `if __name__ == "__main__":` in that file.
# ----------------------------------------------------------------------------
sys.path.insert(0, os.path.join(out_dir, "..", "Part_1_Q_3_4"))
print("Importing the Question 4a time series ...")
import random_phase_timeseries as q4

g, h, Hs, Tp = q4.g, q4.h, q4.Hs, q4.Tp       # 9.81 m/s2, 20 m, 6 m, 8.8 s
D, df, dt = q4.D, q4.df, q4.dt                # 3 h, 1/D, 0.1 s
f, E, a, phi = q4.f, q4.E, q4.a, q4.phi       # frequencies, JONSWAP, amplitudes, phases
omega, k = q4.omega, q4.k                     # angular frequency, wavenumber (dispersion)
t, eta = q4.t, q4.eta                         # time and surface elevation of Q4a
us = q4.u                                     # u at z = 0 from Q4a = velocity near the surface
SEED = 2026                                   # seed used in Q4a
jonswap = q4.jonswap

rho = 1025.0                                  # density of sea water [kg/m^3]

# ----------------------------------------------------------------------------
# Extra linear transfer functions for Question 6 (amplitude per metre of eta)
# ----------------------------------------------------------------------------
H_ub = omega / np.sinh(k * h)                 # u at the bed,      z = -h
H_us = omega / np.tanh(k * h)                 # u at the surface,  z = 0 (as Q4a)
H_ab = omega**2 / np.sinh(k * h)              # du/dt at the bed
H_as = omega**2 / np.tanh(k * h)              # du/dt at the surface
H_dudz = omega * k                            # du/dz at z = 0 (for the extrapolation check)


def wavenumber(om):
    """Same fixed-point iteration as Q4a, for single frequencies (checks below)."""
    kk = om**2 / g
    for _ in range(100):
        kk = om**2 / (g * np.tanh(kk * h))
    return kk


# ----------------------------------------------------------------------------
# Extra time series, same components and phases as eta and u of Q4a
# (direct sum in blocks of time steps, as in Q4a)
# ----------------------------------------------------------------------------
ub, ab, a_s, dudz = (np.zeros_like(t) for _ in range(4))
chunk = 4000
for s in range(0, t.size, chunk):
    th = 2 * np.pi * f[:, None] * t[None, s:s + chunk] + phi[:, None]
    c, sn = np.cos(th), np.sin(th)
    ub[s:s + chunk] = (a * H_ub) @ c
    ab[s:s + chunk] = -(a * H_ab) @ sn        # d/dt cos(w t + phi) = -w sin(w t + phi)
    a_s[s:s + chunk] = -(a * H_as) @ sn
    dudz[s:s + chunk] = (a * H_dudz) @ c
dpdx_b = -rho * ab                            # horizontal pressure gradient at the bed [Pa/m]
us_ext = us + eta * dudz                      # linear extrapolation of u(0) to z = eta (crests only)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def spectral(H):
    """Spectrum of r = H eta and its moments: rms, Tm02, bandwidth."""
    S = H**2 * E
    m = {n: np.sum(f**n * S) * df for n in (0, 1, 2, 4)}
    nu = np.sqrt(m[0] * m[2] / m[1]**2 - 1)          # Longuet-Higgins (1975) bandwidth
    eps = np.sqrt(1 - m[2]**2 / (m[0] * m[4]))        # Cartwright & Longuet-Higgins (1956)
    return S, np.sqrt(m[0]), np.sqrt(m[0] / m[2]), nu, eps


def cycle_peaks(x):
    """Zero up-crossing analysis: positive maximum and negative minimum of every cycle."""
    up = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    mx = np.array([x[i0:i1 + 1].max() for i0, i1 in zip(up[:-1], up[1:])])
    mn = np.array([-x[i0:i1 + 1].min() for i0, i1 in zip(up[:-1], up[1:])])
    return mx, mn, up


def rayleigh_stats(sig, N):
    """Rayleigh peaks with rms sig*sqrt(2): mean of highest 1/3, 1/10, expected maximum."""
    # mean of the highest fraction q of Rayleigh peaks, P(x > X) = exp(-X^2/(2 sig^2))
    def top(q):
        X = sig * np.sqrt(2 * np.log(1 / q))
        xs = np.linspace(X, 12 * sig, 20000)
        pdf = xs / sig**2 * np.exp(-xs**2 / (2 * sig**2))
        return np.trapezoid(xs * pdf, xs) / q
    L = np.sqrt(2 * np.log(N))
    return top(1 / 3), top(1 / 10), sig * (L + 0.5772 / L)


def skew_kurt(x):
    x = x - x.mean()
    return np.mean(x**3) / x.std()**3, np.mean(x**4) / x.std()**4


# ----------------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------------
signals = {
    "u_bed": (ub, H_ub, "m/s"),
    "u_surface": (us, H_us, "m/s"),
    "dudt_bed": (ab, H_ab, "m/s2"),
    "dudt_surface": (a_s, H_as, "m/s2"),
}
lines = []
p = lines.append
p("--- Question 6: linear statistics near the bed (z = -h) and near the surface (z = 0) ---")
p(f"h = {h:g} m, Hm0 = {Hs:g} m, Tp = {Tp:g} s, D = {D/3600:g} h, N_comp = {f.size}, "
  f"seed = {SEED}, rho = {rho:g} kg/m3")
p(f"check: Hm0 of eta series = {4*eta.std():.3f} m (Q4a: 6.001 m)")
p("")
p(f"{'signal':14s} {'rms series':>11s} {'rms spectrum':>13s} {'Tm02':>7s} {'nu':>6s} {'eps':>6s} "
  f"{'skew':>6s} {'kurt':>6s} {'N cyc':>6s} {'mean pk':>8s} {'r1/3':>7s} {'R:r1/3':>7s} "
  f"{'r1/10':>7s} {'R:r1/10':>7s} {'max':>7s} {'R:E[max]':>8s} {'-min':>7s}")
res = {}
for name, (x, H, unit) in signals.items():
    S, rms_sp, Tm02, nu, eps = spectral(H)
    mx, mn, up = cycle_peaks(x)
    N = mx.size
    srt = np.sort(mx)[::-1]
    r13, r110 = srt[:N // 3].mean(), srt[:N // 10].mean()
    R13, R110, Rmax = rayleigh_stats(x.std(), N)
    sk, ku = skew_kurt(x)
    res[name] = dict(x=x, S=S, rms=x.std(), rms_sp=rms_sp, Tm02=Tm02, mx=mx, mn=mn, N=N,
                     r13=r13, r110=r110, R13=R13, R110=R110, Rmax=Rmax, unit=unit)
    p(f"{name:14s} {x.std():11.4f} {rms_sp:13.4f} {Tm02:7.2f} {nu:6.3f} {eps:6.3f} {sk:+6.3f} {ku:6.3f} "
      f"{N:6d} {mx.mean():8.3f} {r13:7.3f} {R13:7.3f} {r110:7.3f} {R110:7.3f} {mx.max():7.3f} "
      f"{Rmax:8.3f} {mn.max():7.3f}   [{unit}]")
p("(R: = Rayleigh value from the rms of the series; E[max] = rms (sqrt(2 ln N) + 0.577/sqrt(2 ln N)))")
p(f"Tm02 of eta (spectral) = {spectral(np.ones_like(f))[2]:.2f} s")

rb = res["dudt_bed"]
p("")
p("--- Horizontal pressure gradient at the bed, dp/dx = -rho du/dt ---")
p(f"rms          = {rho*rb['rms']:8.1f} Pa/m = {rho*rb['rms']/(rho*g):.4f} rho g")
p(f"significant  = {rho*rb['r13']:8.1f} Pa/m = {rho*rb['r13']/(rho*g):.4f} rho g")
p(f"1/10         = {rho*rb['r110']:8.1f} Pa/m = {rho*rb['r110']/(rho*g):.4f} rho g")
mxp = max(rb['mx'].max(), rb['mn'].max())
p(f"maximum |dp/dx| in 3 h = {rho*mxp:8.1f} Pa/m = {mxp/g:.4f} rho g "
  f"(Rayleigh expected max {rho*rb['Rmax']:.1f} Pa/m)")
i_pg = np.argmax(np.abs(dpdx_b))
p(f"occurs at t = {t[i_pg]:.1f} s; there eta = {eta[i_pg]:+.2f} m, u_bed = {ub[i_pg]:+.2f} m/s")

# regular-wave check (H = Hs, T = Tp) and Soulsby-type representative values
kp = wavenumber(np.array([2 * np.pi / Tp]))[0]
p("")
p("--- Checks ---")
p(f"regular wave H = Hs, T = Tp: u_bed = pi H/(T sinh kh) = {np.pi*Hs/(Tp*np.sinh(kp*h)):.3f} m/s, "
  f"u_surface = {np.pi*Hs/(Tp*np.tanh(kp*h)):.3f} m/s, "
  f"dp/dx_bed = rho w^2 a/sinh kh = {rho*(2*np.pi/Tp)**2*Hs/2/np.sinh(kp*h):.0f} Pa/m")
p(f"u_bed significant amplitude 2 u_rms = {2*res['u_bed']['rms']:.3f} m/s "
  f"(= Soulsby 1997 significant bed velocity U_s = 2 sigma_u)")
mxe, _, _ = cycle_peaks(us_ext)
p(f"u at the instantaneous surface by linear extrapolation u(0) + eta du/dz: "
  f"max {mxe.max():.2f} m/s, mean highest 1/3 {np.sort(mxe)[::-1][:mxe.size//3].mean():.2f} m/s "
  f"(Wheeler stretching gives u(z = eta) = u(z = 0): max {res['u_surface']['mx'].max():.2f} m/s)")
# rms velocity over depth
zz = np.linspace(-h, 0, 41)
urms_z = np.array([np.sqrt(np.sum((omega * np.cosh(k * (z + h)) / np.sinh(k * h))**2 * E * df)) for z in zz])
p(f"u_rms(z): bed {urms_z[0]:.3f}, mid-depth {urms_z[20]:.3f}, surface {urms_z[-1]:.3f} m/s")

txt = "\n".join(lines)
print(txt)
with open(os.path.join(out_dir, "q6_output.txt"), "w") as fh:
    fh.write(txt + "\n")

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, name), dpi=200)
    print("saved", name)


def welch_psd(x, dt, nseg):
    win = np.hanning(nseg)
    segs = [x[i:i + nseg] for i in range(0, x.size - nseg + 1, nseg // 2)]
    P = np.mean([np.abs(np.fft.rfft((s_ - s_.mean()) * win))**2 for s_ in segs], axis=0)
    return np.fft.rfftfreq(nseg, dt), 2 * dt * P / np.sum(win**2)


# 1) time series around the largest pressure gradient at the bed
w = (t >= t[i_pg] - 60) & (t <= t[i_pg] + 60)
fig, ax = plt.subplots(3, 1, figsize=(8, 6.5), sharex=True)
ax[0].plot(t[w], eta[w], "C0"); ax[0].set_ylabel("η [m]")
ax[1].plot(t[w], us[w], "C1", label="u near the surface (z = 0)")
ax[1].plot(t[w], ub[w], "C2", label="u near the bed (z = −h)")
ax[1].set_ylabel("u [m/s]"); ax[1].legend(fontsize=8, loc="upper right")
ax[2].plot(t[w], dpdx_b[w] / 1000, "C3"); ax[2].set_ylabel("∂p/∂x at bed [kPa/m]")
ax[2].set_xlabel("t [s]")
for a_ in ax:
    a_.grid(alpha=0.3); a_.axhline(0, color="0.5", lw=0.5)
    a_.axvline(t[i_pg], color="0.5", ls=":", lw=1)
ax[0].set_title("Q6: 120 s around the largest pressure gradient at the bed (linear, 3 h record)")
save(fig, "fig_q6_timeseries.png")

# 2) spectra and rms profile
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for name, col, lab in (("u_surface", "C1", "z = 0"), ("u_bed", "C2", "z = −h")):
    fw, P = welch_psd(res[name]["x"], dt, 4096)
    ax[0].plot(fw, P, col, lw=1, label=f"{lab} (time series)")
    ax[0].plot(f, res[name]["S"], "k--", lw=0.8)
ax[0].plot([], [], "k--", lw=0.8, label="$|H_u|^2 E$")
ax[0].set_xlim(0, 0.55); ax[0].set_xlabel("f [Hz]"); ax[0].set_ylabel("$S_u$ [(m/s)²/Hz]")
ax[0].set_title("Velocity spectra"); ax[0].legend(fontsize=8)
for name, col, lab in (("dudt_surface", "C1", "z = 0"), ("dudt_bed", "C2", "z = −h")):
    fw, P = welch_psd(res[name]["x"], dt, 4096)
    ax[1].plot(fw, P, col, lw=1, label=f"{lab} (time series)")
    ax[1].plot(f, res[name]["S"], "k--", lw=0.8)
ax[1].set_xlim(0, 0.55); ax[1].set_xlabel("f [Hz]"); ax[1].set_ylabel("$S_{\\partial u/\\partial t}$ [(m/s²)²/Hz]")
ax[1].set_title("Acceleration spectra"); ax[1].legend(fontsize=8)
ax[2].plot(urms_z, zz, "k")
ax[2].plot(res["u_bed"]["rms"], -h, "o", color="C2"); ax[2].plot(res["u_surface"]["rms"], 0, "o", color="C1")
ax[2].set_xlabel("$u_{rms}$ [m/s]"); ax[2].set_ylabel("z [m]"); ax[2].set_xlim(0, None)
ax[2].set_title("rms velocity over the depth")
for a_ in ax:
    a_.grid(alpha=0.3)
save(fig, "fig_q6_spectra.png")

# 3) exceedance of the cycle maxima vs Rayleigh, and pdf of the instantaneous values
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
panels = (("u_bed", "C2", "u near the bed", "û [m/s]"),
          ("u_surface", "C1", "u near the surface", "û [m/s]"),
          ("dudt_bed", "C3", "∂p/∂x at the bed", "|∂p/∂x| [kPa/m]"))
for j, (name, col, title, xl) in enumerate(panels):
    r = res[name]
    scale = rho / 1000 if name == "dudt_bed" else 1.0
    for v, mk, lab in ((r["mx"], ".", "positive peaks"), (r["mn"], "x", "negative peaks")):
        vs = np.sort(v)[::-1] * scale
        ax[j].semilogy(vs, np.arange(1, vs.size + 1) / vs.size, mk, color=col, ms=3,
                       alpha=0.8 if mk == "." else 0.4, label=lab)
    xr = np.linspace(0, 1.15 * max(r["mx"].max(), r["mn"].max()) * scale, 200)
    ax[j].semilogy(xr, np.exp(-xr**2 / (2 * (r["rms"] * scale)**2)), "k--", lw=1,
                   label="Rayleigh, $\\exp(-x^2/2\\sigma^2)$")
    ax[j].axvline(r["r13"] * scale, color="0.4", ls=":", lw=1)
    ax[j].text(r["r13"] * scale, 0.5, " 1/3", fontsize=8, color="0.3")
    ax[j].set_ylim(1 / 3000, 1.1); ax[j].set_xlim(0, xr[-1])
    ax[j].set_xlabel(xl); ax[j].set_title(title); ax[j].grid(alpha=0.3, which="both")
ax[0].set_ylabel("exceedance probability"); ax[0].legend(fontsize=8, loc="lower left")
save(fig, "fig_q6_exceedance.png")

plt.show()

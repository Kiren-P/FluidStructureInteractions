
import glob, os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

h = 25.0        # water depth [m]
U = 1.5         # mean current velocity [m/s]
KAPPA = 0.41    # von Karman constant
RHO = 1000.0    # water density [kg/m3]
ZMAX = 0.3      # upper limit of the log layer, as fraction of h (points with z/h <= ZMAX are fitted)

# ---- read data: mean u at every height ----
folder = os.path.dirname(os.path.abspath(__file__))
files = sorted(glob.glob(os.path.join(folder, "Z_*.dat")))  
if not files:
    files = sorted(glob.glob(os.path.join(folder, "Velocity files for FSI", "Z_*.dat")))
if not files:
    raise FileNotFoundError(f"No Z_*.dat files in {folder}")

z, u_mean = [], []
for f in files:
    z_over_h = float(re.search(r"Z_([\d.]+)\.dat", os.path.basename(f)).group(1)) / 100   
    u_t = pd.read_csv(f, sep=r"\s+", header=None, usecols=[1]).to_numpy().ravel()        
    z.append(z_over_h * h)               # z [m]
    u_mean.append((u_t * U).mean())      # [m/s]
z, u_mean = np.array(z), np.array(u_mean)


fit = z <= ZMAX * h
a, b = np.polyfit(np.log(z[fit]), u_mean[fit], 1)
u_star = KAPPA * a
z0 = np.exp(-b / a)
tau_b = RHO * u_star**2
u_pred = a * np.log(z[fit]) + b
r2 = 1 - np.sum((u_mean[fit] - u_pred) ** 2) / np.sum((u_mean[fit] - u_mean[fit].mean()) ** 2)


print(f"{'z/h':>6} {'z [m]':>7} {'u [m/s]':>8}  log fit")
for zi, ui, fi in zip(z, u_mean, fit):
    print(f"{zi/h:6.3f} {zi:7.2f} {ui:8.4f}  {'yes' if fi else '':>7}")

print(f"\nLog-law fit: u = {a:.4f} ln(z) + {b:.4f}   (R^2 = {r2:.4f}, {fit.sum()} points, z/h <= {ZMAX})")
print(f"   u*    = {u_star:.4f} m/s")
print(f"   z0    = {z0:.4f} m,  ks ~ 30 z0 = {30 * z0:.2f} m")
print(f"   tau_b = rho*u*^2 = {tau_b:.2f} Pa")

text_file = open("Output_log_fit.txt", "w")

text_file.write(
    f"""
Log-law fit: u = {a:.4f} ln(z) + {b:.4f}   (R^2 = {r2:.4f}, {fit.sum()} points, z/h <= {ZMAX})
   u*    = {u_star:.4f} m/s
   z0    = {z0:.4f} m,  ks ~ 30 z0 = {30 * z0:.2f} m
   tau_b = rho*u*^2 = {tau_b:.2f} Pa        

    """
)
text_file.close()

fig, ax = plt.subplots(figsize=(6, 4.5))
lnz = np.linspace(np.log(z.min()), np.log(z.max()), 50)
ax.plot(np.log(z[~fit]), u_mean[~fit], "o", mfc="none", label="measured (not fitted)")
ax.plot(np.log(z[fit]), u_mean[fit], "o", label="measured (fitted)")
ax.plot(lnz, a * lnz + b, "-", label=f"fit: u* = {u_star:.3f} m/s, z0 = {z0:.3f} m")
ax.set(xlabel="ln(z)  [z in m]", ylabel="u [m/s]", title="Log-law fit")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout()
plt.savefig(os.path.join(folder, "log_fit_ustar.png"), dpi=200)
plt.show()

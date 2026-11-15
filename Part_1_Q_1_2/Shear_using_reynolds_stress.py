
import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt



folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Velocity files for FSI")

# Find all .dat files
dat_files = sorted(glob.glob(os.path.join(folder_path, "Z_*.dat")))
if not dat_files:
    raise FileNotFoundError(f"No Z_*.dat files found in {folder_path}")

tau_R = []
heights = []

U = 1.5  # mean current velocity [m/s]


# Read and print each file
for file in dat_files:
    # Read the data (whitespace-separated, no header)
    data = pd.read_csv(file, sep=r"\s+", header=None)

    # Assign column names
    data.columns = ["Time", "Streamwise velocity (u)", "Upward velocity (v)"]

    u_mean = data["Streamwise velocity (u)"].mean()
    v_mean = data["Upward velocity (v)"].mean() 

    u_fluc = data["Streamwise velocity (u)"]*U - u_mean*U
    v_fluc = data["Upward velocity (v)"]*U - v_mean*U

    # Calculate Reynolds stress
    reynolds_stress = -1000 * np.mean(u_fluc * v_fluc)

    z = float(os.path.basename(file)[2:-4]) / 100 * 25

    print(reynolds_stress, z)

    tau_R.append(reynolds_stress)
    heights.append(z)

heights, stresses = zip(*sorted(zip(heights, tau_R)))

c1a, c0a = np.polyfit(heights[1:], stresses[1:], 1)   # Fit a linear model to the data, excluding the first point (index 0)
c1b, c0b = np.polyfit(heights[2:], stresses[2:], 1)   # Fit a linear model to the data, excluding the first two points (index 0 and 1)
c1c, c0c = np.polyfit(heights[3:], stresses[3:], 1)   # Fit a linear model to the data, excluding the first three points (index 0, 1, and 2)
c1d, c0d = np.polyfit(heights[4:], stresses[4:], 1)   # Fit a linear model to the data, excluding the first four points (index 0, 1, 2, and 3)

plt.plot(heights, tau_R, marker='o')
plt.plot(heights, c1a * np.array(heights) + c0a, label='Fit (excluding 1st point)', linestyle='--')
plt.plot(heights, c1b * np.array(heights) + c0b, label='Fit (excluding 2nd point)', linestyle='--')
plt.plot(heights, c1c * np.array(heights) + c0c, label='Fit (excluding 3rd point)', linestyle='--')
plt.plot(heights, c1d * np.array(heights) + c0d, label='Fit (excluding 4th point)', linestyle='--')
plt.xlabel('Height (m)')
plt.ylabel('Reynolds Stress (Pa)')
plt.title('Reynolds Stress vs Height')
plt.legend()
plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reynolds_stress_vs_height.png'), dpi=150)
plt.show()


print (f"Reynolds stress fit (excluding first point): tau_R = {c1a:.4f} z + {c0a:.3f}")
print (f"Reynolds stress fit (excluding first two points): tau_R = {c1b:.4f} z + {c0b:.3f}")
print (f"Reynolds stress fit (excluding first three points): tau_R = {c1c:.4f} z + {c0c:.3f}")
print (f"Reynolds stress fit (excluding first four points): tau_R = {c1d:.4f} z + {c0d:.3f}")

text_file = open("Output_reynolds_fit.txt", "w")

text_file.write(
    f"""
    Reynolds stress fit (excluding first point): tau_R = {c1a:.4f} z + {c0a:.3f}
    Reynolds stress fit (excluding first two points): tau_R = {c1b:.4f} z + {c0b:.3f}
    Reynolds stress fit (excluding first three points): tau_R = {c1c:.4f} z + {c0c:.3f}
    Reynolds stress fit (excluding first four points): tau_R = {c1d:.4f} z + {c0d:.3f}

    """
)
text_file.close()
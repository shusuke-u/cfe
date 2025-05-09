import pandas as pd
import matplotlib.pyplot as plt

# File path to the CSV file
file_path = '../../output/result/conf.csv'

# Read the data from the CSV file
data = pd.read_csv(file_path)

# Font size for the plots
fs = 24

# Extract the relevant columns from the data
rep = data['repulsion_factor']
att = data['attraction_factor']
poisson_y = data['poisson_y']
poisson_z = data['poisson_z']
young = data['young']

# First scatter plot: Poisson's ratio (y and z) vs repulsion/attraction
plt.figure(figsize=(8, 6))
plt.scatter(rep / att, poisson_y, color='r', s=300)
plt.scatter(rep / att, poisson_z, edgecolor='b', facecolors='none', s=500, marker='s')

plt.xscale('log')
# plt.xlabel('repulsion / attraction', size=fs)
# plt.ylabel("Poisson's ratio", size=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend(fontsize=fs)
plt.grid(True, which="both", ls="--")
plt.savefig("../res/poisson.pdf")

# Second scatter plot: Young's modulus vs repulsion/attraction
plt.figure(figsize=(8, 6))
plt.scatter(rep / att, young, color='r', s=300)
# plt.xlabel('repulsion / attraction', size=fs)
# plt.ylabel("Young's modulus [Pa]", size=fs)
plt.ylim()
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xscale('log')
plt.grid(True, which="both", ls="--")
plt.savefig("../res/young.pdf")

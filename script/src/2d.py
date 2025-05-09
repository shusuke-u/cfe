import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Define the function to read your data files
def read_data(file_path):
    try:
        data = np.genfromtxt(file_path, delimiter='   ')  # Adjust delimiter if necessary
        if data.ndim == 1:  # Ensure data is at least 2D
            data = data.reshape(-1, 2)
        return data
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return np.array([[0, 0]])

# List all data files and sort by the number in the file name
data_folder = '../../output/distr'  # Replace with your data folder path
data_files = sorted([os.path.join(data_folder, f) for f in os.listdir(data_folder) if f.endswith('.csv')], key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[-1]))

# Check if data files are listed correctly
print(data_files)

# Create the plot and initialization function
fig, ax = plt.subplots()
ax.set_aspect('equal')  # Ensure the aspect ratio is equal
# scat = ax.scatter([], [], c='r', marker='o', s=0.1)
scat = ax.scatter([], [], c='r', marker='o', s=0.2)
text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def init():
    ax.set_xlim(-0.75, 0.75)  # Set according to your data range
    ax.set_ylim(-0.75, 0.75)  # Set according to your data range
    scat.set_offsets(np.empty((0, 2)))
    text.set_text('')
    return scat, text

# Define the animation function
def animate(i):
    data = read_data(data_files[i])
    x = data[:, 2]  # Adjust the column index based on your data
    y = data[:, 3]  # Adjust the column index based on your data
    scat.set_offsets(np.c_[x, y])
    step_number = int(os.path.splitext(os.path.basename(data_files[i]))[0].split('_')[-1])
    text.set_text(f'step={step_number}')
    return scat, text

# Create the animation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(data_files), interval=100, blit=True)

writer = animation.PillowWriter(fps=10)
ani.save('../output/2d.gif', writer=writer, dpi=200)

# Display the animation in Jupyter notebook
from IPython.display import HTML

HTML(ani.to_jshtml())

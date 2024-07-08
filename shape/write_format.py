import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

input_path = 'shape/out/'
output_path = 'raytrace_efficient/vertices/'


radius = 0.005

nx = 150

nums = [25,30,35,80,90,100,110,120,130,140,150,175,200,225,250,275,300,325]

for num in nums:
    file_name = f"{num}-smoothed.csv"
    df = pd.read_csv(input_path+file_name)
    xspace = df['x'].to_numpy()
    yspace = df['y'].to_numpy()

    #rotate and move


    xspace = xspace * -1
    yspace = yspace * -1

    xspace -= np.min(xspace)
    yspace -= np.min(yspace)

    scale = radius/np.max(np.abs(xspace))

    xspace *= scale
    yspace *= scale


    plt.plot(xspace,yspace)
    plt.show()

    yspace_zero = np.zeros_like(yspace)  # y is always zero
    data = np.column_stack((xspace, yspace_zero, yspace))

    # Write the data to the file
    output_name = f"{num}_vertices.txt"
    with open(output_path + output_name, "w") as f:
        f.write(f"{nx} {nx}\n")  # Write nx repeated twice with a space in between
        np.savetxt(f, data, fmt='%.20f')  # Write the data with 20 decimal places
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("./surface_assignments.txt", delim_whitespace=True, header=None)
x = df.to_numpy()[0]

plt.hist(x, bins=max(x)+1, rwidth=0.7)
plt.show()
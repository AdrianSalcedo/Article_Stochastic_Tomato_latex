import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
import matplotlib as mpl

'''
with open('Tomato_data_det.csv', newline='') as File:
    reader = csv.reader(File)
    #for row in reader:
    #   print(row)
print(reader)
'''

df1 = pd.read_csv('Tomato_data_det_R0-1.csv')
df2 = pd.read_csv('Tomato_data_stc_R0-1.csv')

mpl.style.use('ggplot')

ax0 = plt.subplot2grid((2, 3), (0, 0), rowspan=1)
ax0.plot(df1[["t"]],df1[["x"]],color = 'blue', label="deterministic Solution")
ax0.plot(df2[["t"]],df2[["x"]], color ='red', label="stochastic Solution")
ax0.set_ylabel('$Susceptibles plants$')
ax0.grid(True)

ax1 = plt.subplot2grid((2, 3), (0, 1), rowspan=1)
ax1.plot(df1[["t"]],df1[["y"]],color = 'blue', label="deterministic Solution")
ax1.plot(df2[["t"]],df2[["y"]], color ='red', label="stochastic Solution")
ax1.set_xlabel('$t$')
ax1.set_ylabel('$Latents plants$')
ax1.grid(True)

ax2 = plt.subplot2grid((2, 3), (0, 2), rowspan=1)
ax2.plot(df1[["t"]],df1[["z"]],color = 'blue', label="deterministic Solution")
ax2.plot(df2[["t"]],df2[["z"]], color ='red', label="stochastic Solution")
ax2.set_xlabel('$t$')
ax2.set_ylabel('$Infected plants$')
ax2.grid(True)

ax3 = plt.subplot2grid((2, 3), (1, 0), rowspan=1)
ax3.plot(df1[["t"]],df1[["v"]],color = 'blue', label="deterministic Solution")
ax3.plot(df2[["t"]],df2[["v"]], color ='red', label="stochastic Solution")
ax3.set_xlabel('$t$')
ax3.set_ylabel('$Susceptibles vector$')
ax3.grid(True)

ax4 = plt.subplot2grid((2, 3), (1, 1), rowspan=1)
ax4.plot(df1[["t"]],df1[["w"]],color = 'blue', label="deterministic Solution")
ax4.plot(df2[["t"]],df2[["w"]], color ='red', label="stochastic Solution")
ax4.set_xlabel('$t$')
ax4.set_ylabel('$Infected vector$')
ax4.grid(True)


plt.tight_layout()

plt.show()


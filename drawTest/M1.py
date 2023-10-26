from matplotlib import pyplot as plt
import numpy as np
import math
x = np.arange(0,math.pi*2,0.05)
y = np.sin(x)
fig = plt.figure()
ax = fig.add_axes([0,0,0.5,0.5])
ax.plot(x, y)
plt.show()
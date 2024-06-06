#Hermit Purple by Ramirez Garrido
#wave function ploting
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cbook

f = open("/data.txt", 'rt')
data = f.readlines()
x = []
yh = []
p = []

for line in data:
    values = line.split(",")
    x.append(float(values[0]))
    yh.append(float(values[2]))
    p.append(float(values[3]))

plt.style.use('bmh')
plt.plot(x, yh,'#177ba3', label='$\psi_n(x)$')
plt.plot(x,p, '#7f059e', label='$|\psi_n(x)|^2$')
plt.xlabel('$x$')
fig = plt.gcf()
fig.set_facecolor('white')
plt.legend(loc='upper left', labelcolor = 'linecolor')
plt.title(label = 'Wave function n = 5', color = 'gray', fontstyle='italic', fontsize = 11.5)
plt.tight_layout()
!plt.savefig('schodingerCat3.7.png',  transparent=True, dpi=600)

plt.show()

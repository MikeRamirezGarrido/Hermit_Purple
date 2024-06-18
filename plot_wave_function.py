#Hermit Purple by Ramirez Garrido
#wave function ploting
import matplotlib.pyplot as plt

f = open("/content/drive/My Drive/Data/uni_5_data.txt", 'rt')
data = f.readlines()
x = []
y = []
p = []

for line in data:
    values = line.split(",")
    x.append(float(values[0]))
    y.append(float(values[1]))
    p.append(float(values[2]))

plt.style.use('bmh')
plt.plot(x, y,'#7f059e', label='$\psi_n(x)$')
plt.plot(x,p, '#fa7b05', label='$|\psi_n(x)|^2$')
plt.xlabel('$x$  $[m]$')
fig = plt.gcf()
fig.set_facecolor('white')
plt.legend(loc='upper left', labelcolor = 'linecolor')
plt.title(label = 'Wave function n = 5', color = 'gray', fontstyle='italic', fontsize = 11.5)
plt.tight_layout()
plt.savefig('unidim_5_harm.png' ,  transparent=True, dpi=600)
plt.show()

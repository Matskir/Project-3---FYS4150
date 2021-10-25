import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 15})


def read_file(filename):
    infile = open(filename)
    x = []; y = []

    for line in infile:
        elements = line.split()
        x.append(float(elements[0]))
        y.append(float(elements[1]))
    x = np.array(x)
    y = np.array(y)
    infile.close()

    return x,y

x,y = read_file("resonance_test_f04.txt")
x1,y1 = read_file("resonance_test_f07.txt")
x2,y2 = read_file("resonance_test_f01.txt")

plt.figure(0)
plt.plot(y2,x2/100,label='f = 0.1')
plt.plot(y,x/100,label='f = 0.4')
plt.plot(y1,x1/100,label='f = 0.7')
plt.xlabel('$\omega_V$ [$MH_z$]')
plt.ylabel('Fraction of particles trapped at $t = 500\,\,\,\mu m$')
plt.legend()
plt.tight_layout()
plt.savefig('resonance.pdf')

plt.figure(1)
x3,y3 = read_file("resonance_test_f01_zoom.txt")
x4,y4 = read_file("resonance_test_f07_zoom.txt")
x5,y5 = read_file("resonance_test_f04_zoom.txt")
plt.plot(y3,x3/100,label='f = 0.1')
plt.plot(y5,x5/100,label='f = 0.4')
plt.plot(y4,x4/100,label='f = 0.7')
plt.xlabel('$\omega_V$ [$MH_z$]')
plt.ylabel('Fraction of particles trapped at $t = 500\,\,\,\mu m$')
plt.legend()
plt.savefig('zoom_nocoul.pdf')

plt.figure(2)
x6,y6 = read_file("resonance_test_f04_zoom_coul.txt")
plt.plot(y6,x6/100,label='With coulomb forces')
plt.plot(y5,x5/100,label='Without coulomb forces')
plt.xlabel('$\omega_V$ [$MH_z$]')
plt.ylabel('Fraction of particles trapped at $t = 500\,\,\,\mu m$')
plt.legend()
plt.tight_layout()
plt.savefig('zoom_coul.pdf')
plt.show()

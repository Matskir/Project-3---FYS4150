#Imports
import numpy as np
import matplotlib.pyplot as plt

#Function that reads the files
def read_file(filename):
    infile = open(filename)
    x = []; y = []; z = []
    vx = []; vy = []; vz = []

    for line in infile:
        elements = line.split()
        x.append(float(elements[0]))
        y.append(float(elements[1]))
        z.append(float(elements[2]))
        vx.append(float(elements[3]))
        vy.append(float(elements[4]))
        vz.append(float(elements[5]))
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    vx = np.array(vx)
    vy = np.array(vy)
    vz = np.array(vz)
    infile.close()

    return x,y,z,vx,vy,vz

#Define constants and initial values
x_init = 1
z_init = 1
v0 = 1
q = 1
m = 40.08
B0 = 96.5
V0 = 9.65e8
d = 1e4
w0 = q*B0/m
wz = np.sqrt(2*q*V0/(m*d**2))

#Function that computes omega+-
def w(w0,wz,sign):
    if sign=="+":
        return (w0+np.sqrt(w0**2 - 2*wz**2))/2
    elif sign=="-":
        return (w0-np.sqrt(w0**2 - 2*wz**2))/2
    else:
        "Please use either '+' or '-' as sign!"

#Function that computes the analytical solution positions
def compute_analytical(w0,wz,t):
    A_plus = (v0+w(w0,wz,"-")*x_init)/(w(w0,wz,"-")-w(w0,wz,"+"))
    A_min = -(v0+w(w0,wz,"+")*x_init)/(w(w0,wz,"-")-w(w0,wz,"+"))
    i = complex(0,1)
    f = A_plus*np.exp(-i*w(w0,wz,"+")*t) + A_min*np.exp(-i*w(w0,wz,"-")*t)
    x = f.real; y = f.imag
    z = z_init*np.cos(np.sqrt(2*q*V0/(m*d**2))*t)
    return x,y,z

#Function that computes the relative error
def relative_error(x_num, y_num, z_num, x_an, y_an, z_an):
    r_num = np.array([x_num, y_num, z_num])
    r_an = np.array([x_an, y_an, z_an])

    rel_err = np.zeros(len(x_an))
    for i in range(len(rel_err)):
        rel_err[i] = abs(np.linalg.norm(r_num[:,i]-r_an[:,i])/np.linalg.norm(r_num[:,i]))
    return rel_err


#Reading a bunch of data files created with Project3_FYS4150_2.cpp and computes the analytical solution and relative error
x0, y0, z0, vx0, vy0, vz0 = read_file("one_particle_dt0001_RKt.txt")
t0 = np.linspace(0,100,len(x0))
x0_an, y0_an, z0_an = compute_analytical(w0,wz,t0)
rel_err0 = relative_error(x0,y0,z0,x0_an,y0_an,z0_an)

x1, y1, z1, vx1, vy1, vz1 = read_file("one_particle_dt001_RK4.txt")
t1 = np.linspace(0,100,len(x1))
x1_an, y1_an, z1_an = compute_analytical(w0,wz,t1)
rel_err1 = relative_error(x1,y1,z1,x1_an,y1_an,z1_an)

x2, y2, z2, vx2, vy2, vz2 = read_file("one_particle_dt01_RK4.txt")
t2 = np.linspace(0,100,len(x2))
x2_an, y2_an, z2_an = compute_analytical(w0,wz,t2)
rel_err2 = relative_error(x2,y2,z2,x2_an,y2_an,z2_an)

x3, y3, z3, vx3, vy3, vz3 = read_file("one_particle_dt1_RK4.txt")
t3 = np.linspace(0,100,len(x3))
x3_an, y3_an, z3_an = compute_analytical(w0,wz,t3)
rel_err3 = relative_error(x3,y3,z3,x3_an,y3_an,z3_an)

x4, y4, z4, vx4, vy4, vz4 = read_file("one_particle_dt00001_RK4.txt")
t4 = np.linspace(0,100,len(x4))
x4_an, y4_an, z4_an = compute_analytical(w0,wz,t4)
rel_err4 = relative_error(x4,y4,z4,x4_an,y4_an,z4_an)

#Plots the logarithm of the relative errors
plt.figure(0)
plt.plot(t4[1:],np.log(rel_err4[1:]),label='dt=0.00001')
plt.plot(t0[1:],np.log(rel_err0[1:]),label='dt=0.0001')
plt.plot(t1[1:],np.log(rel_err1[1:]),label='dt=0.001')
plt.plot(t2[1:],np.log(rel_err2[1:]),label='dt=0.01')
plt.plot(t3[1:],np.log(rel_err3[1:]),label='dt=0.1')
plt.xlabel('t')
plt.ylabel('$\log(\epsilon)$')
plt.legend()
#plt.savefig('relative_errors.pdf')



#Compute the error convergence of the used method
import numpy.linalg
h = np.array([0.1,0.01,0.001,0.0001,0.00001])

#Make arrays of the radii
r1 = np.sqrt(x3**2 + y3**2 + z3**2)
r1_an = np.sqrt(x3_an**2 + y3_an**2 + z3_an**2)

r2 = np.sqrt(x2**2 + y2**2 + z2**2)
r2_an = np.sqrt(x2_an**2 + y2_an**2 + z2_an**2)

r3 = np.sqrt(x1**2 + y1**2 + z1**2)
r3_an = np.sqrt(x1_an**2 + y1_an**2 + z1_an**2)

r4 = np.sqrt(x0**2 + y0**2 + z0**2)
r4_an = np.sqrt(x0_an**2 + y0_an**2 + z0_an**2)

r5 = np.sqrt(x4**2 + y4**2 + z4**2)
r5_an = np.sqrt(x4_an**2 + y4_an**2 + z4_an**2)

r = np.array([r1,r2,r3,r4,r5],dtype=object)
r_an = np.array([r1_an,r2_an,r3_an,r4_an,r5_an],dtype=object)

#Compute the maximum absolute error
Del_max = np.zeros(len(h))
for i in range(len(h)):
    Del_max[i] = np.max(r_an[i]-r[i])

#Compute the error convergence rate
s = 0
for i in range(1,5):
    s += np.log(Del_max[i]/Del_max[i-1])/np.log(h[i]/h[i-1])
rerr = 1/4*s
print(rerr)

plt.show()

#Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d

#Set font size
matplotlib.rcParams.update({'font.size': 15})

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


#Read data files
x,y,z,vx,vy,vz = read_file("one_particle_dt001_RK4.txt")
x_p1, y_p1, z_p1, vx_p1, vy_p1, vz_p1 = read_file("two_particles_part1_dt001_RK4.txt")
x_p2, y_p2, z_p2, vx_p2, vy_p2, vz_p2 = read_file("two_particles_part2_dt001_RK4.txt")

x_p1_nocol, y_p1_nocol, z_p1_nocol, vx_p1_nocol, vy_p1_nocol, vz_p1_nocol = read_file("two_particles_part1_dt001_RK4_NoCoul.txt")
x_p2_nocol, y_p2_nocol, z_p2_nocol, vx_p2_nocol, vy_p2_nocol, vz_p2_nocol = read_file("two_particles_part2_dt001_RK4_NoCoul.txt")

#Compute analytical solution
t0 = np.linspace(0,100,len(x))
x_an, y_an, z_an = compute_analytical(w0,wz,t0)

#Plots all needed plots for the Result section
fig, ax = plt.subplots(figsize=(8,8))
ax.plot(x,y,label='Numerical solution')
ax.plot(x_an,y_an,'r--',label='Analytical solution')
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('y [$\mu$m]')
ax.legend()
fig.savefig('one_particle_analytical.pdf')

fig1, ax1 = plt.subplots(2,1,figsize=(8,8))
ax1[1].plot(x_p1,vx_p1)
ax1[0].plot(x_p1_nocol, vx_p1_nocol)
ax1[0].set_xlabel('x [$\mu$m]')
ax1[0].set_ylabel('$v_x$ [$\mu$m/$\mu$s]')
ax1[1].set_xlabel('x [$\mu$m]')
ax1[1].set_ylabel('$v_x$ [$\mu$m/$\mu$s]')
ax1[0].set_title('Phase space plot (x,$v_x$) without coulomb forces')
ax1[1].set_title('Phase space plot (x,$v_x$) with coulomb forces')
plt.tight_layout()
fig1.savefig('x_phase.pdf')

fig2, ax2 = plt.subplots(2,1,figsize=(8,8))
ax2[1].plot(y_p1,vy_p1)
ax2[0].plot(y_p1_nocol,vy_p1_nocol)
ax2[0].set_xlabel('y [$\mu$m]')
ax2[0].set_ylabel('$v_y$ [$\mu$m/$\mu$s]')
ax2[1].set_xlabel('y [$\mu$m]')
ax2[1].set_ylabel('$v_y$ [$\mu$m/$\mu$s]')
ax2[0].set_title('Phase space plot (y,$v_y$) without coulomb forces')
ax2[1].set_title('Phase space plot (y,$v_y$) with coulomb forces')
plt.tight_layout()
fig2.savefig('y_phase.pdf')

fig3, ax3 = plt.subplots(2,1,figsize=(8,8))
ax3[1].plot(z_p1,vz_p1)
ax3[0].plot(z_p1_nocol,vz_p1_nocol)
ax3[0].set_xlabel('z [$\mu$m]')
ax3[0].set_ylabel('$v_z$ [$\mu$m/$\mu$s]')
ax3[1].set_xlabel('z [$\mu$m]')
ax3[1].set_ylabel('$v_z$ [$\mu$m/$\mu$s]')
ax3[0].set_title('Phase space plot (z,$v_z$) without coulomb forces')
ax3[1].set_title('Phase space plot (z,$v_z$) with coulomb forces')
plt.tight_layout()

fig3.savefig('z_phase.pdf')

fig4, ax4 = plt.subplots(2,1,figsize=(8,16))
ax4[0].plot(x_p1, y_p1, label='Particle 1')
ax4[0].plot(x_p2, y_p2, label='Particle 2')
ax4[1].plot(x_p1_nocol, y_p1_nocol, label='Particle 1')
ax4[1].plot(x_p2_nocol, y_p2_nocol, label='Particle 2')
ax4[0].set_xlabel('x [$\mu$m]')
ax4[0].set_ylabel('y [$\mu$m]')
ax4[1].set_xlabel('x [$\mu$m]')
ax4[1].set_ylabel('y [$\mu$m]')
ax4[0].set_aspect('equal')
ax4[1].set_aspect('equal')
ax4[0].legend()
ax4[1].legend()
fig4.savefig('two_particles_motion.pdf')

fig5, ax5 = plt.subplots()
ax5.plot(t0,z,label='Numerical solution')
ax5.plot(t0,z_an,'r--',label='Analytical solution')
ax5.set_xlabel('t [$\mu$s]')
ax5.set_ylabel('$z$ [$\mu$m]')
ax5.legend()
plt.tight_layout()
fig5.savefig('z_fig.pdf')

fig6, ax6 = plt.subplots(figsize=(8,8))
ax6.plot(x_p1, y_p1, label='Particle 1')
ax6.plot(x_p2, y_p2, label='Particle 2')
ax6.set_xlabel('x [$\mu$m]')
ax6.set_ylabel('y [$\mu$m]')
ax6.set_aspect('equal')
ax6.legend()
fig6.savefig('coul_two.pdf')

fig7, ax7 = plt.subplots(figsize=(8,8))
ax7.plot(x_p1_nocol, y_p1_nocol, label='Particle 1')
ax7.plot(x_p2_nocol, y_p2_nocol, label='Particle 2')
ax7.set_xlabel('x [$\mu$m]')
ax7.set_ylabel('y [$\mu$m]')
ax7.set_aspect('equal')
ax7.legend()
fig7.savefig('no_coul_two.pdf')

fig8 = plt.figure(9)
ax8 = plt.axes(projection='3d')
ax8.plot3D(x_p1,y_p1,z_p1,label='Particle 1')
ax8.plot3D(x_p2,y_p2,z_p2,label='Particle 2')
ax8.set_xlabel('x [$\mu$m]')
ax8.set_ylabel('y [$\mu$m]')
ax8.set_zlabel('z [$\mu$m]')
ax8.legend()
plt.tight_layout()
fig8.savefig('3d_coul.pdf')

fig9 = plt.figure(10)
ax9 = plt.axes(projection='3d')
ax9.plot3D(x_p1_nocol,y_p1_nocol,z_p1_nocol,label='Particle 1')
ax9.plot3D(x_p2_nocol,y_p2_nocol,z_p2_nocol,label='Particle 2')
ax9.set_xlabel('x [$\mu$m]')
ax9.set_ylabel('y [$\mu$m]')
ax9.set_zlabel('z [$\mu$m]')
ax9.legend()
plt.tight_layout()
fig9.savefig('3d_nocoul.pdf')
plt.show()

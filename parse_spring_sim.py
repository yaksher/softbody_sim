import ProgressBar as PB
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import OdeSolution, DenseOutput
import numpy as np

class LinDense(DenseOutput):
    def __init__(self, t_old, t, y_old, y):
        super().__init__(t_old, t)
        self.y_old = np.array(y_old)
        self.y = np.array(y)
        self.dy = self.y - self.y_old
        self.dt = self.t - self.t_old

    def _call_impl(self, t):
        derv = (t - self.t_old)/self.dt
        return self.dy * derv + self.y_old


f = open("sim_data.txt", "r")
lines = f.readlines()
dims = int(lines[0])
n_objs = int(lines[1])
#obj_types = [int(x) for x in lines[2].split(",")]
obj_colors_line = [float(x) for x in lines[2].split(",")]
obj_colors = [tuple(obj_colors_line[i:i+3]) for i in range(0, 3 * n_objs, 3)]
try:
    sphere_radii = [float(x) for x in lines[3].split(",")]
except ValueError:
    sphere_radii = []
n_spheres = len(sphere_radii)
n_points = n_objs - n_spheres
data = []
for line in lines[4:]:
    data.append([float(x) for x in line.split(",")])

t_s = []
y_s = []
for line in data:
    t_s.append(line[-1])
    y_s.append(np.array(line[:-1]))

introps = []
for i in range(len(t_s) - 1):
    introps.append(LinDense(t_s[i], t_s[i + 1], y_s[i], y_s[i + 1]))

sol = OdeSolution(t_s, introps)

t0 = t_s[0]  # initial time
tf = t_s[-1] # final time
nframes = int(50*(tf - t0))
t = np.linspace(t0, tf, nframes)

sph_count = 15
sphere_render_angles = [(phi, theta) for phi in np.arange(0, sph_count)/sph_count * np.pi for theta in np.arange(0, 2 * sph_count)/sph_count * np.pi]
x_sphere = np.array([np.sin(a[0])*np.cos(a[1]) for a in sphere_render_angles])
y_sphere = np.array([np.sin(a[0])*np.sin(a[1]) for a in sphere_render_angles])
z_sphere = np.array([np.cos(a[0]) for a in sphere_render_angles])
def drawSphere(xCenter, yCenter, zCenter, r):
    x = r*x_sphere + xCenter
    y = r*y_sphere + yCenter
    z = r*z_sphere + zCenter
    return (list(x), list(y), list(z))

"""colors = []
for obj in obj_types:
    if obj == 0:
        colors.append((1, 0, 0))
    elif obj == 1:
        colors += len(sphere_render_angles) * [(0, 1, 1)]
    elif obj == 2:
        colors.append((0, 0, 1))
    else:
        colors += len(sphere_render_angles) * [(0, 1, 0)]"""
colors = obj_colors

frame0 = []
for d in range(dims):
    frame0.append(list(y_s[0])[d:dims * n_points:dims])
for d in range(dims):
    for s in range(n_spheres):
        frame0[d] += drawSphere(*y_s[0][(n_points + s)*dims:(n_points + s + 1)*dims], sphere_radii[s])[d]

if dims == 2:
    fig = plt.figure()
    ax = plt.axes(xlim=(-11.0, 11.0), ylim=(-11.0, 11.0))
    ax.set_aspect('equal', adjustable='box')
    graph, = ax.plot(*frame0, color=colors, lw=5.0)
    #graph.set_color(colors)
    ax.grid()
elif dims == 3:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x_lim = max([abs(x) for x in frame0[0]]) + 1
    y_lim = max([abs(y) for y in frame0[1]]) + 1
    z_lim = max([abs(z) for z in frame0[2]]) + 1
    lim = max(x_lim, y_lim, z_lim)
    ax.set_xlim3d(-lim,lim)
    ax.set_ylim3d(-lim,lim)
    ax.set_zlim3d(-lim/2,lim)
    graph = ax.scatter(*frame0, c = colors)
    #graph.set_color(colors)



def animate(i):
    pb2.count()
    row = sol(t[i])
    frame = []
    for d in range(dims):
        frame.append(list(row)[d:dims * n_points:dims])
    for d in range(dims):
        for s in range(n_spheres):
            frame[d] += drawSphere(*row[(n_points + s)*dims:(n_points + s + 1)*dims], sphere_radii[s])[d]
    if dims == 2:
        graph.set_data(*frame)
        return graph,
    if dims == 3:
        graph._offsets3d = tuple(frame)
        return ax,

interval = 20 #ms
nframes = int(tf/ (interval* 1e-3))
pb2 = PB.ProgressBar(nframes + 1, bar_length = 10)
ani = animation.FuncAnimation(fig, animate, frames=nframes, interval=interval, blit=False)
#plt.show()
ani.save("c_spring_ivp.mp4")

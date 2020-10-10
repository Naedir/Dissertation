

import qutip as q
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation as ani
import numpy as np

## create an up bit:

up = (q.Qobj([[1 + 1j],[-1j]])).unit()

## add a magnetic field Hamiltonian:
t = np.linspace(0,100,1000)

H_constants = 1
H = H_constants/2 * q.sigmaz()

res = q.mesolve(H, up, t)


figr, ax = plt.subplots()
ax = Axes3D(fig = figr)#, animated = True)
sphere=q.Bloch(axes=ax, fig = figr)
sphere.render(fig=figr, axes=ax)
sphere.add_states(up)
# sphere.show()
plt.pause(0.1)

for i in range(len(res.states)):
    sphere.clear()
    sphere.add_states(res.states[i])
    sphere.make_sphere()
    plt.pause(0.01)


def ini():
    sphere.vector_color=("r")
    return ax

def animate(i):
    sphere.clear()
    sphere.add_states(res.states[i])
    sphere.make_sphere()
    return ax, sphere


ani.FuncAnimation(figr, animate, frames=t, init_func=ini, repeat=False, blit = False)





# def anim(i):
#     b.clear()
#     b.add_states(i)
#     b.make_sphere()
#     return ax



# ##solve the time evolution of the up vector






# # an.FuncAnimation(b, result.states, 10)
# plt.ion()
# fig = plt.figure()

# ax = Axes3D(fig, animated = True)

# b = qt.Bloch(fig=fig, axes=ax)
# b.render(fig=fig, axes=ax)

# b.show()
# plt.pause(0.1)


# # for i in range(len(result.states)):
# #     b.clear()
# #     b.add_states(result.states[i])
# #     b.make_sphere()
# #     b.show()

# # b.add_states(up)
# # b.show()
# result = qt.mesolve(H, up, time)

# for i in range(len(result.states)):
#     anim(result.states[i])
#     fig.pause(0.1)
# # fig.show()
# # an.FuncAnimation(fig, anim(time), frames = time, blit=True)#, repeat=False)
# # anim(time)


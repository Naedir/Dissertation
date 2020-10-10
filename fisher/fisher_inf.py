import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
T1 = 2*10**(-3)
t = np.linspace(0.000,10,100)*10**(-3)


def inf(T1, t):
    b = np.e**(-t/T1)

    return (t**2)/((T1**4)*(b-1))

def infDer(T1,t):
    b = np.e**(-t/T1)

    return (-1/(T1**4)) * ((2*t*T1*(b-1)+(t**2)*b)/(T1*((b-1)**2)))


# T1, t = np.meshgrid(T1, t)

T1 = 2*10**(-3)
t = np.linspace(0.000,10,100)*10**(-3)

I = (t**2)/((T1**4)*(np.e**(t/T1)-1))
plt.plot(t,I)
plt.show()


# plt.plot(t,der)
# plt.plot([0, 0], [min(der), max(der)], '--', linewidth = 0.7)
# plt.plot([min(t), max(t)], [0, 0], '--', linewidth = 0.7)
# plt.legend(['Fisher Info', "derivative w.r.t. 't'"])
# plt.xlabel('time/t')
# plt.title('T1 arbitrarily chosen to be 1')
# plt.tight_layout()
# plt.show()



# der = (-1/(T1**4)) * ((2*t*T1*(b-1)+(t**2)*b)/(T1*((b-1)**2)))


# fig = plt.figure()
# ax = fig.gca(projection='3d')

# surf = ax.plot_surface(T1, t, I, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# ax.set_xlabel('T1')
# ax.set_ylabel('t')
# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.tight_layout()
# plt.show()


# fig, x = plt.subplots(2)

# x[0].plot(t,I)
# x[1].plot(t,der)

# x[0].set_ylabel('Fisher Information')
# x[1].set_ylabel('Derivative of \nFisher Information w.r.t. t')

# x[0].set_xlabel('time / t')
# x[1].set_xlabel('time / t')

# x[0].set_title('T1 = 1')

# x[0].plot([min(t), max(t)], [0, 0], '--', linewidth = 0.7)
# x[1].plot([min(t), max(t)], [0, 0], '--', linewidth = 0.7)

# x[0].plot([0, 0], [min(I), max(I)], '--', linewidth = 0.7)
# x[1].plot([0, 0], [min(der), max(der)], '--', linewidth = 0.7)
# # plt.plot([min(t), max(t)],[0, 0])
# # plt.plot(t,I)
plt.plot(t,I)
# plt.plot(t,der)
# plt.plot([0, 0], [min(der), max(der)], '--', linewidth = 0.7)
# plt.plot([min(t), max(t)], [0, 0], '--', linewidth = 0.7)
plt.legend(['Fisher Info', "derivative w.r.t. 't'"])
plt.xlabel('time/t')
plt.title('T1 arbitrarily chosen to be 1')
plt.tight_layout()
plt.show()

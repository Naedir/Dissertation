import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

t = np.linspace(1, 50,500)*(10**(-3))
T1 = np.linspace(1, 10, 500)*(10**(-3))

#X is t, and Y is T1

X, Y = np.meshgrid(t, T1)

def ticks(lists, number):
    n = len(lists)
    skip = int(n/number)
    arr = [None]*number
    ind = 0
    for i in range(number):
        arr[i] = lists[ind]
        ind += skip
    return arr


I = (X**2)/((Y**4)*(np.e**(X/Y)-1))

I = np.log10(I)

fig, ax = plt.subplots(1,1)

img = ax.imshow(I)

ax.set_xticks(list(np.linspace(0,len(t),10)))
ax.set_xticklabels(np.round(ticks(t*1000, 10), 2))
ax.set_yticks(list(np.linspace(0,len(T1),10)))
ax.set_yticklabels(np.round(ticks(T1*1000, 10), 2))

plt.show()


# fig = plt.figure()
# ax = fig.gca(projection='3d')

# surf = ax.plot_surface(X, Y, I, cmap='plasma', linewidth=0, antialiased=False)
# ax.set_xticklabels(np.round(t*1000,2))
# ax.set_yticklabels(np.round(T1*1000,2))
# ax.set_xlabel('$T1*10^{-6}$')
# ax.set_ylabel('$t*10^{-6}$')
# ax.set_title('Fisher Information in log$_{10}$ scale')
# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.tight_layout()
# plt.show()

# plt.imshow(I)
# plt.colorbar(shrink=0.5, aspect=5)
# # plt.xticklabels(np.round(t*1000,2))
# # plt.yticklabels(np.round(T1*1000,2))
# plt.xlabel('$T1*10^{-3}$')
# plt.ylabel('$t*10^{-3}$')
# plt.title('Fisher Information in log$_{10}$ scale')
# plt.show()

# # plt.hist2d(t,T1)#,[len(t), len(T1)])

# # ax.set_yticklabels(vegetables)

# plt.show()

# # fig = plt.figure()
# # ax = fig.gca(projection='3d')

# # surf = ax.plot_surface(Y, X, I, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# # ax.set_xlabel('T1')
# # ax.set_ylabel('t')
# # fig.colorbar(surf, shrink=0.5, aspect=5)

# # plt.tight_layout()
# # plt.show()
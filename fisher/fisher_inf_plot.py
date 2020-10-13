import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

def fisher_inf(t, T1):

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

    I2 = np.log10(I)

    fig, ax = plt.subplots(1,2)

    img = ax[0].imshow(I)

    img = ax[1].imshow(I2)

    ax[0].set_xticks(list(np.linspace(0,len(t),10)))
    ax[0].set_xticklabels(np.round(ticks(t*1000, 10), 2), rotation = 45)
    ax[0].set_yticks(list(np.linspace(0,len(T1),10)))
    ax[0].set_yticklabels(np.round(ticks(T1*1000, 10), 2))
    ax[0].set_title('Fisher Information')
    ax[0].set_ylabel('T1')
    ax[0].set_xlabel('t')

    ax[1].set_xticks(list(np.linspace(0,len(t),10)))
    ax[1].set_xticklabels(np.round(ticks(t*1000, 10), 2), rotation = 45)
    ax[1].set_yticks(list(np.linspace(0,len(T1),10)))
    ax[1].set_yticklabels(np.round(ticks(T1*1000, 10), 2))
    ax[1].set_title('Fisher Information \nin a $log_{10}$ scale')
    ax[1].set_ylabel('T1')
    ax[1].set_xlabel('t')

    fig.colorbar(img, ax = ax.ravel().tolist())
    # fig.tight_layout()

    plt.show()

if __name__ == '__main__':
    t = np.linspace(1, 50,500)*(10**(-3))
    T1 = np.linspace(1, 10, 500)*(10**(-3))

    fisher_inf(t, T1)
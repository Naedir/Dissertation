from matplotlib import pyplot as pyplot
import numpy as np
from experiment import *
import pandas as pd
import sys, os
from subprocess import call
         
pathname = os.path.dirname(sys.argv[0])        


n = 1000 #number of experiments in each run
reps = 1000   #number of times the experiments will be repeated


matrix_a = [] #this will store all the results
matrix_b = []
matrix_c = []


times = np.linspace(1e-9,20e-3,n)

for i in range(reps):
        

    a = experiment(10e-3, 5e-3, times)

    a.loop_1()
    matrix_a.append(a.T1_estimates)
    # a = experiment(10e-3, 5e-3, times)
    a.loop_2(n)
    matrix_b.append(a.T1_estimates)
    # a = experiment(10e-3, 5e-3, times)
    a.loop_3(n)
    matrix_c.append(a.T1_estimates)
    call('clear')
    print("\n\n")
    print(f'{i/reps} %')
    print("\n\n")

matrix_a = np.mean(matrix_a, axis = 0)
matrix_b = np.mean(matrix_b, axis = 0)
matrix_c = np.mean(matrix_c, axis = 0)


fig, ax = plt.subplots(3)

def plot_ax(i, data):
    ax[i].plot(data)
    ax[i].plot(np.linspace(1,n+1,n), [a.T1]*n)

    ax[i].set_ylabel('T1')
    ax[i].set_xlabel('experiment')

plot_ax(0, matrix_a)
plot_ax(1, matrix_b)
plot_ax(2, matrix_c)

ax[0].set_title('approach A')
ax[0].legend(['T1 estimate', 'true T1'])
    
ax[1].set_title('approach B')

ax[1].legend(['T1 estimate', 'true T1'])

ax[2].set_title('approach C')

ax[2].legend(['T1 estimate', 'true T1'])


fig.suptitle(f'Comparison of 3 different approaches averaged afrer {reps} data sets')
fig.tight_layout()

plt.show()

if input('save results? (y/n):\n')=='y':
    import os 
    import sys
    path1 = os.getcwd() + '/data/matrix_1a.csv'
    path2 = os.getcwd() + '/data/matrix_1b.csv'
    path3 = os.getcwd() + '/data/matrix_1c.csv'

    pd.DataFrame(matrix_a).to_csv(path1)
    pd.DataFrame(matrix_b).to_csv(path2)
    pd.DataFrame(matrix_c).to_csv(path3)
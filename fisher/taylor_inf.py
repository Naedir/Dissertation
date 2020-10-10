from taylor import taylor
from sympy import *
import matplotlib.pyplot as plt
import numpy as np

# Factorial function



t, T1 = symbols('t T1')

f = t**2/(T1**4*(exp(t/T1) - 1))
f_lam = lambdify([T1, t], f)

f_taylor = taylor(f, t, T1, 2)

#get values for range of t and a few T1

#T1 = {1ms, 5ms, 10ms, 15ms}
#t = {0:50 us}


T1_nums = [1e-3, 5e-3, 10e-3, 15e-3]
t_nums = np.linspace(1e-6,50e-3, 50)

res = [None]*len(T1_nums)
res2 = [None]*len(T1_nums)
ind1 = 0
for i in T1_nums:
    ind2 = 0
    y = [None]*len(t_nums)
    y2 = [None]*len(t_nums)
    for k in t_nums:
        y[ind2] = f_taylor.evalf(subs={T1: i, t:k})#([(T1, i), (t, k)]))
        y2[ind2] = f_lam(i, k)
        ind2 += 1
    res[ind1] = y
    res2[ind1] = y2
    ind1 += 1


fig, ax = plt.subplots(2,2)

ind = 0
for i in range(len(ax)):
    for k in range(len(ax[i])):
        ax[i][k].plot(t_nums, res2[ind])
        ax[i][k].plot(t_nums, res[ind])
        ax[i][k].plot([1.462117157 * T1_nums[ind], 1.462117157 * T1_nums[ind]], [min(res[ind]), max(res[ind])], '--', color = 'green')
        ax[i][k].set_title(f"T1 = {T1_nums[ind]}")
        ax[i][k].legend(['func', 'approx. func', 'maximum'])

        ind += 1


fig.tight_layout()
fig.show()
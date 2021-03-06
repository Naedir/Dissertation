import taylor
import sympy as sp
from sympy import *
import matplotlib.pyplot as plt
import numpy as np

'''

Plotting fisher information at given T1 times and a range of t times
along with a taylor expansion around t = T1

a maximum is shown using a dashed line
'''

def max_fisher_info(t_nums, T1_nums):
    t, T1 = symbols('t T1')

    f = (t**2*(1+sp.exp(-t/T1)))/(2*T1**4*(sp.exp(t/T1)+1)**2) + (t**2*(1-sp.exp(-t/T1)))/(2*T1**4*(sp.exp(t/T1)-1)**2)
    f_lam = lambdify([T1, t], f)

    f_taylor = taylor(f, t, T1, 2)

    res = [None]*len(T1_nums)
    res2 = [None]*len(T1_nums)
    ind1 = 0
    for i in T1_nums:
        ind2 = 0
        y = [None]*len(t_nums)
        y2 = [None]*len(t_nums)
        for k in t_nums:
            y[ind2] = f_taylor.evalf(subs={T1: i, t:k})
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
            ax[i][k].plot([0.73425356207799* T1_nums[ind], 0.73425356207799 * T1_nums[ind]], [min(res[ind]), max(res[ind])], '--', color = 'green')
            ax[i][k].set_title(f"T1 = {T1_nums[ind]}")
            ax[i][k].legend(['func', 'approx. func', 'maximum'])

            ind += 1


    fig.tight_layout()
    plt.show()
def factorial(n):
    if n < 2:
        return 1
    else:
        return n * factorial(n-1)

def taylor(func, variable, a, order):
    f1 = func.subs(variable, a)
    f_der = {}
    for i in range(1, order+1):
        f_der[f"order{i}"] = (1/(factorial(i))) * ((diff(func, variable, i)).subs(variable, a))*((variable - a)**i)

    result = 0
    # if order > 1:
    for i in list(f_der.values()):
        result += i

    return f1 + result


if __name__ == '__main__':
    T1_nums = [1e-3, 5e-3, 10e-3, 15e-3]
    t_nums = np.linspace(1e-6,50e-3, 50)

    max_fisher_info(t_nums, T1_nums)

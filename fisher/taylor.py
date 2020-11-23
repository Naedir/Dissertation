from matplotlib import pyplot as plt
import numpy as np
from sympy import *
import time as tm



class taylor_exp:

    """
    A class that calculates a Taylor expansion of a given function

    Attributes:
        fun: function (symppy function fype)
        order: order of the approsimation
        around: around which point to evaluate the function
        legend_lim: set a limit to which order to draw the legend (useful when making animations)
        variable: function variable (usually 'x') as a sympy variable type
        xrange: range of x values to evaluate the function at. list or numpy.array type

    Methods:
        plot_func: plots a function along with it's Taylor expansion plots
        anim_func: makes an animated plot of a function along with its Taylor expansion

    """
    def __init__(self, fun, order, around, variable, xrange, legend_lim = False):
        self.fun = fun
        self.order = order
        self.around = around
        self.legend_lim = legend_lim
        self.variable = variable
        self.xrange = xrange
        init_printing(use_unicode=True)


    def factorial(self, n):
        if n < 2:
            return 1
        else:
            return n * factorial(n-1)
    
    def eval(self):
        lam_f = lambdify(self.variable, self.fun, 'numpy')

        leg = [self.fun]
        fig, ax = plt.subplots(1,1)
        ind2 = 0
        for i in range(1, self.order + 1):
            f_taylor = self.taylor(i)
            y = [None]*len(self.xrange)
            ind = 0
            for i in self.xrange:
                y[ind] = f_taylor.subs(self.variable, i)
                ind +=1
            ind2 += 1


    def taylor(self):
        f1 = self.fun.subs(self.variable, self.around)
        f_der = {}
        for i in range(1, self.order+1):
            f_der[f"order{i}"] = (1/(self.factorial(i))) * ((diff(self.fun, self.variable, i)).subs(self.variable, self.around))*((self.variable - self.around)**i)

        result = 0
        for i in list(f_der.values()):
            result += i

        self.solution = f1 + result
        return f1 + result

    def plot_func(self):
        lam_f = lambdify(self.variable, self.fun, 'numpy')

        leg = [f"${latex(self.fun)}$"]
        fig, ax = plt.subplots(1,1)
        ax.set_title("Taylor expansion of ${}$ around {}".format(latex(self.fun), self.around))
        ax.plot(self.xrange, lam_f(self.xrange))#,self.variable)#, 'x', animated = True)
        ind2 = 0
        for i in range(1, self.order + 1):
            self.order = i
            f_taylor = self.taylor()
            y = [0]*len(self.xrange)
            ind = 0
            for i in self.xrange:
                y[ind] = f_taylor.subs(self.variable, i)
                ind +=1
            ax.plot(self.xrange, y)
            ind2 += 1
            n = "order {}".format(ind2)
            if self.legend_lim == False:
                leg.append(n)

            else:
                if ind2<=self.legend_lim:
                    leg.append(n)
                    ax.legend(leg)
                else:
                    pass
            
            ax.legend(leg)        
        fig.show()
    
    def anim_func(self):

        lam_f = lambdify(self.variable, self.fun, 'numpy')

        leg = [self.fun]
        plt.ion()
        fig, ax = plt.subplots(1,1)
        ax.set_title("Taylor expansion of {} around {}".format(self.fun, self.around))
        ax.plot(self.xrange, lam_f(self.xrange),'x')#, 'x', animated = True)

        plt.pause(0.5)
        plt.block = True

        ind2 = 0

        for i in range(1, self.order + 1):
            self.order = i
            f_taylor = self.taylor()
            y = [None]*len(self.xrange)
            ind = 0
            for i in self.xrange:
                y[ind] = f_taylor.subs(self.variable, i)
                ind +=1

            ax.plot(self.xrange, y)

            ind2 += 1
            n = "order {}".format(ind2)
            if self.legend_lim == False:
                leg.append(n)

            else:
                if ind2<=self.legend_lim:
                    leg.append(n)
                    ax.legend(leg)
                else:
                    pass
            
            ax.legend(leg)    
            ax.set_title("Taylor expansion of {} around {}".format(self.fun, self.around))                
            plt.pause(1)


### here just the taylor expansion function:


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

import sympy as sp

t, T1 = sp.symbols('t T1')
t_1 = 0.001
f = t**2/(T1**4*(exp(2*t/T1) - 1))
# g = t**3
# f = f.subs(t,t_1)
ts = np.linspace(0.1e-3, 50e-3,1000)
# ts = np.linspace(-5,5,100)
# a = taylor_exp(f, 0, t, T1, ts)

def t(ord):
    t, T1 = sp.symbols('t T1')
    t_1 = 0.001
    f = t**2/(T1**4*(exp(2*t/T1) - 1)) 

    g = taylor(f,T1, t, ord)

    f = f.subs(t,t_1)
    g = g.subs(t,t_1)

    fl = lambdify(T1, f)
    gl = lambdify(T1, g)

    return ord, fl, gl
# a.plot_func()
a, b, c = t(3)
plt.plot(ts, b(ts))
plt.plot(ts, c(ts), '--')
plt.legend(['FI', f'approximation order 3'])
a, b, c = t(4)
# plt.plot(b(ts))
plt.plot(ts, c(ts), '--')
plt.legend(['FI', f'approximation order 4'])
a, b, c = t(5)
# plt.plot(b(ts))
plt.plot(ts, c(ts), '--')
plt.legend(['FI', f'approximation order 3', f'approximation order 4', f'approximation order 5'])
plt.xlabel('T1')
plt.title(f't = {t_1}')
plt.ylim([-100, 500000])
plt.show()


# taylor_exp()
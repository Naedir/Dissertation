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
        return f1 + result

    def plot_func(self):
        lam_f = lambdify(self.variable, self.fun, 'numpy')

        leg = [self.fun]
        fig, ax = plt.subplots(1,1)
        ax.set_title("Taylor expansion of {} around {}".format(self.fun, self.around))
        ax.plot(self.xrange, lam_f(self.xrange),'x')#, 'x', animated = True)
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

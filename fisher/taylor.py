from matplotlib import pyplot as plt
import numpy as np
from sympy import *
import time as tm

init_printing(use_unicode=True)


# Taylor::
# f ~~ f(a) + f'(a)(x-a) + f''(a)/2 * (x-a)^2


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
    x = symbols('x')

    f = x**3 + 5 + sin(x)*cos(x)

    nums = np.linspace(-1*np.pi,1*np.pi,50)
    lam_f = lambdify(x, f, 'numpy')
    around = 1

    leg = [f]
    plt.ion()
    fig, ax = plt.subplots(1,1)
    # fig.show()


    ax.set_title("Taylor expansion of {} around {}".format(f, around))


    ax.plot(nums, lam_f(nums),'x')#, 'x', animated = True)

    ax.set_ylim([-5, 20])
    ax.set_xlim([min(nums), max(nums)])

    plt.pause(0.5)
    plt.block = True

    ord = 15
    ind2 = 0




    for i in range(1, ord + 1):
        

        f_taylor = taylor(f, x, around, i)
        y = [None]*len(nums)
        ind = 0
        for i in nums:
            y[ind] = f_taylor.subs(x, i)
            ind +=1

        ax.plot(nums, y)#, animated = True)

        ind2 += 1
        n = "order {}".format(ind2)
        if ind2<=4:
            leg.append(n)
            ax.legend(leg)

        else:
            leg = []
        plt.pause(0.1)


        
        # fig.draw()

    plt.ioff()

    ax.legend(leg)    
    ax.set_title("Taylor expansion of {} around {}".format(f, around))

# fig.show()
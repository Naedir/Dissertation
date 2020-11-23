import numpy as np
import matplotlib.pyplot as plt
from operator import add
from subprocess import call


class experiment:
    '''This class simulates an experiment performed on NV centres. NV centre qubit is coupled to a magnetic field which causes it to precess around z axis.
    The qubit is a subject to decoherence which affects the amount of information that's possible to extract by performing the measurement.

    Properties:
        
        T1: true value of T1
        T1_est: current estimation of a decoherence time
        T1_span: span of T1 times for T1 estimation
        sigma: used as a standard deviation for T1_estimation
        mes_time: optimal time for taking a measurement which will (hopefully) maximise amount of information
        time_span: span of times for taking measurements

    Methods:
        
        gauss: creates a gaussian distribution around estimation of T1 with a given standard deviation sigma
        loop_1: simulates an experiment. Fisher information maximisation is NOT used in this experiment, rather it is a simple sweep of measurements at given list of times
        loop_1_stat: makes statistics of performing loop_1 experiment n times

        loop_2: simulates an experiment n times. At each time Fisher information is being calculated and time with maximum FI is being used for the next measurement.
        loop_2_stat: makes statistics of performing loop_1 experiment n_stat times

        loop_3: simulates an experiment n times. At each measurement the oprimal time for taking the measurement is calculated using the approximate formula found using Taylor expansion.
        loop_3_stat: makes statistics of performing loop_1 experiment n_stat times


        F_I: calculates Fisher infirmation given the T1 estimation and a range of times (self.time_span)
        
    '''
    def __init__(self, T1_est, sigma, time_span, T1_span = np.linspace(1e-9,200e-3,500), T1 = 20e-3):

        self.T1 = T1    #true value of T1
        self.T1_span = T1_span
        self.time_span = time_span
        self.T1_stats = self.uniform_prior()
        self.T1_est = np.mean(self.T1_stats)    #current estimation of T1
        self.T1_estimates = []  #all estimates of T1
        self.T1_estimates.append(T1_est)
        
        self.sigma = sigma
        self.mes_time = 0.73425356207799 * self.T1_est  #optimal time for taking the next measurement
        
        
        self.mes_true_time = 1e-9  #keep a track of the actual time of the measurement

    def gauss(self):
        return (1/(self.sigma * np.sqrt(2 * np.pi)) * np.exp( - (self.T1_span - self.T1_est)**2 / (2 * self.sigma**2)))
    def uniform_prior(self):
        pr = [1]*len(self.T1_span)
        re_normalised = [i/sum(pr) for i in pr]
        return re_normalised

    def loop_1(self): #we just sweep t non-adaptively from a min value to a max value
        #reinitialize T1_estimation
        self.T1_est = np.mean(self.T1_stats)
        self.T1_estimates = [self.T1_est]
        self.T1_stats = self.uniform_prior()
        #the results of the measurement will be rither 1 or 0
        outcomes = [None]*len(self.time_span)    #store outcomes of the measurement

        stats_matrix = [None]*len(self.time_span)

        for i in range(len(self.time_span)):
            #calculate the probability of getting a spin up (p(m=0)):
            self.mes_true_time = self.time_span[i]
            prob = self.p_0(self.time_span[i], self.T1)
            #get the outcome of the measurement using a random number and a p(m=0):
            if np.random.random() < prob:
                outcomes[i] = 0
            else:
                outcomes[i] = 1

            #update T1 using bayesian update
            self.bayesian_update(outcomes[i])
            self.T1_est = self.T1_span[np.nanargmax(self.T1_stats)]
            self.T1_estimates.append(self.T1_est)
            stats_matrix[i] = self.T1_stats

        return outcomes, stats_matrix

    def loop_1_stat(self, n):   #this will just use loop_1 and repeat it n times to make some statistics
        stats = self.loop_1()
        for i in range(n - 1):
            stats = list(map(add, stats, self.loop_1()))

        return [number / n for  number in stats]    #normalise the results by dividing each bin by 'n'

    def loop_2(self, n):    #at each point, you compute FI numerically and find its maximum
        #reinitialize T1_estimation
        self.T1_est = np.mean(self.T1_span)
        self.T1_estimates = [self.T1_est]
        self.T1_stats = self.uniform_prior()
        outcomes = [None]*n    #store outcomes of the measurement
        stats_matrix = [None]*n
        for i in range(n):
            index = np.nanargmax(self.F_I()) #finds time at which fisher information is at maximum
            self.mes_true_time = self.T1_span[index]
            #calculate the probability of getting a spin up (p(m=0)):
            prob = self.p_0(self.mes_true_time, self.T1)
            #get the outcome of the measurement using a random number and a p(m=0):
            if np.random.random() < prob:
                outcomes[i] = 0
            else:
                outcomes[i] = 1
            #update T1 using bayesian update
            self.bayesian_update(outcomes[i])
            self.T1_est = self.T1_span[np.nanargmax(self.T1_stats)]
            self.T1_estimates.append(self.T1_est)
            stats_matrix[i] = self.T1_stats
        return outcomes, stats_matrix

    def loop_2_stat(self, n, n_stat):   #this will just use loop_2 and repeat it n_stat times to make some statistics
        stats = self.loop_2(n)
        for i in range(n - 1):
            stats = list(map(add, stats, self.loop_2(n)))
        return [number / n for  number in stats] 

    def loop_3(self, n):    #use my formula for the optimal time of the measurement found using Taylor's expansion
        #reinitialize T1_estimation
        self.T1_est = np.mean(self.T1_stats)
        self.T1_estimates = [self.T1_est]
        self.T1_stats = self.uniform_prior()
        outcomes = [None]*n    #store outcomes of the measurement
        stats_matrix = [None]*n
        for i in range(n):
            time = 0.73425356207799 * self.T1_est
            self.mes_true_time = time
            #calculate the probability of getting a spin up (p(m=0)):
            prob = self.p_0(time, self.T1)
            #get the outcome of the measurement using a random number and a p(m=0):
            if np.random.random() < prob:
                outcomes[i] = 0
            else:
                outcomes[i] = 1
            
            #update T1 using bayesian update
            self.bayesian_update(outcomes[i])
            self.T1_est = self.T1_span[np.nanargmax(self.T1_stats)]
            self.T1_estimates.append(self.T1_est)
            stats_matrix[i] = self.T1_stats
        return outcomes, stats_matrix

    def loop_3_stat(self, n, n_stat):   #this will just use loop_2 and repeat it n_stat times to make some statistics
        #reset the T1 distribution:
        self.T1_stats = self.gauss()
        stats = self.loop_3(n)
        for i in range(n - 1):
            stats = list(map(add, stats, self.loop_3(n)))
        return [number / n for  number in stats] 

    def F_I(self):  #calculate Fisher information at a given T1
        t = self.T1_span
        return (t**2)*(1+np.exp(-t/self.T1_est))/(2*self.T1_est**4*(np.exp(t/self.T1_est)+1)**2) + (t**2)*(1-np.exp(-t/self.T1_est))/(2*self.T1_est**4*(np.exp(t/self.T1_est)-1)**2)

    def num_integrate(self, fun, span, var):
        #function is a prior estimate multiplied by a Fisher infirmation
        fun = self.stats * self.F_I()
        # step = (span[len(span)-1] - span[0])/len(span)
        # cumsum = 0
        # intgr = [None] * len(span)
        # ind = 0
        # for i in span:
        #     cumsum = cumsum + step * fun.subs(var, i)
        #     intgr[ind] = cumsum
        #     ind += 1
        # return intgr

    def p_0(self, t, T1):
        return ((1 - np.exp(-t/T1))/2)

    def p_1(self, t, T1):
        return ((1 + np.exp(-t/T1))/2)

    def bayesian_update(self, m):

        measurement_time = self.mes_true_time
        ind = 0
        if m == 0:
            for i in self.T1_span:
                self.T1_stats[ind] = self.T1_stats[ind]*self.p_0(measurement_time, i)
                ind+=1
        elif m == 1:
            for i in self.T1_span:
                self.T1_stats[ind] = self.T1_stats[ind]*self.p_1(measurement_time, i)
                ind+=1
        else:
            raise ValueError("Measurement outcome has to be 1 or 0")

        # #normalize data:
        norm = sum(self.T1_stats)
        self.T1_stats = [i/norm for i in self.T1_stats]

if __name__ == '__main__':  #some plots, delete if necessary

    from sympy import *


    # def num_integrate(fun, span, var):
    #     if span[0] != 0:
    #         st = np.linspace(0, span[0], 10)
    #         step = st[2] - st[1]
    #         cumsum = 0
    #         for i in st:
    #             cumsum += step * fun.subs(var, i)
    #     else:
    #         cumsum = 0
    #     step = span[2] - span[1]
    #     intgr = [None] * len(span)
    #     ind = 0
    #     for i in span:
    #         cumsum += step * fun.subs(var, i)
    #         intgr[ind] = cumsum
    #         ind += 1
    #     return cumsum

    

    # t, T1 = symbols('t, T1')

    # f = t**2/(T1**4*(exp(2*t/T1)-1))
    # ts = np.linspace(1e-9, 200e-3, 100)
    # d = num_integrate(f, ts, T1)
    # lam_d = lambdify(t,d)
    # plt.plot(ts, lam_d(ts))
    # plt.show()
    # s = np.linspace(0,4,100)

    # g = integrate(f,x)
    # d = num_integrate(f,s,x)


    # d = trapezoidal(y,-2, 2, 100)
    

    # lam_f = lambdify(x,f)
    # lam_g = lambdify(x,g)
    
    # # d = np.cumsum(lam_f(s))*(s[1] - s[0])

    # plt.plot(s, lam_f(s))
    # plt.plot(s, lam_g(s))
    # plt.plot(s, d, '--')
    # plt.show()
    # # plt.close("all")
    
    n = 100 #number of experiments
    n_stats = 100#number of repeats of the experiment for the statistics
    n_ts = 100
    times = np.linspace(1e-9,200e-3,n)
    Ts = np.linspace(1e-4, 80e-3, n_ts)
    # # Ts = [20e-3]

    r1 = [np.array([0.0]*(n+1))]*n_stats
    r2 = [np.array([0.0]*(n+1))]*n_stats
    r3 = [np.array([0.0]*(n+1))]*n_stats
    matrix_a = [] #this will store all the results
    matrix_b = []
    matrix_c = []
    i = 0
    for ins in Ts:
        
        pr = (i/len(Ts))*100
        call('clear')
        print("\n\n")
        print(f'{pr} %')
        print("\n\n")

        for j in range(n_stats):
                
            a = experiment(10e-3, 5e-3, times, T1 = ins)
            a.loop_1()
            matrix_a.append(a.T1_estimates)
            r1[j] += np.abs(a.T1_estimates - a.T1)**2
            call('clear')
            print("\n\n")
            print(f'{pr} %')
            print("\n\n")

            a = experiment(10e-3, 5e-3, times, T1 = ins)
            a.loop_2(n)
            matrix_b.append(a.T1_estimates)
            r2[j] += np.abs(a.T1_estimates - a.T1)**2
            call('clear')
            print("\n\n")
            print(f'{pr} %')
            print("\n\n")

            a = experiment(10e-3, 5e-3, times, T1 = ins)
            a.loop_3(n)
            matrix_c.append(a.T1_estimates)
            r3[j] += np.abs(a.T1_estimates - a.T1)**2
            call('clear')
            print("\n\n")
            print(f'{pr} %')
            print("\n\n")
        i+=1

    # import numpy as np
    # from matplotlib import pyplot as plt
    # from matplotlib import animation
    # if 0:
            
    #     # First set up the figure, the axis, and the plot element we want to animate
    #     # fig = plt.figure()
    #     # ax = plt.axes()
    #     # fps = 30
    #     # name = 'approach_A_1.mp4'
    #     # def anim1(i):
    #     #     ax.clear()
    #     #     ax.text(0.045, 0.02, f'exp: {int(i/100)}*100')
    #     #     ax.text(0.045, 0.015, f"error = {round(r1[i], 7)}")
    #     #     ax.plot([a.T1]*2, [0, 0.023])
    #     #     ax.plot(a.T1_span, st[i])
    #     #     ax.set_ylim([0, 0.023])
    #     #     return ax
        
    #     # frames = np.linspace(0,n - 1,n, dtype=int)
    #     # anim = animation.FuncAnimation(fig, anim1, frames=frames, blit=False)
    #     # anim.save(name, fps = fps)
    #     # print("saved anim 1")

    #     # fig = plt.figure()
    #     # ax = plt.axes()
    #     # name = 'approach_B_2.mp4'
    #     # def anim1(i):
    #     #     ax.clear()
    #     #     ax.text(0.045, 0.02, f'exp: {int(i/100)}*100')
    #     #     ax.text(0.045, 0.015, f"error = {round(r2[i], 7)}")
    #     #     ax.plot([a.T1]*2, [0, 0.023])
    #     #     ax.plot(a.T1_span, st2[i])
    #     #     ax.set_ylim([0, 0.023])
    #     #     return ax
        
    #     # frames = np.linspace(0,n - 1,n, dtype=int)
    #     # anim = animation.FuncAnimation(fig, anim1, frames=frames, blit=False)
    #     # anim.save(name, fps = fps)
    #     # print("saved anim 2")

    #     # fig = plt.figure()
    #     # ax = plt.axes()
    #     # name = 'approach_C_1.mp4'
    #     # def anim1(i):
    #     #     ax.clear()
    #     #     ax.text(0.045, 0.02, f'exp: {int(i/100)}*100')
    #     #     ax.text(0.045, 0.015, f"error = {round(r3[i],7)}")
    #     #     ax.plot([a.T1]*2, [0, 0.023])
    #     #     ax.plot(a.T1_span, st3[i])
    #     #     ax.set_ylim([0, 0.023])
    #     #     return ax
        
    #     # frames = np.linspace(0,n - 1,n, dtype=int)
    #     # anim = animation.FuncAnimation(fig, anim1, frames=frames, blit=False)
    #     # anim.save(name, fps = fps)
    #     # print("saved anim 3")
    #     pass

    plt.close("all") 
    r1 = np.mean(r1, axis = 0)[1::]
    r2 = np.mean(r2, axis = 0)[1::]
    r3 = np.mean(r3, axis = 0)[1::]
    r1 = r1/n_ts
    r2 = r2/n_ts
    r3 = r3/n_ts
    # r1 = r1[1::]
    # r2 = r2[1::]
    # r3 = r3[1::]
    matrix_a = np.mean(matrix_a, axis = 0)[1::]
    matrix_b = np.mean(matrix_b, axis = 0)[1::]
    matrix_c = np.mean(matrix_c, axis = 0)[1::]

    fig, ax = plt.subplots()
    # ax.plot(np.linspace(1,n+1,n),[a.T1]*n); ax.plot(np.linspace(1,n+1,n), matrix_a, color = 'red'); ax.plot(np.linspace(1,n+1,n), matrix_b); ax.plot(np.linspace(1,n+1,n), matrix_c, linestyle='--')
    # plt.title(f'Comparison of 3 different approaches averaged afrer {n_stats} data sets')
    # plt.xlabel('Experiment')
    # plt.ylabel('T1')
    # ax.legend(['true T1', 'approach A estimate', 'approach B estimate','approach C estimate'])
    # plt.show()

    ax.plot(np.linspace(1,n+1,n), [0]*n); ax.plot(np.linspace(1,n+1,n), r1, color = 'red'); ax.plot(np.linspace(1,n+1,n), r2); ax.plot(np.linspace(1,n+1,n), r3, linestyle = '--')
    # plt.ylabel('error')
    # plt.xlabel('experiment number')
    # plt.title('T1 est. error')
    ax.legend(['0', 'approach A','approach B','approach C'])
    ax.set_xlabel("measurement number")
    ax.set_ylabel("Absolute error")
    ax.set_title("Error for a range if T1 times")
    # # plt.ylim([0, 0.0009])     
    plt.show()


    # # if input('save results? (y/n):\n')=='y':
    #     # import os 
    #     # import sys
    #     # import pandas as pd
    #     # path1 = '/home/naedir/Desktop/Project/github/Krzysztof Skrzypczak/Simulations' + '/data/matrix_T1_dist_a.csv'
    #     # path2 = '/home/naedir/Desktop/Project/github/Krzysztof Skrzypczak/Simulations' + '/data/matrix_T1_dist_b.csv'
    #     # path3 = '/home/naedir/Desktop/Project/github/Krzysztof Skrzypczak/Simulations' + '/data/matrix_T1_dist_c.csv'

    #     # pd.DataFrame(r1).to_csv(path1)
    #     # pd.DataFrame(r2).to_csv(path2)
    #     # pd.DataFrame(r3).to_csv(path3)


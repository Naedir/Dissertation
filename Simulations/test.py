import numpy as np
import matplotlib.pyplot as plt
from operator import add



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
    def __init__(self, T1_est, sigma, time_span, T1_span = np.linspace(0,50e-3,1000), T1 = 20e-3):

        self.T1 = T1    #true value of T1
        self.T1_est = T1_est    #current estimation of T1
        self.T1_estimates = []  #all estimates of T1
        self.T1_estimates.append(T1_est)
        self.T1_span = T1_span
        self.sigma = sigma
        self.mes_time = 0.73425356207799 * self.T1_est  #optimal time for taking the next measurement
        self.time_span = time_span
        self.T1_stats = self.gauss()
        self.mes_true_time = 0  #keep a track of the actual time of the measurement

    def gauss(self):
        return (1/(self.sigma * np.sqrt(2 * np.pi)) * np.exp( - (self.T1_span - self.T1_est)**2 / (2 * self.sigma**2)))

    def loop_1(self): #we just sweep t non-adaptively from a min value to a max value
        #reinitialize T1_estimation
        self.T1_est = self.T1_estimates[0]
        self.T1_estimates = [self.T1_est]
        #the results of the measurement will be rither 1 or 0
        outcomes = [None]*len(self.time_span)    #store outcomes of the measurement

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
        return outcomes

    def loop_1_stat(self, n):   #this will just use loop_1 and repeat it n times to make some statistics
        stats = self.loop_1()
        for i in range(n - 1):
            stats = list(map(add, stats, self.loop_1()))

        return [number / n for  number in stats]    #normalise the results by dividing each bin by 'n'

    def loop_2(self, n):    #at each point, you compute FI numerically and find its maximum
        #reinitialize T1_estimation
        self.T1_est = self.T1_estimates[0]
        self.T1_estimates = [self.T1_est]
        self.T1_stats = self.gauss()
        outcomes = [None]*n    #store outcomes of the measurement
        for i in range(n):
            index = np.nanargmax(self.F_I()) #finds time at which fisher information is at maximum
            self.mes_true_time = self.time_span[index]
            #calculate the probability of getting a spin up (p(m=0)):
            prob = self.p_0(self.time_span[index], self.T1)
            #get the outcome of the measurement using a random number and a p(m=0):
            if np.random.random() < prob:
                outcomes[i] = 0
            else:
                outcomes[i] = 1
            #update T1 using bayesian update
            self.bayesian_update(outcomes[i])
            self.T1_est = self.T1_span[np.nanargmax(self.T1_stats)]
            self.T1_estimates.append(self.T1_est)
        return outcomes

    def loop_2_stat(self, n, n_stat):   #this will just use loop_2 and repeat it n_stat times to make some statistics
        stats = self.loop_2(n)
        plt.show()
        for i in range(n - 1):
            stats = list(map(add, stats, self.loop_2(n)))
        return [number / n for  number in stats] 

    def loop_3(self, n):    #use my formula for the optimal time of the measurement found using Taylor's expansion
        #reinitialize T1_estimation
        self.T1_est = self.T1_estimates[0]
        self.T1_estimates = [self.T1_est]
        self.T1_stats = self.gauss()
        outcomes = [None]*n    #store outcomes of the measurement
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
        return outcomes

    def loop_3_stat(self, n, n_stat):   #this will just use loop_2 and repeat it n_stat times to make some statistics
        #reset the T1 distribution:
        self.T1_stats = self.gauss()
        stats = self.loop_3(n)
        for i in range(n - 1):
            stats = list(map(add, stats, self.loop_3(n)))
        return [number / n for  number in stats] 

    def F_I(self):  #calculate Fisher information at a given T1
        t = self.time_span
        return (t**2)*(1+np.exp(-t/self.T1_est))/(2*self.T1_est**4*(np.exp(t/self.T1_est)+1)**2) + (t**2)*(1-np.exp(-t/self.T1_est))/(2*self.T1_est**4*(np.exp(t/self.T1_est)-1)**2)

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

        
        #normalize data:
        norm = sum(self.T1_stats)
        self.T1_stats = [i/norm for i in self.T1_stats]
        if ind%50 == 0:
            plt.plot(self.T1_span, self.T1_stats)
            plt.pause(0.01)

print('\n\n\nRunning Test\n\n')

if __name__ == '__main__':  #some plots, delete if necessary



    plt.close("all")
    
    n = 100 #number of experiments
    n_stat = 1#number of repeats of the experiment for the statistics

    n_ax = [i+1 for i in range(n)]
    times = np.linspace(1e-9,20e-3,n)
    a = experiment(10e-3, 5e-3, times)

    fig, ax = plt.subplots(3)

    def plot_ax(i, data):

        ax[i].plot(np.linspace(1,n+1,n), [a.T1]*n)
        ax[i].plot(data)

        ax[i].set_ylabel('T1')
        ax[i].set_xlabel('experiment')
    
    a.loop_1()
    plot_ax(0, a.T1_estimates)
    a.loop_2(n)
    plot_ax(1, a.T1_estimates)
    a.loop_3(n)
    plot_ax(2, a.T1_estimates)
    
    ax[0].set_title('approach A')
    ax[0].legend(['true T1','T1 estimate'])
        
    ax[1].set_title('approach B')
    
    ax[1].legend(['true T1','T1 estimate'])
    
    ax[2].set_title('approach C')
    
    ax[2].legend(['true T1','T1 estimate'])
    

    fig.suptitle('Comparison of 3 different approaches')
    fig.tight_layout()
    fig.show()
    plt.show()
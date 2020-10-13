

import qutip as q
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation as ani
import numpy as np
import matplotlib as mp


class anim_pars():
    """
    Parent class for anim_te and anim_bloch

    it has all the necessary parameters and methods for anim_bloch and ani_bloch_te subclasses.

    Parameters:

        state: qubit state. has to be two dimentional so that it can be represented on a Bloch sphere. normalisation is not required, but recommended

        time: array of numbers that will be a time parameter in the time evolution of the qubit it must be a list type

        name: string with a file extension (e.g. movie.mp4). 

        fps: optional parameter responsible for frames per second. if none given, default is set to 50.
    
    Methods:

        animate_bloch(): creates the animation and saves it in a working directory.


    """
    def __init__(self, state, time, name, fps = 50):
        self.state = state
        self.time = time
        self.name = name
        self.fps = fps

        self.res = []        
        self.figr, self.ax = plt.subplots()
        self.ax = Axes3D(fig = self.figr)#, animated = True)
        self.sphere=q.Bloch(axes=self.ax, fig = self.figr)
        self.sphere.render(fig=self.figr, axes=self.ax)
        self.sphere.add_states(self.state)
        self.txt_pos = (q.sigmaz() - q.sigmay())*0.49
    

    def __animate(self, i):
        pos = self.time.index(i)
        self.sphere.clear()
        self.sphere.add_states(self.res[pos])
        self.sphere.add_annotation(self.txt_pos, "{:.2e} s".format(i))
        self.sphere.make_sphere()
        return self.ax

    def animate_bloch(self):
        anim = ani.FuncAnimation(self.figr, self.__animate, frames=self.time)

        anim.save(self.name, fps = self.fps)

class anim_bloch_te(anim_pars):
    """
    This class uses qutip and matplotlib to animate a time evolution of a qubit. a Hamiltonian must be given so that it can be solved in time using a Master Equation

    Parameters:

        state: qubit state. has to be two dimentional so that it can be represented on a Bloch sphere. normalisation is not required, but recommended

        H: Hamiltonian for the system.

        time: array of numbers that will be a time parameter in the time evolution of the qubit it must be a list type

        name: string with a file extension (e.g. movie.mp4). 

        fps: optional parameter responsible for frames per second. if none given, default is set to 50.
    
    Methods:

        animate_bloch(): creates the animation and saves it in a working directory.

    """
    def __init__(self, state, time, name, H, fps = 50):
        super().__init__(state, time, name, fps = 50)
        self.H = H
    
        #resolve the time evolution in time:
        self.res = q.mesolve(self.H, self.state, self.time)
        self.res = self.res.states
        
class anim_bloch(anim_pars):

    """
    This class uses qutip and matplotlib to animate a time evolution of a qubit.
    This class does NOT take a Hamiltonian, instead it requires a list of states solved in time.

    Parameters:

        state: qubit state. has to be two dimentional so that it can be represented on a Bloch sphere. normalisation is not required, but recommended

        state_list: list of qubit states solved in time

        time: array of numbers that will be a time parameter in the time evolution of the qubit it must be a list type

        name: string with a file extension (e.g. movie.mp4). 

        fps: optional parameter responsible for frames per second. if none given, default is set to 50.
    
    Methods:

        animate_bloch(): creates the animation and saves it in a working directory.

    """

    def __init__(self, state, time, name, state_list, fps = 50):
        super().__init__(state, time, name, fps = 50)

        #resolve the time evolution in time:
        self.res = state_list



if __name__ == '__main__':
        
    up = q.basis(2,0)
    t = list(np.linspace(0,10,50))
    H = q.sigmax()

    res = q.mesolve(H, up, t)


    up = up*up.dag() 
    
    a = anim_bloch_te(up, t, "vid.mp4", H)

    # a = anim_bloch(up, t, "vid.mp4", res.states)

    a.animate_bloch()

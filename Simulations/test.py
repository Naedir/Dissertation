

import qutip as q
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation as ani
import numpy as np
import matplotlib as mp



class anim_bloch():
    """
    This class uses qutip and matplotlib to animate a time evolution of a qubit.

    Parameters:

        state: qubit state. has to be two dimentional so that it can be represented on a Bloch sphere. normalisation is not required, but recommended

        H: Hamiltonian for the system.

        time: array of numbers that will be a time parameter in the time evolution of the qubit

        name: string with a file extension (e.g. movie.mp4). 

        fps: optional parameter responsible for frames per second. if none given, default is set to 50.
    
    Methods:

        animate_bloch(): creates the animation and saves it in a working directory.

    """
    def __init__(self, state, H, time, name, fps = 50):
        self.state = state
        self.H = H
        self.time = time
        self.name = name
        self.fps = fps
    
        #resolve the time evolution in time:
        self.res = q.mesolve(self.H, self.state, self.time)
        
        self.figr, self.ax = plt.subplots()
        self.ax = Axes3D(fig = self.figr)#, animated = True)
        self.sphere=q.Bloch(axes=self.ax, fig = self.figr)
        self.sphere.render(fig=self.figr, axes=self.ax)
        self.sphere.add_states(self.state)
    

    def __animate(self, i):
        self.sphere.clear()
        self.sphere.add_states(self.res.states[i])
        self.sphere.make_sphere()
        self.figr.canvas.draw()
        self.figr.canvas.flush_events()
        return self.figr

    def animate_bloch(self):
        anim = ani.FuncAnimation(self.figr, self.__animate, frames=len(self.time))

        anim.save(self.name, fps = self.fps)
        plt.ioff()
        plt.close('all')

if __name__ =='__main__':


    # # #   This is just an example

    # # #   in this case a qubit is in a magnetic field B

    # # #   The values of B and time are arbitrary, as this is just an example


    ## create an up bit:

    up = (q.Qobj([[1 + 1j],[-1j]])).unit()

    ## add a magnetic field Hamiltonian:
    t = np.linspace(0,20,200)

    H_constants = 1
    H = H_constants/2 * q.sigmaz()

    res = q.mesolve(H, up, t)

    a = anim_bloch(up, H,t,'video.mp4', fps = 30)
    a.animate_bloch()
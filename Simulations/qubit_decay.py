from anim_bloch import  *

#################################
## Time dependant Hamiltonian (H = H_0 + H_1), where H_1 = f(t) * H_1_const 
## unitary evolution can be solved as follows:
## create a python function that returns H_1 coefficients of the form:
##
##      def H_f(t, args):
##          return f(t)
##
## The overall Hamiltonian must be written as a list as follows:
## 
## H = [H_0, [H_1_const, H_f]]

## This Hamiltonian 'H' can be used in a qutip.mesolve() as normal
#################################


####
# Let the relaxation hamiltonian be H(t) = 1/T1 (|0> + |1>), and H(t) = f(t) * H_2 which as a unitary will give an exponential decay
# Therefore H = H_0 + H_t


import time as tm
t1 = tm.time()

## create an up bit:

qubit = (q.Qobj([[1 + 1j],[-1j]])).unit()


up = q.basis(2,1)
down = q.basis(2,0)
# ud = q.sigmax() + q.sigmay()
ud = q.sigmam()
## add a magnetic field Hamiltonian:
t = list(np.linspace(0,10,100))
print(ud)
B_z = 1
H_2 = 1 * ud
T1 = 0.1

def f(t, args):
    return np.e**(t/T1)

H_0 = B_z * q.sigmaz()
H = [H_0, [H_2, f]]

a = anim_bloch(qubit, H_0 ,t,'video1.mp4', fps = 20)
a.animate_bloch()
print(tm.time() - t1)
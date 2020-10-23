from anim_bloch import  *
import sympy as sp

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

qubit = q.basis(2,0)

def bl(state):
    b = q.Bloch()
    b.add_states(state)
    b.make_sphere()
    b.show()
    plt.show()


up = q.basis(2,1)
down = q.basis(2,0)
# qubit = (up+down).unit()

t = np.linspace(0,20,1000)
T1 = 1

Bz = 0.1
H = Bz * q.sigmaz()
hbar = 1

relax = [None] * (len(t))
relax[0] = qubit
stat = [None] * (len(t)-1)
# ex = -1j * (H * (t[2] - t[1])/hbar)
# unitary = ex.expm()
con = 4
con2 = 0.15
for i in range(len(t)):
    if i >= 1:
        stat[i-1] = q.Qobj([[1, (np.exp(1j*con*t[i]))*float(sp.exp(-con2*t[i]/T1))], [np.exp(-1j*con*t[i])*float(sp.exp(-con2*t[i]/T1)), 1]])*1/2

        # matrix = unitary*(stat[i-1]*stat[i-1].dag())*unitary.dag()
        relax[i] = stat[i-1]

        # relax[i] = unitary*relax[i-1]

print(relax[2])
print(relax[0])
print(relax[len(relax) - 1])
a = anim_bloch(qubit, list(t[1::]), "relaxation3.mp4", relax, fps = 60)
a.animate_bloch()

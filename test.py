import qutip as q
import matplotlib.pyplot as plt



b = q.Bloch()


up = q.basis(2,0)
down = 4*q.basis(2,1)


phi = (up * down.dag())
phi = phi.unit()


b.add_states(up)

b.show()
#plt.show()
a = input('\n')
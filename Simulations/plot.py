import pandas as pd
import matplotlib.pyplot as plt

a = pd.read_csv('data/matrix_a.csv', header=None, usecols = [1])
a = a[1].values.tolist()[1::]

b = pd.read_csv('data/matrix_b.csv', header=None, usecols = [1])
b = b[1].values.tolist()[1::]

c = pd.read_csv('data/matrix_c.csv', header=None, usecols = [1])
c = c[1].values.tolist()[1::]

y = [0.02]*len(c)

plt.plot(y); plt.plot(a); plt.plot(b); plt.plot(c, linestyle='--')
plt.title('Comparison of 3 different approaches averaged afrer 1000 data sets')
plt.xlabel('Experiment')
plt.ylabel('T1')
plt.legend(['true T1', 'approach A estimate', 'approach B estimate','approach C estimate'])
plt.show()
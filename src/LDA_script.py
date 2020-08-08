import numpy as np
import matplotlib.pyplot as plt
x=np.loadtxt('./datas.txt')[:,0]
y1=np.loadtxt('./datas.txt')[:,1]
y2=np.loadtxt('./datas.txt')[:,2]
y3=np.loadtxt('./datas.txt')[:,3]
y4=np.loadtxt('./datas.txt')[:,4]
y5=np.loadtxt('./datas.txt')[:,5]
y6=np.loadtxt('./datas.txt')[:,6]
plt.plot(x,y1,'.-b',label='Re=64')
plt.plot(x,y2,'.-g',label='Re=118')
plt.plot(x,y3,'.-r',label='Re=157')
plt.plot(x,y4,'.-',label='Re=245')
plt.plot(x,y5,'.-',label='Re=260')
plt.plot(x,y6,'.-',label='Re=474')
plt.xlabel('Phase')
plt.ylabel('Vis.Pentr.Depth')
plt.legend()
#plt.savefig('./fig1.pdf')
plt.show()



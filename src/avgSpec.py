import matplotlib.pyplot as plt
import numpy as np

for i in range(5680,9980,20):
        data=np.loadtxt(' dns_data_n'+str(i)+'.txt')
        plt.loglog(data[:,1],data[:,0],label='itr='+str(i))

plt.legend()
plt.show()


avg=np.zeros((1,128))
for i in range(5680,9980,20):
        data=np.loadtxt(' dns_data_n'+str(i)+'.txt')
        plt.loglog(data[:,1],data[:,0],label='itr='+str(i))
        avg=avg+np.sum(data[:,0])
        
        
plt.legend()
plt.show()



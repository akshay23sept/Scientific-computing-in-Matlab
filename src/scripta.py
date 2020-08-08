import numpy as np
import matplotlib.pyplot as plt
x=np.loadtxt('./pivlda.dat')[:-6,0]
y1=np.loadtxt('./pivlda.dat')[:-6,1]
y2=np.loadtxt('./pivlda.dat')[:-6,2]
y3=np.loadtxt('./pivlda.dat')[:,3]
y4=np.loadtxt('./pivlda.dat')[:,4]
y5=np.loadtxt('./pivlda.dat')[:,5]
y6=np.loadtxt('./pivlda.dat')[:,6]
plt.plot(x,y1,'.-',label='Re=245(LDA)')
plt.plot(y3,y4,'.-',label='Re=240(PIV)')
plt.plot(x,y2,'.-',label='Re=474(LDA)')
plt.plot(y5,y6,'.-',label='Re=460(PIV)')
plt.xlabel('Phase (in pi radians)')
plt.ylabel('Vis.Pentr.Depth')
plt.legend()
plt.show()

####################
import numpy as np
import matplotlib.pyplot as plt
x=np.loadtxt('./pivfinal.txt')[:,0]
y1=np.loadtxt('./pivfinal.dat')[:,1]
y2=np.loadtxt('./pivfinal.dat')[:,2]
y3=np.loadtxt('./pivfinal.dat')[:,3]
y4=np.loadtxt('./pivfinal.dat')[:,4]
y5=np.loadtxt('./pivfinal.dat')[:,5]
y6=np.loadtxt('./pivfinal.dat')[:,6]
y7=np.loadtxt('./pivfinal.dat')[:,7]
y8=np.loadtxt('./pivfinal.dat')[:,8]
y9=np.loadtxt('./pivfinal.dat')[:,9]
y10=np.loadtxt('./pivfinal.dat')[:,10]
y11=np.loadtxt('./pivfinal.dat')[:,11]
y12=np.loadtxt('./pivfinal.dat')[:,12]
y13=np.loadtxt('./pivfinal.dat')[:,13]
y14=np.loadtxt('./pivfinal.dat')[:,14]
y15=np.loadtxt('./pivfinal.dat')[:,15]
y16=np.loadtxt('./pivfinal.dat')[:,16]
y17=np.loadtxt('./pivfinal.dat')[:,17]
plt.plot(x,y1,'.-',label='Re=205')
plt.plot(y2,y3,'.-',label='Re=240')
plt.plot(y4,y5,'.-',label='Re=272')
plt.plot(y6,y7,'.-',label='Re=302')
plt.plot(y8,y9,'.-',label='Re=336')
plt.plot(y10,y11,'.-',label='Re=375')
plt.plot(y12,y13,'.-',label='Re=418')
plt.plot(y14,y15,'.-',label='Re=445')
plt.plot(y16,y17,'.-',label='Re=466')
plt.xlabel('Phase (in pi radians)')
plt.ylabel('Vis.Pentr.Depth')
plt.legend()
plt.show()
################################



import numpy as np
import matplotlib.pyplot as plt
x=np.loadtxt('./clarcy.dat')[:,0]
y1=np.loadtxt('./clarcy.dat')[:,1]
plt.show()

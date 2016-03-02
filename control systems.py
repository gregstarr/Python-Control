import control
import matplotlib.pyplot as plt
import numpy
import time

def analyze(openloop):
    freq=numpy.arange(1,90000)
    mag,phase,omega=control.bode(openloop,omega=freq,dB=True,deg=True,Plot=False)
    pm=phase[numpy.argmin(numpy.abs(mag))]+180
    gm=-mag[numpy.argmin(numpy.abs(phase+180))]

    closedloop=control.feedback(openloop)
    t=numpy.arange(100)
    t,y=control.step_response(closedloop,t)
    overshoot=numpy.amax(y)-1
    settime=len(y)-numpy.argmax(numpy.abs(y[::-1]-1)>.1)

    return pm,gm,overshoot,settime
    

T = 1/30000     #sampling period
kp=10**8*T*T*.5 #plant gain
leadpole=-.5    #lead compensator pole
plant=control.tf([kp,kp],[1,-2,1],T)

cnt=control.tf([1,-1],[1,leadpole],T)
ol=control.series(plant,cnt)
poles,k=control.rlocus(ol,Plot=False)
p1 = numpy.real(poles[:,0])**2+numpy.imag(poles[:,0])**2<1
p2 = numpy.real(poles[:,1])**2+numpy.imag(poles[:,1])**2<1
stablek = numpy.logical_and(p1,p2)

data=[]
for i in numpy.nonzero(stablek)[0]:
    t1=time.time()
    cnt=control.tf([k[i],-k[i]],[1,leadpole],T)
    ol=control.series(plant,cnt)
    a=analyze(ol)
    data.append(a)
    print(a)
    print(time.time()-t1)

"""
plt.figure(1)
plt.title('Root Locus Plot')
plt.plot(numpy.real(poles[:,0]),numpy.imag(poles[:,0]))
plt.plot(numpy.real(poles[:,1]),numpy.imag(poles[:,1]))
plt.plot(numpy.real(poles[:,2]),numpy.imag(poles[:,2]))
cos=[numpy.cos(i) for i in numpy.arange(0,2*numpy.pi,.02*numpy.pi)]
sin=[numpy.sin(j) for j in numpy.arange(0,2*numpy.pi,.02*numpy.pi)]
plt.plot(cos,sin)
plt.grid(True)
plt.axis([-2,2,-2,2])
"""



"""
cl=control.feedback(ol)

plt.figure(2)
mag,phase,omega=control.bode(ol,omega=freq,dB=True,deg=True)

t=numpy.arange(100)
t,y=control.step_response(cl,T=t)
plt.figure(3)
plt.plot(t,y[:len(t)])
plt.show()

char=numpy.array([ol.minreal().den[0][0][0],ol.minreal().den[0][0][1]+ol.minreal().num[0][0][0],ol.minreal().den[0][0][2]+ol.minreal().num[0][0][1]])

"""

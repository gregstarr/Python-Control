import control
import matplotlib.pyplot as plt
import numpy
import time

#this function does long division up to a specified number of times
def divide(num,den,limit): 
    answer=[]
    i=0
    n=(len(den)-len(num))*[0]+list(num)
    d=list(den)

    while(i<limit):
        co=n[0]/d[0]
        answer.append(co)
        n=[n[j]-co*d[j] for j in range(len(n))]
        n.pop(0)
        n.append(0)
        i+=1

    return numpy.array(answer)

#this function takes an open loop system and returns phase margin, gain margin,
#phase crossover frequency, gain crossover frequency, overshoot on step response,
#and the 10% settling time
def analyze(openloop):
    freq=numpy.arange(1,10000)
    mag,phase,omega=control.bode(openloop,omega=freq,dB=True,deg=True,Hz=True,Plot=False)
    gcf=numpy.argmin(numpy.abs(mag))
    pm=phase[gcf]+180
    pcf=numpy.argmin(numpy.abs(phase+180))
    gm=-mag[pcf]

    closedloop=control.feedback(openloop)
    n=numpy.convolve(closedloop.num[0][0],[1,0])
    d=numpy.convolve(closedloop.den[0][0],[1,-1])
    
    y=divide(n,d,10000)
    overshoot=numpy.amax(y)-1
    settime=(len(y)-numpy.argmax(numpy.abs(y[::-1]-1)>.1))*T

    return pm,gm,pcf,gcf,overshoot,settime
    

T = 1/30000     #sampling period
kp=10**8*T*T*.5 #plant gain
plant=control.tf([kp,kp],[1,-2,1],T)

poletests=numpy.arange(-.1,-1,-.1)
data=numpy.zeros((9,50,8))

for p,pole in enumerate(poletests):
    data[p,:,0]=pole
    t1=time.time()
    print("pole: {}".format(pole))
    cnt=control.tf([1,-1],[1,pole],T)
    ol=control.series(plant,cnt)
    poles,k=control.rlocus(ol,Plot=False)
    p1 = numpy.real(poles[:,0])**2+numpy.imag(poles[:,0])**2<1
    p2 = numpy.real(poles[:,1])**2+numpy.imag(poles[:,1])**2<1
    stablek = numpy.logical_and(p1,p2)

    for i in numpy.nonzero(stablek)[0]:
        cnt=control.tf([k[i],-k[i]],[1,pole],T)
        ol=control.series(plant,cnt)
        a=analyze(ol)
        data[p,i,1]=k[i]
        data[p,i,2:]=a
        
    print("{} seconds".format(time.time()-t1))

st=.01
os=.2
pm=40
gm=5
gc=500*2*numpy.pi

stMask=data[:,:,7]<st
osMask=data[:,:,6]<os
pmMask=data[:,:,2]>pm
gmMask=data[:,:,3]>gm
gcMask=data[:,:,5]>gc

mask=numpy.logical_and.reduce((stMask,osMask,pmMask,gmMask,gcMask))

#i picked the fifth combination in data[mask]
kc=data[mask][9][1]
pole=data[mask][9][0]
cnt=control.tf([kc,-kc],[1,pole],T)
openloop=control.series(cnt,plant)
closedloop=control.feedback(openloop)


plt.figure(1)
freq=numpy.arange(1,10000)
mag,phase,omega=control.bode(openloop,omega=freq,dB=True,deg=True)


t,y=control.step_response(closedloop)
plt.figure(2)
plt.title("Step Response")
plt.xlabel("Sample Periods")
plt.plot(y[:100])
plt.grid(True)
plt.axis([-1,100,0,1.5])
plt.show()



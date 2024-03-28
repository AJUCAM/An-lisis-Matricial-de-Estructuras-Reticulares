import sympy as sy
import numpy as np
import matplotlib.pyplot as plt



E,I0,L,x,A0,r,xA,xB,xC,h0A,h0B,h0C,xi,xiA,xiB,xiC,I = sy.symbols('E,I_0,L,x,A_0,r,x_A,x_B,x_C,h_0A,h_0B,h_0C,xi,xi_A,xi_B,xi_C,I')

E = sy.Rational(2.2e7,1)

LA = int(3)
LB = int(4)
LC = int(3)


qvA = -60 + xA*40/3 #entre 0  y 3m
qvB = -20 - xB*40/4 #entre 0  y 4m
qvC = -60 + xC*50/3 #entre 0  y 3m

bA = sy.Rational(9,20)
bB = sy.Rational(9,20)
bC = sy.Rational(9,20)

h0A = sy.Rational(11,20)
h0B = sy.Rational(2,5) 
h0C = sy.Rational(3,5) 

rA = sy.Rational(-3,20)
rB = sy.Rational(1,5)
rC = sy.Rational(-3,20)

fhA = h0A + rA*xA/LA
fhB = h0B + rB*xB/LB
fhC = h0C + rC*xC/LC

IA = sy.Rational(1,12)*bA*h0A**3
IB = sy.Rational(1,12)*bA*h0B**3
IC = sy.Rational(1,12)*bA*h0C**3

IxA = sy.Rational(1,12)*bA*fhA**3
IxB = sy.Rational(1,12)*bB*fhB**3
IxC = sy.Rational(1,12)*bC*fhC**3


#%% ELEMENTO A
I=IxA
L=LA
x=xA
xi=xiA


I3=sy.integrate(x/(E*I),(x,0,x))
I4=sy.integrate(1/(E*I),(x,0,x))
I1=sy.integrate(I3,(x,0,x))
I2=sy.integrate(I4,(x,0,x))

C2=sy.simplify(-I1+x*(I2+I3)-x**2*I4)
C3=sy.simplify(x*I4-I3)
C5=sy.simplify(-I1.subs({x:L})+I1+x*(I2.subs({x:L})-I2-(L-x)*I4)+(L-x)*I3)
C6=sy.simplify(x*(I4.subs({x:L})-I4)-(I3.subs({x:L})-I3))

k22A=-I4.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k55A=k22A
k25A=-k22A
k52A=-k22A

k23A=-I3.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k32A=k23A
k35A=-k23A
k53A=-k23A

k26A=-I2.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k62A=k26A
k56A=-k26A
k65A=-k26A

k33A=(I1.subs({x:L})-I3.subs({x:L})*L)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k36A=-I1.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k63A=k36A
k66A=(I1.subs({x:L})-I2.subs({x:L})*L)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))

kVigaA = sy.zeros(4,4)	

kVigaA[0,0]=k22A
kVigaA[0,1]=k23A
kVigaA[0,2]=k25A
kVigaA[0,3]=k26A

kVigaA[1,0]=k32A
kVigaA[1,1]=k33A
kVigaA[1,2]=k35A
kVigaA[1,3]=k36A

kVigaA[2,0]=k52A
kVigaA[2,1]=k53A
kVigaA[2,2]=k55A
kVigaA[2,3]=k56A

kVigaA[3,0]=k62A
kVigaA[3,1]=k63A
kVigaA[3,2]=k65A
kVigaA[3,3]=k66A

kVigaA=sy.N(kVigaA)


psi2A=sy.simplify((-I4.subs({x:L})*I1+I3.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))+1)
psi3A=sy.simplify(((I2.subs({x:L})-I4.subs({x:L})*L)*I1+(-I1.subs({x:L})+I3.subs({x:L})*L)*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))+x)
psi5A=sy.simplify((I4.subs({x:L})*I1-I3.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L})))
psi6A=sy.simplify((-I2.subs({x:L})*I1+I1.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L})))

Gyy1A=sy.simplify(sy.expand(C2*psi2A.subs({xA:xiA})+C3*psi3A.subs({xA:xiA})))
Gyy2A=sy.simplify(sy.expand(C5*psi5A.subs({xA:xiA})+C6*psi6A.subs({xA:xiA})))


#%% ELEMENTO B

I=IxB
L=LB
x=xB
xi=xiB


I3=sy.integrate(x/(E*I),(x,0,x))
I4=sy.integrate(1/(E*I),(x,0,x))
I1=sy.integrate(I3,(x,0,x))
I2=sy.integrate(I4,(x,0,x))

C2=sy.simplify(-I1+x*(I2+I3)-x**2*I4)
C3=sy.simplify(x*I4-I3)
C5=sy.simplify(-I1.subs({x:L})+I1+x*(I2.subs({x:L})-I2-(L-x)*I4)+(L-x)*I3)
C6=sy.simplify(x*(I4.subs({x:L})-I4)-(I3.subs({x:L})-I3))

k22B=-I4.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k55B=k22B
k25B=-k22B
k52B=-k22B

k23B=-I3.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k32B=k23B
k35B=-k23B
k53B=-k23B

k26B=-I2.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k62B=k26B
k56B=-k26B
k65B=-k26B

k33B=(I1.subs({x:L})-I3.subs({x:L})*L)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k36B=-I1.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k63B=k36B
k66B=(I1.subs({x:L})-I2.subs({x:L})*L)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))

kVigaB = sy.zeros(4,4)

kVigaB[0,0]=k22B
kVigaB[0,1]=k23B
kVigaB[0,2]=k25B
kVigaB[0,3]=k26B

kVigaB[1,0]=k32B
kVigaB[1,1]=k33B
kVigaB[1,2]=k35B
kVigaB[1,3]=k36B

kVigaB[2,0]=k52B
kVigaB[2,1]=k53B
kVigaB[2,2]=k55B
kVigaB[2,3]=k56B

kVigaB[3,0]=k62B
kVigaB[3,1]=k63B
kVigaB[3,2]=k65B
kVigaB[3,3]=k66B

kVigaB=sy.N(kVigaB)


psi2B=sy.simplify((-I4.subs({x:L})*I1+I3.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))+1)
psi3B=sy.simplify(((I2.subs({x:L})-I4.subs({x:L})*L)*I1+(-I1.subs({x:L})+I3.subs({x:L})*L)*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))+x)
psi5B=sy.simplify((I4.subs({x:L})*I1-I3.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L})))
psi6B=sy.simplify((-I2.subs({x:L})*I1+I1.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L})))

Gyy1B=sy.simplify(C2*psi2B.subs({xB:xiB})+C3*psi3B.subs({xB:xiB}))
Gyy2B=sy.simplify(C5*psi5B.subs({xB:xiB})+C6*psi6B.subs({xB:xiB}))





#%% ELEMENTO C

I=IxC
L=LC
x=xC
xi=xiC

I3=sy.integrate(x/(E*I),(x,0,x))
I4=sy.integrate(1/(E*I),(x,0,x))
I1=sy.integrate(I3,(x,0,x))
I2=sy.integrate(I4,(x,0,x))

C2=sy.simplify(-I1+x*(I2+I3)-x**2*I4)
C3=sy.simplify(x*I4-I3)
C5=sy.simplify(-I1.subs({x:L})+I1+x*(I2.subs({x:L})-I2-(L-x)*I4)+(L-x)*I3)
C6=sy.simplify(x*(I4.subs({x:L})-I4)-(I3.subs({x:L})-I3))

k22C=-I4.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k55C=k22C
k25C=-k22C
k52C=-k22C

k23C=-I3.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k32C=k23C
k35C=-k23C
k53C=-k23C

k26C=-I2.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k62C=k26C
k56C=-k26C
k65C=-k26C


k33C=(I1.subs({x:L})-I3.subs({x:L})*L)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k36C=-I1.subs({x:L})/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))
k63C=k36C
k66C=(I1.subs({x:L})-I2.subs({x:L})*L)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))

kVigaC = sy.zeros(4,4)

kVigaC[0,0]=k22C
kVigaC[0,1]=k23C
kVigaC[0,2]=k25C
kVigaC[0,3]=k26C

kVigaC[1,0]=k32C
kVigaC[1,1]=k33C
kVigaC[1,2]=k35C
kVigaC[1,3]=k36C

kVigaC[2,0]=k52C
kVigaC[2,1]=k53C
kVigaC[2,2]=k55C
kVigaC[2,3]=k56C

kVigaC[3,0]=k62C
kVigaC[3,1]=k63C
kVigaC[3,2]=k65C
kVigaC[3,3]=k66C

kVigaC=sy.N(kVigaC)


psi2C=sy.simplify((-I4.subs({x:L})*I1+I3.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))+1)
psi3C=sy.simplify(((I2.subs({x:L})-I4.subs({x:L})*L)*I1+(-I1.subs({x:L})+I3.subs({x:L})*L)*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L}))+x)
psi5C=sy.simplify((I4.subs({x:L})*I1-I3.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L})))
psi6C=sy.simplify((-I2.subs({x:L})*I1+I1.subs({x:L})*I2)/(I1.subs({x:L})*I4.subs({x:L})-I2.subs({x:L})*I3.subs({x:L})))

Gyy1C=sy.simplify(sy.expand(C2*psi2C.subs({xC:xiC})+C3*psi3C.subs({xC:xiC})))
Gyy2C=sy.simplify(sy.expand(C5*psi5C.subs({xC:xiC})+C6*psi6C.subs({xC:xiC})))


#%% DESARROLLO
fEmpA = sy.zeros(4,1)

fEmpA[0,0]=-sy.re(sy.N(sy.integrate(psi2A*qvA,(xA,0,LA))))
fEmpA[1,0]=-sy.re(sy.N(sy.integrate(psi3A*qvA,(xA,0,LA))))
fEmpA[2,0]=-sy.re(sy.N(sy.integrate(psi5A*qvA,(xA,0,LA))))
fEmpA[3,0]=-sy.re(sy.N(sy.integrate(psi6A*qvA,(xA,0,LA))))

fEmpB = sy.zeros(4,1)

fEmpB[0,0]=-sy.re(sy.N(sy.integrate(psi2B*qvB,(xB,0,LB))))
fEmpB[1,0]=-sy.re(sy.N(sy.integrate(psi3B*qvB,(xB,0,LB))))
fEmpB[2,0]=-sy.re(sy.N(sy.integrate(psi5B*qvB,(xB,0,LB))))
fEmpB[3,0]=-sy.re(sy.N(sy.integrate(psi6B*qvB,(xB,0,LB))))

fEmpC = sy.zeros(4,1)

fEmpC[0,0]=-sy.re(sy.N(sy.integrate(psi2C*qvC,(xC,0,LC))))
fEmpC[1,0]=-sy.re(sy.N(sy.integrate(psi3C*qvC,(xC,0,LC))))
fEmpC[2,0]=-sy.re(sy.N(sy.integrate(psi5C*qvC,(xC,0,LC))))
fEmpC[3,0]=-sy.re(sy.N(sy.integrate(psi6C*qvC,(xC,0,LC))))



SFYPA=fEmpA[0,0]+fEmpA[2,0]+sy.integrate(qvA,(xA,0,LA))
SFYPB=fEmpB[0,0]+fEmpB[2,0]+sy.integrate(qvB,(xB,0,LB))
SFYPC=fEmpC[0,0]+fEmpC[2,0]+sy.integrate(qvC,(xC,0,LC))

#%% SOL RIGIDEZ
kRes = sy.zeros(6,6)


kRes[0,0] = kVigaA[1,1]
kRes[0,1] = kVigaA[1,3]
kRes[0,2] = 0
kRes[0,3] = 0
kRes[0,4] = 0
kRes[0,5] = 0

kRes[1,0] = kVigaA[3,1]
kRes[1,1] = kVigaA[3,3] + kVigaB[1,1]
kRes[1,2] = kVigaB[1,2]
kRes[1,3] = kVigaB[1,3]
kRes[1,4] = 0
kRes[1,5] = 0

kRes[2,0] = 0
kRes[2,1] = kVigaB[2,1]
kRes[2,2] = kVigaB[2,2] + kVigaC[0,0]
kRes[2,3] = kVigaB[2,3]
kRes[2,4] = kVigaC[0,1]
kRes[2,5] = kVigaC[0,3]

kRes[3,0] = 0
kRes[3,1] = kVigaB[3,1]
kRes[3,2] = kVigaB[3,2]
kRes[3,3] = kVigaB[3,3]
kRes[3,4] = 0
kRes[3,5] = 0

kRes[4,0] = 0
kRes[4,1] = 0
kRes[4,2] = kVigaC[1,0]
kRes[4,3] = 0
kRes[4,4] = kVigaC[1,1]
kRes[4,5] = kVigaC[1,3]

kRes[5,0] = 0
kRes[5,1] = 0
kRes[5,2] = kVigaC[3,0]
kRes[5,3] = 0
kRes[5,4] = kVigaC[3,1]
kRes[5,5] = kVigaC[3,3]


fEmpRes = sy.zeros(6,1)

fEmpRes[0,0] = fEmpA[1,0]
fEmpRes[1,0] = fEmpA[3,0] + fEmpB[1,0]
fEmpRes[2,0] = fEmpB[2,0] + fEmpC[0,0]
fEmpRes[3,0] = fEmpB[3,0]
fEmpRes[4,0] = fEmpC[1,0]
fEmpRes[5,0] = fEmpC[3,0]

fNodRes = sy.zeros(6,1)
fNodRes[0,0] = 0
fNodRes[1,0] = 0
fNodRes[2,0] = 0
fNodRes[3,0] = 0
fNodRes[4,0] = 0
fNodRes[5,0] = 0


Sol = sy.linsolve([kRes,fNodRes-fEmpRes])
v1 = 0
theta1 = Sol.args[0][0]
v2 = 0 
theta2 = Sol.args[0][1]
v3 = Sol.args[0][2]
theta3b = Sol.args[0][3]
theta3c = Sol.args[0][4]
v4 = 0
theta4 = Sol.args[0][5]

#%% DESPLAZAMIENTOS NODALES

SolMat = sy.zeros(6,1)
SolMat[0,0] = theta1
SolMat[1,0] = theta2
SolMat[2,0] = v3
SolMat[3,0] = theta3b
SolMat[4,0] = theta3c
SolMat[5,0] = theta4

Mat = sy.zeros(3,6)

Mat[0,0] = kVigaA[0,1]
Mat[0,1] = kVigaA[0,3]
Mat[0,2] = 0
Mat[0,3] = 0
Mat[0,4] = 0
Mat[0,5] = 0

Mat[1,0] = kVigaA[2,1]
Mat[1,1] = kVigaA[2,3] + kVigaB[0,1]
Mat[1,2] = kVigaB[0,2]
Mat[1,3] = kVigaB[0,3]
Mat[1,4] = 0
Mat[1,5] = 0

Mat[2,0] = 0
Mat[2,1] = 0
Mat[2,2] = kVigaC[2,0]
Mat[2,3] = 0
Mat[2,4] = kVigaC[2,1]
Mat[2,5] = kVigaC[2,3]

Matf = sy.zeros(3,1)

Matf[0,0] = fEmpA[0,0]
Matf[1,0] = fEmpA[2,0] + fEmpB[0,0]
Matf[2,0] = fEmpC[2,0]

MatR = Mat*SolMat + Matf


#%% CAMPO A
vAh=psi2A*v1 + psi3A*theta1 + psi5A*v2 + psi6A*theta2
vAf=sy.integrate(sy.nsimplify(Gyy2A*qvA.subs({xA:xiA})),(xiA,0,xA))+sy.integrate(sy.nsimplify(Gyy1A*qvA.subs({xA:xiA})),(xiA,xA,LA))
vA = sy.expand(sy.N(vAh+vAf))

#%% CAMPO B
vBh=psi2B*v2 + psi3B*theta2 + psi5B*v3 + psi6B*theta3b
vBf=sy.integrate(sy.nsimplify(Gyy2B*qvB.subs({xB:xiB})),(xiB,0,xB))+sy.integrate(sy.nsimplify(Gyy1B*qvB.subs({xB:xiB})),(xiB,xB,LB))
vB = sy.expand(sy.N(vBh+vBf))
#%% CAMPO C
vCh=psi2C*v3 + psi3C*theta3c + psi5C*v4 + psi6C*theta4
vCf=sy.integrate(sy.nsimplify(Gyy2C*qvC.subs({xC:xiC})),(xiC,0,xC))+sy.integrate(sy.nsimplify(Gyy1C*qvC.subs({xC:xiC})),(xiC,xC,LC))
vC = sy.expand(sy.N(vCh+vCf))
#%%
VA=sy.re(sy.N(-sy.diff(E*IxA*sy.diff(vA,xA,2),xA,1)))
MA=sy.re(sy.N(+E*IxA*sy.diff(vA,xA,2)))

VB=sy.re(sy.N(-sy.diff(E*IxB*sy.diff(vB,xB,2),xB,1)))
MB=sy.re(sy.N(+E*IxB*sy.diff(vB,xB,2)))

VC=sy.re(sy.N(-sy.diff(E*IxC*sy.diff(vC,xC,2),xC,1)))
MC=sy.re(sy.N(+E*IxC*sy.diff(vC,xC,2)))


FY1 = sy.N(-VA.subs({xA:0}))
FY2 = sy.N(VA.subs({xA:LA})-VB.subs({xB:0}))
FY4 = sy.N(VC.subs({xC:LC}))

CA1 = sy.re(sy.N(FY1 + FY2 + FY4 + sy.integrate(qvA,(xA,0,LA)) + sy.integrate(qvB,(xB,0,LB)) + sy.integrate(qvC,(xC,0,LC))))


#%% CAMPO DE DESPLAZAMIENTOS
N=100

xAp=np.linspace(0,LA,N)
xBp=np.linspace(0,LB,N)
xCp=np.linspace(0,LC,N)

vAr=np.zeros([N])
vBr=np.zeros([N])
vCr=np.zeros([N])


for i in range(N):
    vAr[i]=sy.re(vA.subs({xA:xAp[i]}))
    vBr[i]=sy.re(vB.subs({xB:xBp[i]}))
    vCr[i]=sy.re(vC.subs({xC:xCp[i]}))
    
plt.figure(3,figsize=(18,6))
plt.plot(xAp,vAr,color='g',linewidth = 4.0)
plt.plot(LA+xBp,vBr,color='r',linewidth = 4.0)
plt.plot(LA+LB+xCp,vCr,color='b',linewidth = 4.0)
plt.plot([0,LA+LB+LC],[0,0],color='k')
plt.plot([0,3,7,10],[0,0,0,0],'ko')
plt.xlabel(r'$x$ [m]',fontsize=16)
plt.ylabel(r'$Desplazamiento$ [m]',fontsize=16)
plt.tick_params(labelsize=16)
plt.legend(["Elemento A","Elemento B","Elemento C"], ncol = 3)
plt.grid('on')

#%% CAMPO DE CORTANTE
N=100

xAp=np.linspace(0,LA,N)
xBp=np.linspace(0,LB,N)
xCp=np.linspace(0,LC,N)

vAv=np.zeros([N])
vBv=np.zeros([N])
vCv=np.zeros([N])

vAm=np.zeros([N])
vBm=np.zeros([N])
vCm=np.zeros([N])

for i in range(N):
    vAv[i]=sy.re(VA.subs({xA:xAp[i]}))
    vBv[i]=sy.re(VB.subs({xB:xBp[i]}))
    vCv[i]=sy.re(VC.subs({xC:xCp[i]}))
    
    vAm[i]=sy.re(MA.subs({xA:xAp[i]}))
    vBm[i]=sy.re(sy.N(MB.subs({xB:xBp[i]})))
    vCm[i]=sy.re(MC.subs({xC:xCp[i]}))
    
plt.figure(3,figsize=(18,6))
plt.plot(xAp,vAv,color='g',linewidth = 4.0)
plt.plot(LA+xBp,vBv,color='r',linewidth = 4.0)
plt.plot(LA+LB+xCp,vCv,color='b',linewidth = 4.0)
plt.plot([0,LA+LB+LC],[0,0],color='k')
plt.plot([0,3,7,10],[0,0,0,0],'ko')
plt.xlabel(r'$x$ [m]',fontsize=16)
plt.ylabel(r'$V$ [KN]',fontsize=16)
plt.legend(["Elemento A","Elemento B","Elemento C"], ncol = 3)
plt.tick_params(labelsize=16)
plt.grid('on')


#%% CAMPO DE MOMENTOS
N=100

xAp=np.linspace(0,LA,N)
xBp=np.linspace(0,LB,N)
xCp=np.linspace(0,LC,N)

vAv=np.zeros([N])
vBv=np.zeros([N])
vCv=np.zeros([N])

vAm=np.zeros([N])
vBm=np.zeros([N])
vCm=np.zeros([N])

for i in range(N):
    vAv[i]=sy.re(VA.subs({xA:xAp[i]}))
    vBv[i]=sy.re(VB.subs({xB:xBp[i]}))
    vCv[i]=sy.re(VC.subs({xC:xCp[i]}))
    
    vAm[i]=sy.re(MA.subs({xA:xAp[i]}))
    vBm[i]=sy.re(sy.N(MB.subs({xB:xBp[i]})))
    vCm[i]=sy.re(MC.subs({xC:xCp[i]}))


plt.figure(3,figsize=(18,6))
plt.plot(xAp,vAm,color='g',linewidth = 4.0)
plt.plot(LA+xBp,vBm,color='r',linewidth = 4.0)
plt.plot(LA+LB+xCp,vCm,color='b',linewidth = 4.0)
plt.plot([0,LA+LB+LC],[0,0],color='k')
plt.plot([0,3,7,10],[0,0,0,0],'ko')
plt.xlabel(r'$x$ [m]',fontsize=16)
plt.ylabel(r'$M$ [KN-m]',fontsize=16)
plt.legend(["Elemento A","Elemento B","Elemento C"], ncol = 3)
plt.tick_params(labelsize=16)
plt.gca().invert_yaxis()
plt.grid('on')

#%%DESPLAZAMIENTOS
ValvA = []
for i in range (0,11):
    vAp = sy.re(sy.N(-vA.subs({xA:(i*LA/10)})))
    ValvA.append(vAp)
vVA = np.array(ValvA).reshape(-1,1)
    


ValvB = []
for i in range (0,11):
    vBp = sy.re(sy.N(-vB.subs({xB:(i*LB/10)})))
    ValvB.append(vBp)
vVB = np.array(ValvB).reshape(-1,1)
    
                
ValvC = []
for i in range (0,11):
    vCp = sy.re(sy.N(-vC.subs({xC:(i*LC/10)})))
    ValvC.append(vCp)
vVC = np.array(ValvC).reshape(-1,1)

vVA




#%% CORTANTES
ValVA = []
for i in range (0,11):
    VAp = sy.re(sy.N(-VA.subs({xA:(i*LA/10)})))
    ValVA.append(VAp)
VVA = np.array(ValVA).reshape(-1,1)
    


ValVB = []
for i in range (0,11):
    VBp = sy.re(sy.N(-VB.subs({xB:(i*LB/10)})))
    ValVB.append(VBp)
VVB = np.array(ValVB).reshape(-1,1)
    
                
ValVC = []
for i in range (0,11):
    VCp = sy.re(sy.N(-VC.subs({xC:(i*LC/10)})))
    ValVC.append(VCp)
VVC = np.array(ValVC).reshape(-1,1)

#%%MOMENTOS
ValMA = []
for i in range (0,11):
    MAp = sy.re(sy.N(-MA.subs({xA:(i*LA/10)})))
    ValMA.append(MAp)
MMA = np.array(ValMA).reshape(-1,1)
    


ValMB = []
for i in range (0,11):
    MBp = sy.re(sy.N(-MB.subs({xB:(i*LB/10)})))
    ValMB.append(MBp)
MMB = np.array(ValMB).reshape(-1,1)
    
                
ValMC = []
for i in range (0,11):
    MCp = sy.re(sy.N(-MC.subs({xC:(i*LC/10)})))
    ValMC.append(MCp)
MMC = np.array(ValMC).reshape(-1,1)


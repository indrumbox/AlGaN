import math

h1=0.02
h0=0.01
h=h0/1


x=[0,0,0,0,0,0,0,0,0,0]
y=[0,0,0,0,0,0,0,0,0,0]

def C1Koeff(f):
    g = -0.4570245289 * f + 1.151745336
    return g

def A1Koeff(f):
    return 0.465

def A2Koeff(f):
    g = -9.189137089 *f + 7.669893919
    return g
def C2Koeff(f):
    g = 27.26910615 * f - 23.03122697
    return g

def A3Koeff(f):
    g = 0.465/2
    return g

def C3Koeff(f):
    g = 4.141167261 - 0.0121604643 * f - 1.450673214 * pow(f,2) + 3.291528196 * pow(f,3) - 4.355129034 * pow(f,4)
    return g
    
    
def Peresechenie(A1,C1, A2,C2):
    det0=-A2+A1
    det1=-C1*A2+C2*A1
    det2=C2-C1
    x0=det2/det0
    y0=det1/det0
    return x0,y0

def Lagrvypukl(A1,C1,A2,C2,z):
    x0,y0=Peresechenie(A1,C1,A2,C2)   
    x[5]=x0
    x[1]=x0-h0*4
    x[2]=x0-h0*3
    x[3]=x0-h0*2
    x[4]=x0-h0
    x[6]=x0+h0
    x[7]=x0+h0*2
    x[8]=x0+h0*3
    x[9]=x0+4*h0

    y[1]=A1*x[1]+C1
    y[9]=A2*x[9]+C2
    y[5]=0.2*y[1]+0.8*y[9]
    y[3]=0.3*y[1]+0.7*y[5]
    y[2]=0.4*y[1]+0.6*y[3]
    y[4]=0.6*y[5]+0.4*y[3]
    y[7]=0.4*y[5]+0.6*y[9]
    y[8]=0.4*y[7]+0.6*y[9]
    y[6]=0.4*y[5]+0.6*y[7]

    j=10
    while j<10:
        print(y[j])
        j=j+1

    j0=1
    L=0
    l=[1,1,1,1,1,1,1,1,1,1,1]
    while j0<10:
        i0=0
        while i0<10:
            if i0!=j0:
                l[j0]=l[j0]*(z-x[i0])/(x[j0]-x[i0])
            i0=i0+1
                
        L=y[j0]*l[j0]+L
        j0=j0+1
    if 1:
        return L


def Lagrvognut(A1,C1,A2,C2,z):
    x0,y0=Peresechenie(A1,C1,A2,C2)   
    
    x[5]=x0
    x[1]=x0-h1*4
    x[2]=x0-h1*3
    x[3]=x0-h1*2
    x[4]=x0-h1
    x[6]=x0+h1
    x[7]=x0+h1*2
    x[8]=x0+h1*3
    x[9]=x0+4*h1

    y[1]=A1*x[1]+C1
    y[9]=A2*x[9]+C2
    y[5]=0.8*y[1]+0.2*y[9]
    y[3]=0.7*y[1]+0.3*y[5]
    y[2]=0.6*y[1]+0.4*y[3]
    y[4]=0.4*y[5]+0.6*y[3]
    y[7]=0.6*y[5]+0.4*y[9]
    y[8]=0.6*y[7]+0.4*y[9]
    y[6]=0.6*y[5]+0.4*y[7]

    j=10
    while j<10:
        print(y[j])
        j=j+1

    j0=1
    L=0
    l=[1,1,1,1,1,1,1,1,1,1,1]
    while j0<10:
        i0=0
        while i0<10:
            if i0!=j0:
                l[j0]=l[j0]*(z-x[i0])/(x[j0]-x[i0])
            i0=i0+1
                
        L=y[j0]*l[j0]+L
        j0=j0+1
    return L
          

def Linear(A,C,x):
    return A*x+C

def Equation(x1,y1,x2,y2): 
    det0=x1*y2-x2*y1
    det1=-y2+y1
    det2=-x1+x2
    A=det1/det0
    B=det2/det0
    A=-A/B
    C=-1/B
    return (A, C)


def Eg(x1):
    return (3.42*(1-x1)+6.13*x1-1*x1*(1-x1))

#print(Eg(0.11),Eg(0.2), Eg(0.5))

A1011,C1011=Equation(2.9245,2.602,3.115,2.6989)
A2011,C2011=Equation(3.48,3.3,3.5849,4)


A1020,C1020=Equation(2.679,2.3,3.0566,2.4771)
A2020,C2020=Equation(3.6154,3.3,3.7358,4)


A1050,C1050=Equation(3.3773,2.4771,2,1.845)
A2050,C2050=Equation(4.1698,3.4771,4.3396,4)

def Pe4at(f):
    A1=A1Koeff(f)
    A2=A2Koeff(f)
    C1=C1Koeff(f)
    C2=C2Koeff(f)
    A3=A3Koeff(f)
    C3=C3Koeff(f)
    
    x0,y0=Peresechenie(A1Koeff(f),C1Koeff(f),A2Koeff(f),C2Koeff(f))
    x1,y1=Peresechenie(A2Koeff(f),C2Koeff(f),A3Koeff(f),C3Koeff(f))
    w=3
    
    while Linear(A3,C3,w)<5.04:
        if w<=(x0-4*h1):
            print (w,pow(10,Linear(A1,C1,w)))
        if (w>(x0-4*h1))&(w<(x0+4*h1)):
            print (w,pow(10,Lagrvognut(A1,C1,A2,C2,w)))
        if (w>(x0+4*h1))&(w<(x1-4*h0)):
            print (w,pow(10,Linear(A2,C2,w)))
        if (w>(x1-4*h0))&(w<(x1+4*h0)):
            print (w,pow(10,Lagrvypukl(A2,C2,A3,C3,w)))
        if w>=(x1+4*h0):
            print(w, pow(10,Linear(A3,C3,w)))
        w=w+h

    return 0


def AbsorpKoeff(f,lambd):
    A1=A1Koeff(f)
    A2=A2Koeff(f)
    C1=C1Koeff(f)
    C2=C2Koeff(f)
    A3=A3Koeff(f)
    C3=C3Koeff(f)

    x1,y1=Peresechenie(A2,C2,0,5)
    x0,y0=Peresechenie(A1,C1,A2,C2)

    hnyu=(6.626E-34)*3E8/lambd/1.6E-19
    w=hnyu

    if w<=(x0-4*h1):
        return pow(10,Linear(A1,C1,w))
    if (w>(x0-4*h1))&(w<(x0+4*h1)):
        return pow(10,Lagrvognut(A1,C1,A2,C2,w))
    if (w>=(x0+4*h1))&(w<(x1-4*h0)):
        return pow(10,Linear(A2,C2,w))
    if (w>=(x1-4*h0))&(w<(x1+4*h0)):
        return pow(10,Lagrvypukl(A2,C2,A3,C3,w))
    if w>=(x1+4*h0):
        if Linear(A3,C3,w)<5.04:
            return pow(10,Linear(A3,C3,w))
        else:
            return pow(10,5.04)
#Pe4at(0.4)


def ReturnArray(f):
    step=2
    lam=[365]
    alfa=[AbsorpKoeff(f,lam[0]*1E-9)]
    b=0
    while lam[b]>260:
        if AbsorpKoeff(f,lam[b]*1E-9)<pow(10, 5.04):
            lam.append(lam[b]-step)
            b=b+1
            alfa.append(AbsorpKoeff(f,lam[b]*1E-9))
        else:
            break
    return (lam, alfa)
            
    
#print(ReturnArray(0.3))

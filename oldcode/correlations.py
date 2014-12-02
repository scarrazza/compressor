#!/usr/bin/python

import lhapdf
from numpy import matrix
from numpy import linalg

#### compute C matrix for pdf with x array
def f(C,pdf,x,q):
    nrep = len(pdf)
    for i in range(len(x)):
        for j in range(len(x)):
            f = 0
            I = f*len(x)+i
            J = f*len(x)+j

            a = b = ab = 0            
            aa = bb = 0
            for r in range(nrep):
                a += pdf[r].xfxQ(f,x[i],q) 
                aa+= pdf[r].xfxQ(f,x[i],q)**2
                b += pdf[r].xfxQ(f,x[j],q)
                bb+= pdf[r].xfxQ(f,x[j],q)**2
                ab+= pdf[r].xfxQ(f,x[i],q)*pdf[r].xfxQ(f,x[j],q)
            a /= nrep
            b /= nrep            
            ab/= nrep
            sa = (aa/nrep - a*a)**0.5
            sb = (bb/nrep - b*b)**0.5

            C[i][j] = (ab-a*b)/(sa*sb)

#### compute C matrix for pdf with x array
def f2(C,pdf,x,q):
    nrep = len(pdf)
    for i in range(-3,4):
        for j in range(-3,4):
            a = b = ab = 0            
            aa = bb = 0
            for r in range(nrep):
                for xi in range(len(x)):
                    a += pdf[r].xfxQ(i,x[xi],q)                
                    b += pdf[r].xfxQ(j,x[xi],q)
                    ab+= pdf[r].xfxQ(i,x[xi],q)*pdf[r].xfxQ(j,x[xi],q)
            a /= nrep
            b /= nrep            
            ab/= nrep

            sa = sb = 0            
            for r in range(nrep):
                s1 = s2 = 0
                for xi in range(len(x)):                    
                    s1 += (pdf[r].xfxQ(i,x[xi],q))
                    s2 += (pdf[r].xfxQ(j,x[xi],q))
                sa += (s1-a)**2
                sb += (s2-b)**2
            sa /= nrep
            sb /= nrep
            sa = sa**0.5
            sb = sb**0.5

            C[i][j] = (ab-a*b)/(sa*sb)

##### PDF sets
pset = lhapdf.getPDFSet("1000rep")
pdf = []
for i in range(1,1001):
    pdf.append(pset.mkPDF(i))

pset2= lhapdf.getPDFSet("1000rep_compressed_50")
pdf2 = []
pdf3 = []
for i in range(1,51):
    pdf2.append(pset2.mkPDF(i))

for i in range(101,201):
    pdf3.append(pset.mkPDF(i))

q = 1.0
x = [
    1.0974987654930569E-005,
    1.5922827933410941E-005,
    2.3101297000831580E-005,
    3.3516026509388410E-005,
    4.8626015800653536E-005,
    7.0548023107186455E-005,
    1.0235310218990269E-004,
    1.4849682622544667E-004,
    2.1544346900318823E-004,
    3.1257158496882353E-004,
    4.5348785081285824E-004,
    6.5793322465756835E-004,
    9.5454845666183481E-004,
    1.3848863713938717E-003,
    2.0092330025650459E-003,
    2.9150530628251760E-003,
    4.2292428743894986E-003,
    6.1359072734131761E-003,
    8.9021508544503934E-003,
    1.2915496650148829E-002,
    1.8738174228603830E-002,
    2.7185882427329403E-002,
    3.9442060594376556E-002,
    5.7223676593502207E-002,
    8.3021756813197525E-002,
    0.10000000000000001,
    0.11836734693877551,
    0.13673469387755102,
    0.15510204081632653,
    0.17346938775510204,
    0.19183673469387758,
    0.21020408163265308,
    0.22857142857142856,
    0.24693877551020407,
    0.26530612244897961,
    0.28367346938775512,
    0.30204081632653063,
    0.32040816326530613,
    0.33877551020408170,
    0.35714285714285710,
    0.37551020408163271,
    0.39387755102040811,
    0.41224489795918373,
    0.43061224489795924,
    0.44897959183673475,
    0.46734693877551026,
    0.48571428571428565,
    0.50408163265306127,
    0.52244897959183678,
    0.54081632653061229,
    0.55918367346938780,
    0.57755102040816331,
    0.59591836734693870,
    0.61428571428571421,
    0.63265306122448983,
    0.65102040816326534,
    0.66938775510204085,
    0.68775510204081625,
    0.70612244897959175,
    0.72448979591836737,
    0.74285714285714288,
    0.76122448979591839,
    0.77959183673469379,
    0.79795918367346941,
    0.81632653061224492,
    0.83469387755102042,
    0.85306122448979593,
    0.87142857142857133,
    0.88979591836734695,
    0.90816326530612246 
    ]

Cprior = [[0 for i in range(len(x))] for j in range(len(x))]
Ccompr = [[0 for i in range(len(x))] for j in range(len(x))]
Crando = [[0 for i in range(len(x))] for j in range(len(x))]

"""
f(Cprior,pdf,x,q)
f(Ccompr,pdf2,x,q)
f(Crando,pdf3,x,q)
"""
f2(Cprior,pdf,x,q)
f2(Ccompr,pdf2,x,q)
f2(Crando,pdf3,x,q)

# compute erf
ERFcomp = 0
ERFrando = 0
for i in range(len(x)):
    for j in range(len(x)):
        ERFcomp += (Ccompr[i][j]-Cprior[i][j])**2
        ERFrando+= (Crando[i][j]-Cprior[i][j])**2

print "\n***********************"
print " Summary: compressed", len(pdf2), "replicas"
print "          random", len(pdf3), "replicas"

print "ERFcomp", ERFcomp
print "ERFrando", ERFrando

A = matrix( Cprior )
B = matrix( Ccompr )
C = matrix( Crando )
#print "Prior trace: ", float((A.I*A).trace()) / len(x)
#print "Comp trace: ", float((A.I*B).trace()) / len(x)
#print "Rand trace: ", float((A.I*C).trace()) / len(x)

evA, wA = linalg.eig(A)
evB, wB = linalg.eig(B)
evC, wC = linalg.eig(C)

erfB = 0
erfC = 0
for i in range(len(x)):
    erfB += (evA[i] - evB[i])**2
    erfC += (evA[i] - evC[i])**2

print "ERFcomp eigval:", erfB
print "ERFrand eigval:", erfC

import os
import numpy as np
import scipy

#chosenPath = "/Users/nihao/Dropbox/Documents/Career/Calculs/twistedBilayer/GBN/KondaRelaxedStructures/finalGBNFigureMapUsingMassFactor/55x55ID_usingExtep_QUIP_substrate/"
chosenPath = os.getcwd()
os.chdir(chosenPath)

# -------- Step 1 ---------- #

name = "generateInitWithBN.xyz"

f = open(name, "r")
g = open(name + ".temp", "w")

lineN = 0
for line in f:
    lineN = lineN + 1
    a = line.split()
    if (lineN == 1):
        shiftX = float(a[0])
    if (lineN == 2):
        shiftY = float(a[1])
        shiftXY = float(a[0])
    if a[0]=="C":
        g.write("{} {} {} {}\n".format(3, a[1], a[2], a[3]))
    #elif a[0]=="N":
    elif (a[0]=="N"):# and float(a[3]) > 46):
        g.write("{} {} {} {}\n".format(1, a[1], a[2], a[3]))
    #elif a[0]=="B":
    elif (a[0]=="B"):# and float(a[3]) > 46):
        g.write("{} {} {} {}\n".format(2, a[1], a[2], a[3]))
    #elif (len(a)==4 and float(a[3])<46):
    #    pass
    else:
        g.write(line)

f.close()
g.close()

data = np.genfromtxt(name + ".temp", skip_header=4)

layerIndex = data[:,0]
xdata = data[:,1]
ydata = data[:,2]
zdata = data[:,3]

xdataTop0 = xdata[layerIndex==3]
ydataTop0 = ydata[layerIndex==3]
zdataTop0 = zdata[layerIndex==3]

xdataBottom0 = xdata[layerIndex!=3]
ydataBottom0 = ydata[layerIndex!=3]
zdataBottom0 = zdata[layerIndex!=3]


# -------- Step 2 ---------- #

os.chdir(chosenPath)

name = "generateWithBN.xyz"

f = open(name, "r")
g = open(name + ".temp", "w")

lineN = 0
for line in f:
    lineN = lineN + 1
    a = line.split()
    if (lineN == 1):
        shiftX = float(a[0])
    if (lineN == 2):
        shiftY = float(a[1])
        shiftXY = float(a[0])
    if a[0]=="C":
        g.write("{} {} {} {}\n".format(3, a[1], a[2], a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftXY, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftXY, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftXY-shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftXY+shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftXY-shiftX, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftXY+shiftX, float(a[2])-shiftY, a[3]))
    #elif a[0]=="N":
    elif (a[0]=="N"):# and float(a[3]) > 46):
        g.write("{} {} {} {}\n".format(1, a[1], a[2], a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])-shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])+shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])+shiftXY, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])-shiftXY, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])+shiftXY-shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])+shiftXY+shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])-shiftXY-shiftX, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(1, float(a[1])-shiftXY+shiftX, float(a[2])-shiftY, a[3]))
    #elif a[0]=="B":
    elif (a[0]=="B"):# and float(a[3]) > 46):
        g.write("{} {} {} {}\n".format(2, a[1], a[2], a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])-shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])+shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])+shiftXY, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])-shiftXY, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])+shiftXY-shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])+shiftXY+shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])-shiftXY-shiftX, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(2, float(a[1])-shiftXY+shiftX, float(a[2])-shiftY, a[3]))
    #elif (len(a)==4 and float(a[3])<46):
    #    pass
    else:
        g.write(line)

f.close()
g.close()

data = np.genfromtxt(name + ".temp", skip_header=4)

layerIndex = data[:,0]
xdata = data[:,1]
ydata = data[:,2]
zdata = data[:,3]

xdataTop = xdata[layerIndex==3]
ydataTop = ydata[layerIndex==3]
zdataTop = zdata[layerIndex==3]

xdataBottom = xdata[layerIndex!=3]
ydataBottom = ydata[layerIndex!=3]
zdataBottom = zdata[layerIndex!=3]

npts = 100
x = np.linspace(min(xdataTop),max(xdataTop),npts)
y = np.linspace(min(ydataTop),max(ydataTop),npts)

xdiscr = x
ydiscr = y

from scipy.interpolate import griddata

zdataTopInt=griddata((xdataTop,ydataTop),zdataTop,(x[None,:], y[:,None]), method='cubic')
zdataBottomInt=griddata((xdataBottom,ydataBottom),zdataBottom,(x[None,:], y[:,None]), method='cubic')

zdataDiffInt = (zdataTopInt-zdataBottomInt) # interlayer distances


# -------- Step 3 ---------- #
zdataDiscr = zdataDiffInt

xVec = []
yVec = []
zVec = []

i = 0
for xval in xdiscr:
    j = 0
    for yval in ydiscr:
        xVec.append(xval)
        yVec.append(yval)
        zVec.append(zdataDiscr[j,i])
        j = j+1
    i = i+1
    
        
xVect = np.array(xVec)[~np.isnan(np.array(zVec))]
yVect = np.array(yVec)[~np.isnan(np.array(zVec))]
zVect = np.array(zVec)[~np.isnan(np.array(zVec))]

f_interp = scipy.interpolate.griddata((xVect, yVect), zVect, (xdataTop, ydataTop))#,method="nearest")
#plt.figure(figsize=(20,10))
#im = plt.scatter(xdataTop, ydataTop, c=f_interp, vmin=3.2, vmax=3.5) 
#colorbar(im)







# -------- Step 4 ---------- #
f_interp = scipy.interpolate.griddata((xVect, yVect), zVect, (xdataTop0, ydataTop0))#,method="nearest")
#plt.figure(figsize=(20,10))
#plt.scatter(xdataTop0, ydataTop0, c=f_interp, vmin=3.22, vmax=3.5)

os.chdir(chosenPath)
g = open("interlayerDistances.dat","w")

for el in f_interp:
    g.write(str(el) + "\n")

g.close()

#A = np.genfromtxt("generate.xyz", skip_header=4)
A = np.genfromtxt("generate.xyz", skip_header=4)
xdata = A[:,1]
ydata = A[:,2]

B = np.genfromtxt("interlayerDistances.dat")

#plt.figure(figsize=(20,10))
#plt.scatter(xdata, ydata, c=B)


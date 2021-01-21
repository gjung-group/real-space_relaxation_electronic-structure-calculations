import os
import numpy as np
import scipy
#import matplotlib.pyplot as plt

#def latexify(fig_width=None, fig_height=None, columns=1):
#    """Set up matplotlib's RC params for LaTeX plotting.
#    Call this before plotting a figure.
#
#    Parameters
#    ----------
#    fig_width : float, optional, inches
#    fig_height : float,  optional, inches
#    columns : {1, 2}
#    """
#
#    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
#
#    # Width and max height in inches for IEEE journals taken from
#    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf
#
#    assert(columns in [1,2])
#
#    if fig_width is None:
#        fig_width = 3.39 if columns==1 else 6.9 # width in inches
#
#    if fig_height is None:
#        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
#        fig_height = fig_width*golden_mean # height in inches
#
#    #MAX_HEIGHT_INCHES = 7.0
#    #if fig_height > MAX_HEIGHT_INCHES:
#    #    print("WARNING: fig_height too large:" + fig_height + 
#    #          "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
#    #    fig_height = MAX_HEIGHT_INCHES
#
#    params = {'backend': 'ps',
#              'text.latex.preamble': [r'\usepackage{gensymb}'],
#              'axes.labelsize': 8, # fontsize for x and y labels (was 10)
#              'axes.titlesize': 8,
#              'font.size': 8, # was 10
#              'legend.fontsize': 8, # was 10
#              'xtick.labelsize': 8,
#              'ytick.labelsize': 8,
#              'text.usetex': True,
#              'figure.figsize': [fig_width,fig_height],
#              'font.family': 'serif',
#              'figure.autolayout': True,
#              'lines.linewidth' : 1.0
#    }
#
#    matplotlib.rcParams.update(params)
#
#
#def format_axes(ax):
#    
#    SPINE_COLOR = 'gray'
#
#    for spine in ['top', 'right']:
#        ax.spines[spine].set_visible(False)
#
#    for spine in ['left', 'bottom']:
#        ax.spines[spine].set_color(SPINE_COLOR)
#        ax.spines[spine].set_linewidth(0.5)
#
#    ax.xaxis.set_ticks_position('bottom')
#    ax.yaxis.set_ticks_position('left')
#
#    for axis in [ax.xaxis, ax.yaxis]:
#        axis.set_tick_params(direction='out', color=SPINE_COLOR)

chosenPath = os.getcwd()
os.chdir(chosenPath)

# -------- Step 1 ---------- #

name = "generateInitWithBN.xyz"

f = open(name, "r")
g = open(name + ".temp", "w")

lineN = 0
# threshold = 50 (choose zcoordinate between layer 1 and layer 2 of hBN if more than one hBN layer)

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
    elif (a[0]=="N"):# and float(a[3]) > threshold): # uncomment if more than 1 hBN layer
        g.write("{} {} {} {}\n".format(1, a[1], a[2], a[3]))
    elif (a[0]=="B"):# and float(a[3]) > threshold): # uncomment if more than 1 hBN layer
        g.write("{} {} {} {}\n".format(2, a[1], a[2], a[3]))
    #elif (len(a)==4 and float(a[3])<46): # uncomment if more than 1 hBN layer
    #    pass # uncomment if more than 1 hBN layer
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

# ------------ displacements ------------

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
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftX, a[2], a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftXY, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftXY, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftXY-shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])+shiftXY+shiftX, float(a[2])+shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftXY-shiftX, float(a[2])-shiftY, a[3]))
        g.write("{} {} {} {}\n".format(3, float(a[1])-shiftXY+shiftX, float(a[2])-shiftY, a[3]))
        #g.write("{} {} {} {}\n".format(3, float(a[1])+shiftY, a[2], a[3]))
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
        #g.write("{} {} {} {}\n".format(1, float(a[1])+shiftY, a[2], a[3]))
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
        #g.write("{} {} {} {}\n".format(2, float(a[1])+shiftY, a[2], a[3]))
    #elif (len(a)==4 and float(a[3])<46):
    #    pass
    else:
        g.write(line)

f.close()
g.close()

initArray = np.genfromtxt("generateInitWithBN.xyz.temp",skip_header=4)
  
initAtomType   = initArray[:,0]
initX          = initArray[:,1]
initY          = initArray[:,2]
initZ          = initArray[:,3]

xdataTop = initX[layerIndex==3]
ydataTop = initY[layerIndex==3]
zdataTop = initZ[layerIndex==3]

npts = 100
x = np.linspace(min(xdataTop),max(xdataTop),npts)
y = np.linspace(min(ydataTop),max(ydataTop),npts)

xdiscr = x
ydiscr = y

finalArray = np.genfromtxt("generateWithBN.xyz.temp",skip_header=4)

finalAtomType   = finalArray[:,0]
finalX          = finalArray[:,1]
finalY          = finalArray[:,2]
finalZ          = finalArray[:,3]

dispArray = np.sqrt(np.square(finalX-initX) + np.square(finalY-initY))# + np.square(finalZ-initZ))
dispArrayX = finalX-initX
dispArrayY = finalY-initY

dispArrayTop = dispArray[[finalAtomType==3.0]]
dispArrayTopX = dispArrayX[[finalAtomType==3.0]]
dispArrayTopY = dispArrayY[[finalAtomType==3.0]]
finalXTop = finalX[[finalAtomType==3.0]]
finalYTop = finalY[[finalAtomType==3.0]]
finalZTop = finalZ[[finalAtomType==3.0]]
initXTop = initX[[finalAtomType==3.0]]
initYTop = initY[[finalAtomType==3.0]]
initZTop = initZ[[finalAtomType==3.0]]

dispArrayBottom = dispArray[[finalAtomType!=3.0]]
dispArrayBottomX = dispArrayX[[finalAtomType!=3.0]]
dispArrayBottomY = dispArrayY[[finalAtomType!=3.0]]
finalXBottom = finalX[[finalAtomType!=3.0]]
finalYBottom = finalY[[finalAtomType!=3.0]]
finalZBottom = finalZ[[finalAtomType!=3.0]]
initXBottom = initX[[finalAtomType!=3.0]]
initYBottom = initY[[finalAtomType!=3.0]]
initZBottom = initZ[[finalAtomType!=3.0]]

# Write out the files for the displacement of carbon atoms
# Distance (sqrt(dx2+dy2))

with open('displacementsTop.txt', 'w') as f1:
    #for x, y, disp in zip(finalXTop, finalYTop, dispArrayTop):
    for x, y, disp in zip(initXTop, initYTop, dispArrayTop):
        f1.write(str(x) + "   " + str(y) + "   " + str(disp) + "\n")
        
with open('displacementsTopX.txt', 'w') as f1:
    #for x, y, disp in zip(finalXTop, finalYTop, dispArrayTopX):
    for x, y, disp in zip(initXTop, initYTop, dispArrayTopX):
        f1.write(str(x) + "   " + str(y) + "   " + str(disp) + "\n")
        
with open('displacementsTopY.txt', 'w') as f1:
    #for x, y, disp in zip(finalXTop, finalYTop, dispArrayTopY):
    for x, y, disp in zip(initXTop, initYTop, dispArrayTopY):
        f1.write(str(x) + "   " + str(y) + "   " + str(disp) + "\n")
        
 
# -------- Step 3 ---------- #

uBNxVecS = []
uBNyVecS = []
uBNVecS = []

for index, (x1, y1) in enumerate(zip(initXTop, initYTop)): # Looping over C atoms
    if ((index//9) % 2 == 0): # assign the same value for each sublattice
        minDist = 1000
        for index2, (x2, y2, x2f, y2f) in enumerate(zip(initXBottom, initYBottom, finalXBottom, finalYBottom)): # Looping over BN atoms
            if ((index2//9) % 2 == 0): # only compare to B atoms
                dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)
                if (dist < minDist):
                    minDist = dist
                    minx2 = x2
                    miny2 = y2
                    minx2f = x2f
                    miny2f = y2f
        uBNxVecS.append(minx2f-minx2) 
        uBNyVecS.append(miny2f-miny2) 
        uBNVecS.append(np.sqrt((minx2f-minx2)**2 + (miny2f-miny2)**2))
    else:
        pass # Removing B sublattice contributions
        

# -------- Step 4 ---------- #

dispArrayTopSL = [el for index, el in enumerate(dispArrayTop) if ((index//9) % 2 == 0)] 
dispArrayTopXSL = [el for index, el in enumerate(dispArrayTopX) if ((index//9) % 2 == 0)] 
dispArrayTopYSL = [el for index, el in enumerate(dispArrayTopY) if ((index//9) % 2 == 0)] 

# -------- Step 5 ---------- #

uG_uBNS = np.array(dispArrayTopSL) - uBNVecS
uG_uBN_xS = np.array(dispArrayTopXSL)- uBNxVecS
uG_uBN_yS = np.array(dispArrayTopYSL) - uBNyVecS

# -------- Step 5 ---------- #

name = "displacementsTopX.txt"

data = np.genfromtxt(name)#, dtype=(np.str, np.float, np.float, np.float))

xdata = data[:,0]
ydata = data[:,1]
zdata = data[:,2]

xdataSL = [el for index, el in enumerate(xdata) if ((index//9) % 2 == 0)] 
ydataSL = [el for index, el in enumerate(ydata) if ((index//9) % 2 == 0)] 

num = 1000

#x = np.linspace(min(xdataTop[0::2]),max(xdataTop[0::2]),num)
#y = np.linspace(min(ydataTop[0::2]),max(ydataTop[0::2]),num)
x = np.linspace(min(xdataTop),max(xdataTop),num)
y = np.linspace(min(ydataTop),max(ydataTop),num)

xdiscr2X = x
ydiscr2X = y

uG_uBN_xS_int=griddata((xdataSL,ydataSL),uG_uBN_xS,(x[None,:], y[:,None]), method='linear')

os.chdir(chosenPath)

name = "displacementsTopY.txt"

data = np.genfromtxt(name)#, dtype=(np.str, np.float, np.float, np.float))

#xdata = data[:,0][0::2]
#ydata = data[:,1][0::2]
zdata = data[:,2][0::2]


x = np.linspace(min(xdataTop),max(xdataTop),num)
y = np.linspace(min(ydataTop),max(ydataTop),num)

xdiscr2Y = x
ydiscr2Y = y

uG_uBN_yS_int=griddata((xdataSL,ydataSL),uG_uBN_yS,(x[None,:], y[:,None]), method='cubic')

os.chdir(chosenPath)

name = "displacementsTop.txt"

data = np.genfromtxt(name)#, dtype=(np.str, np.float, np.float, np.float))

#xdata = data[:,0][0::2]
#ydata = data[:,1][0::2]
zdata = data[:,2][0::2]


x = np.linspace(min(xdataTop),max(xdataTop),num)
y = np.linspace(min(ydataTop),max(ydataTop),num)

xdiscr2 = x
ydiscr2 = y

uG_uBNS_int=griddata((xdataSL,ydataSL),uG_uBNS,(x[None,:], y[:,None]), method='cubic')




# -------- Step 6 ---------- #

#latexify(columns=1, fig_height=1.4)
#import matplotlib.gridspec as gridspec
#
#nrow = 1
#ncol = 3
#
##g1 = gridspec.GridSpec(nrow, ncol)#, height_ratios=[6,1])
##g1.update(wspace=0.0, hspace=0.2) # set the spacing between axes.
#
#fig, ((ax0, ax1, ax2)) = plt.subplots(nrow, ncol, sharex='col', sharey='row')
#
#ax0 = subplot(g1[0])
#ax1 = subplot(g1[1])
#ax2 = subplot(g1[2])
#
#im = ax0.imshow(uG_uBNS_int,origin="lower",aspect="auto")#,extent=[min(xdata), max(xdata), min(ydata), max(ydata)],cmap="gist_ncar")#,extent=[400,800,0,700])
#im = ax1.imshow(uG_uBN_xS_int,origin="lower",aspect="auto")#,extent=[min(xdata), max(xdata), min(ydata), max(ydata)],cmap="gist_ncar")#,extent=[400,800,0,700])
#im = ax2.imshow(uG_uBN_yS_int,origin="lower",aspect="auto")#,extent=[min(xdata), max(xdata), min(ydata), max(ydata)],cmap="gist_ncar")#,extent=[400,800,0,700])
#
#ax0.set_yticklabels([])
#ax1.set_yticklabels([])
#
#ax0.axis('off')
#ax1.axis('off')
#ax2.axis('off')
#
#ax0.set_title(r"$u_G - u_{BN}$")
#ax1.set_title(r"$u_G^x - u_{BN}^x$")
#ax2.set_title(r"$u_G^y - u_{BN}^y$")
#
#ax0.set_aspect('equal')
#ax1.set_aspect('equal')
#ax2.set_aspect('equal')
#
#os.chdir("/Users/nihao/Desktop/")
#fig.savefig("dispTest.png",dpi=400,bbox_inches="tight")
#
#
# -------- Step 7 ---------- #

zdataDiscr = uG_uBN_xS_int

xVec = []
yVec = []
zVec = []

i = 0
for xval in xdiscr2X:
    j = 0
    for yval in ydiscr2Y:
        xVec.append(xval)
        yVec.append(yval)
        zVec.append(zdataDiscr[j,i])
        j = j+1
    i = i+1

xVect2 = np.array(xVec)[~np.isnan(np.array(zVec))]
yVect2 = np.array(yVec)[~np.isnan(np.array(zVec))]
zVect2 = np.array(zVec)[~np.isnan(np.array(zVec))]

f_interp2X = scipy.interpolate.griddata((xVect2, yVect2), zVect2, (xdataTop0[0::2], ydataTop0[0::2]))#,method="nearest")
#plt.figure(figsize=(20,10))
#plt.scatter(xdataTop0[0::2], ydataTop0[0::2], c=f_interp2X)#, vmin=3.22, vmax=3.5) 

zdataDiscr = uG_uBN_yS_int

xVec = []
yVec = []
zVec = []

i = 0
for xval in xdiscr2X:
    j = 0
    for yval in ydiscr2Y:
        xVec.append(xval)
        yVec.append(yval)
        zVec.append(zdataDiscr[j,i])
        j = j+1
    i = i+1

xVect2 = np.array(xVec)[~np.isnan(np.array(zVec))]
yVect2 = np.array(yVec)[~np.isnan(np.array(zVec))]
zVect2 = np.array(zVec)[~np.isnan(np.array(zVec))]

f_interp2Y = scipy.interpolate.griddata((xVect2, yVect2), zVect2, (xdataTop0[0::2], ydataTop0[0::2]))#,method="nearest")
#plt.figure(figsize=(20,10))
#plt.scatter(xdataTop0[0::2], ydataTop0[0::2], c=f_interp2Y)#, vmin=3.22, vmax=3.5) 

zdataDiscr = uG_uBNS_int

xVec = []
yVec = []
zVec = []

i = 0
for xval in xdiscr2X:
    j = 0
    for yval in ydiscr2Y:
        xVec.append(xval)
        yVec.append(yval)
        zVec.append(zdataDiscr[j,i])
        j = j+1
    i = i+1

xVect2 = np.array(xVec)[~np.isnan(np.array(zVec))]
yVect2 = np.array(yVec)[~np.isnan(np.array(zVec))]
zVect2 = np.array(zVec)[~np.isnan(np.array(zVec))]

f_interp2 = scipy.interpolate.griddata((xVect2, yVect2), zVect2, (xdataTop0[0::2], ydataTop0[0::2]))#,method="nearest")
#plt.figure(figsize=(20,10))
#plt.scatter(xdataTop0[0::2], ydataTop0[0::2], c=f_interp2)#, vmin=3.22, vmax=3.5) 

os.chdir(chosenPath)
g = open("displacements.txt","w")

for elX, elY, el in zip(f_interp2X, f_interp2Y, f_interp2):
    g.write(str(elX) + "  " + str(elY) + "  " + str(el) + "\n")
    g.write(str(elX) + "  " + str(elY) + "  " + str(el) + "\n")
     
g.close()

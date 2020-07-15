from __future__ import division	#Python 2 compatibility
import pylab as pl
import numpy as np
import sys
import random

from numpy.polynomial.legendre import leggauss
from scipy.ndimage.filters import gaussian_filter 
from termcolor import colored, cprint


##NOTE:
# Sections with #! are commented out because they probably can be deleted safely, but left for a while to make sure it is true
# 12 Jan 2017


#!def dist(point1,point2):
#!	return pl.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)

#!def gauss(x,mu,sigma):
#!	return 1/(sigma*pl.sqrt(2*pl.pi)) * pl.exp(-0.5*( (x-mu)/sigma)**2)

##Calculate integration points and weights for gaussian integration
#
#@param N	Number of integration points
#@param a	Start value of integration domain
#@param b	End value of integration domain
#
#@return		arrays with integration points (p) and weights (w)
def calcGauss(N,a,b):
	p,w = leggauss(N)
	p = (b-a)/2.*p + (b+a)/2.
	w = (b-a)/2.*w
	return p,w

## Class to contain and manipulate coordinates and radii for an ensemble of cylinders
#
# The purpose of this class is to be able to perform structure factor calculations on 
# a set of 2D coordinates corresponding to the centre of the cylinders. The cylinder radii
# are used for collision checks.
class cylinderEnsemble():
	## Constructor for cylinderEnsemble
	# Note that k and phi as standard is discretized using gaussian integration points. This can be changed
	#	by calling setUniformIntegrationPoints
	#
	#@param klim	Tupple with min/max wave vector (k) values for which the structure factor is calculated
	#@param Nk		Number of discretization points of the range of k
	#@param Nphi	Number of discretization points of phi angle for the structure factor calculations
	#@param bbox	Bounding box for the cylinders. Used for visualisation and filling, but not collision check
	#
	#@return			cylinderEnsemble object
	def __init__(self,klim,Nk,Nphi,bbox):
		self.klim = klim
		self.Nk = Nk
		self.Nphi = Nphi
		self.bbox = bbox

		self.setGaussianIntegrationPoints()
		self.lk,self.wk = calcGauss(Nk,klim[0],klim[1])
		self.lphi,self.wphi = calcGauss(Nphi,0,2*pl.pi)
		self.lphi = self.lphi[pl.newaxis].T
		self.wphi = self.wphi[pl.newaxis].T

		#!self.lk = pl.linspace(klim[0],klim[1],Nk)
		#!self.lphi = pl.linspace(0,2*pl.pi,Nphi)[pl.newaxis].T
		#Initialization of other variables:
		self.Ncyl = 0
		self.Nx = 2*self.Ncyl
		self.lCyl = pl.array([])

	## Change k and phi to be equidistant arrays
	#
	# This makes visualiasation easier, but may give less precise integration.
	#
	# @return none
	def setUniformIntegrationPoints(self):
		#Lazy way to initialise arrays...
		self.setGaussianIntegrationPoints()
		self.lk[:] = pl.linspace(self.klim[0],self.klim[1],self.Nk)
		self.wk[:] = self.lk[1]-self.lk[0]
		self.lphi[:,0] = pl.linspace(0,2*pl.pi,self.Nphi)
		self.wphi[:] = self.lphi[1]-self.lphi[0]

	## Change k and phi arrays to match Gaussian integration points
	#
	# This makes visualiasation more difficult, but may give more precise integration.
	#
	# @return none
	def setGaussianIntegrationPoints(self):
		self.lk,self.wk = calcGauss(self.Nk,self.klim[0],self.klim[1])
		self.lphi,self.wphi = calcGauss(self.Nphi,0,2*pl.pi)
		self.lphi = self.lphi[pl.newaxis].T
		self.wphi = self.wphi[pl.newaxis].T


	## Update the list of cylinders
	#
	#	This function overwrites whatever list of cylinders the object may contain with the function input
	#
	#@param x	1D array containing first all x coordinates for the cylinders and then all y coordinates
	#@param r	(scalar) radius of all cylinders, (array) indiviual radii for all cylinder coordinates
	#
	#@return none
	#@see UpdateCylinders,addCylinders
	def UpdateCylinderFromX(self,x,r):
		lCyl = pl.zeros((x.size//2,3))
		lCyl[:,0] = x[:x.size//2]
		lCyl[:,1] = x[x.size//2:]
		#Lazy implementation so that users
		# don't have to specify r as array
		try:
			lCyl[:,2] = r
		except:
			lCyl[:,2] = [r]*(x.size//2)
		self.updateCylinders(lCyl)

	## Get cylinder coordinates as 1D array
	#
	#@return	x	1D array with cylinder x coordinates followed by cylinder y coordinates. Can be fed back
	#				into UpdateCylinderFromCoordinates
	#@see UpdateCylinderFromX,getCylinders
	def getx(self):
		x = pl.zeros(self.Nx)
		x[:self.Ncyl] = self.lCyl[:,0]
		x[self.Ncyl:] = self.lCyl[:,1]
		return x

	## Get Nx3 array of (x,y,r) values of all N cylinders in container
	#
	#@return		Nx3 array of (x,y,r) values
	def getCylinders(self):
		return self.lCyl

	## Save coordinates and radii to file in txt format
	#
	#@param FOUT	Filename to save coordinates and radii to (first line is a header)
	#@return none
	def saveToTxt(self,FOUT):
		f = open(FOUT,'w')
		f.write("x-pos\ty-pos\tradius\n")
		for elm in self.lCyl:
			f.write("{:.5f}\t{:.5f}\t{:.5f}\n".format(elm[0],elm[1],elm[2]))
		f.close()

	## Replace list of cylinders and radii with new content
	#
	#@param	lCyl		Nx3 array of (x,y,r) values for N new cylinders
	#@return none
	#@see UpdateCylinderFromX,addCylinders
	def updateCylinders(self,lCyl):
		self.lCyl = pl.array(lCyl)
		self.Ncyl = lCyl.shape[0]
		self.Nx = 2*self.Ncyl

	## Add cylinders to the current ensemble	
	#
	# As updateCylinders, but just adding instead of replacing
	#
	#@param lCylAdded		Nx3 array of (x,y,r) values for N new cylinders to be added
	#@see updateCylinders
	def addCylinders(self,lCylAdded):
		if len(lCylAdded) > 0:
			total = pl.concatenate((self.lCyl,lCylAdded))
			self.updateCylinders(total)


	## Method to fill up cylinder voids within the bounding box with cylinders, both for Monodisperse and Polydisperse r. 
	# Clearance between cylinders is allowed for through rclearance.
	#
	#@param r
	#@param rclearance
	#


	def fillUpVoidsPoly(self,r,rsigma,rclearance=None,ff=pl.inf):
		print("Filling up voids starting from "+str(self.Ncyl)+" cylinders...")
		nmesh = 256
		X,Y = pl.meshgrid(pl.linspace(self.bbox[0],self.bbox[1],nmesh),
								pl.linspace(self.bbox[0],self.bbox[1],nmesh))

		#res.imag >0  indicates the occupied space (sphere radius + new sphere radius)
		#res.real indicates summed distance to spheres
		res = pl.zeros(X.shape,dtype='complex')
		aOccupied=0

		for elm in self.lCyl: #maybe put elm[1+2]
			res += pl.sqrt(pl.sqrt((X-elm[0])**2+(Y-elm[1])**2) -((elm[2]+elm[2])*rclearance)+0.j)
			aOccupied+=pl.pi*elm[2]**2
			
		lCylAdd = []

		#try to calculate res_try before the while loop
		#saving all xy and then try to fit a pl.normal in all of them

		while res.imag.min() < 1e-6 and (aOccupied/(((-self.bbox[0]+self.bbox[1])**2)))<ff:
		#while (aOccupied/(((-self.bbox[0]+self.bbox[1])**2)))<ff:	  #it runs until there isn't an imaginary part (i.e. voids to fill) or the ff desired is reached 
			iteration=0
			iteration2=0


			while True: #it picks a random positions and checks if it is possible to place a disk there (if res is real we are safe)
				i1=random.randrange(0, nmesh,1)
				i2=random.randrange(0, nmesh,1)		
				g=res[i1][i2]
				#print(iteration)
				#print(g)
				if (g.imag>1e-6):
					#print("complex")
					iteration+=1
				else:
					#print("real")
					break
					
			idx=(i1,i2)
			x = X[idx]
			y = Y[idx]

			#res += pl.sqrt(pl.sqrt((X-x)**2+(Y-y)**2) -(r+rclearance)+0.j)
			#lCylAdd += [ [x,y,r] ]
			
			#this part if poly
			while True: #once identified the void here we try to fit it with one radius from the normal distribution 

				def normDist(r,sigma): #try to have an always positive normal dist
					x=pl.normal(r,rsigma)
					while x<0:
						x=pl.normal(r,rsigma)
					return x
				rnew=normDist(r,rsigma)
				res_try = pl.zeros(X.shape,dtype='complex')
				
				#print(iteration2)
				for elm in self.lCyl:
					res_try += pl.sqrt(pl.sqrt((elm[0]-x)**2+(elm[1]-y)**2) -((rnew+elm[2])*rclearance)+0.j)
				if (res_try.imag.min() < 1e-6):
					cprint("Okay, we can add it!","green")
					print(aOccupied/(((-self.bbox[0]+self.bbox[1])**2)))
					res += pl.sqrt(pl.sqrt((X-x)**2+(Y-y)**2) -((rnew+rnew)*rclearance)+0.j)
					lCylAdd = [ [x,y,rnew] ]
					aOccupied+=pl.pi*rnew**2
					self.addCylinders(lCylAdd)
					print(self.Ncyl)
					print(aOccupied/(((-self.bbox[0]+self.bbox[1])**2)))
					break
				else:
					cprint("It overlaps","red")
					print(aOccupied/(((-self.bbox[0]+self.bbox[1])**2)))
	
					iteration2+=1

		print("Voids filled so that we have "+str(self.Ncyl)+" cylinders.")
		print(ff)
		print(aOccupied/(((-self.bbox[0]+self.bbox[1])**2)))

		#exit()

	def fillUpVoidsMono(self,r,rclearance=None,fnumber=pl.inf):
		print("Filling up voids starting from "+str(self.Ncyl)+" cylinders...")
		nmesh = 256
		X,Y = pl.meshgrid(pl.linspace(self.bbox[0],self.bbox[1],nmesh),
								pl.linspace(self.bbox[0],self.bbox[1],nmesh))
		#res.imag >0  indicates the occupied space (sphere radius + new sphere radius)
		#res.real indicates summed distance to spheres
		res = pl.zeros(X.shape,dtype='complex')
		aOccupied=0
		if rclearance == None:
			rclearance = r
		for elm in self.lCyl:
			res += pl.sqrt(pl.sqrt((X-elm[0])**2+(Y-elm[1])**2) - (elm[2]+r*rclearance)+0.j)
			aOccupied+=pl.pi*elm[2]**2

		lCylAdd = []

		while res.imag.min() < 1e-6 and (self.Ncyl+len(lCylAdd))<fnumber: #it runs until there isn't an imaginary part (i.e. voids to fill) or the ff desired is reached 
		#while (self.Ncyl+len(lCylAdd))<fnumber:
			res.real[0,0]=pl.inf
			res.real[res.imag>1e-6] = pl.inf
			iteration=0

			while True: #it picks a random positions and checks if it is possible to place a disk there (if res is real we are safe)
				i1=random.randrange(0, nmesh,1)
				i2=random.randrange(0, nmesh,1)
				g=res[i1][i2]
				#print(iteration)
				#print(g)
				if (g.imag>1e-6):
					#print("complex")
					iteration+=1
				else:
					idx=(i1,i2)
					break

			while False: #it places the disk where in the minimum of res (it's faster put the particles start clustering and the optimisation for some SFs does not work)
				idx = pl.argmin(res.real)
				idx = pl.unravel_index(idx,res.real.shape)
				print(iteration)
				iteration+=1
				break
	
			x = X[idx]
			y = Y[idx]

			res += pl.sqrt(pl.sqrt((X-x)**2+(Y-y)**2) -(r+r*rclearance)+0.j)
			lCylAdd += [ [x,y,r] ]
			aOccupied+=pl.pi*r**2
			print(self.Ncyl) #starting disks
			print(len(lCylAdd)) #added disks
			#print(res.imag.min())

		self.addCylinders(lCylAdd)
		print("Voids filled so that we have "+str(self.Ncyl)+" cylinders.")

		cprint('Final ff','cyan')
		print(aOccupied/(((-self.bbox[0]+self.bbox[1])**2)))


	## By calling this method, all cylinders are temporarily assumed to have a fixed radius (fixedR)
	#	and colloding cylinders are then removed from the set of cylinders stored by the object. The 
	#	fixed radius approach was taken to ease implementation and speed. If variable radius is needed,
	#	users are urged to implement it themselves (which should be straight forward using the distMatSq
	#  if speed is not an issue).
	#
	#@param fixedR	Sets the radius that all cylinders are assumed to have for the collision check
	#
	#@return Nothing, but updates lCyl hold by the object
	#
	#@see distMatSq
	def removeCollisions(self):

		#Collision matrix (1 if i,j are colliding)
		cMat = self.distMatSq(onlyUpperTriangle=True)
		d=self.radMatSq(onlyUpperTriangle=True)
				
		cMat[cMat>d] = 0
		cMat[cMat>0] = 1	#Set remaining elements to 1
	
		removalList = []
		print ("Elements before collision removal",self.lCyl.shape[0])
		while True:	
			nCollisions = cMat.sum(0)+cMat.sum(1)
			#Remove (first occurence of) particle with max collision number
			idx = pl.argmax(nCollisions)
			if nCollisions[idx] == 0:
				break
			cMat[idx,:] = 0
			cMat[:,idx] = 0
			removalList += [idx]

		lCyl = pl.delete(self.lCyl,pl.s_[removalList],axis=0)
		self.updateCylinders(lCyl)
		print ("Elements after collision removal",self.lCyl.shape[0])

	## Plot histogram of the distances from each element to their six nearest neighbours
	#
	#@param show	If False, then plt.show() is not called
	def sixNearestHistogram(self,show=True):
		dMatSq = self.distMatSq()
		lDistMean = []
		lClosest = []
		for i in range(dMatSq.shape[0]):
			dist = dMatSq[i,:]
			sixClosestMean = pl.sqrt(sorted(dist)[1:7]).mean()
			lDistMean += [sixClosestMean]
			lClosest += [pl.sqrt(sorted(dist)[1])]

		#print lClosest
		print (pl.array(lClosest).mean())

		xmax = 1
		bins = 80
		pl.hist(pl.array(lDistMean),bins=bins,range=(0,xmax),normed=True)
		pl.xlim(0,xmax)
		if show:
			pl.show()


	##Returns squared distances between cylinders as a matrix
	#
	#@param onlyUpperTriangle	If this value is True, then the values are stored as an upper triangle matrix
	#		only
	#
	#@return A matrix where the entries (i,j) and (j,i) contains the distance between cylinder i and j. If 
	#				onlyUpperTriangle is set to True, then all values in the lower triangle are set to zero.
	def distMatSq(self,onlyUpperTriangle=False):
		X,Y = pl.meshgrid(self.lCyl[:,0],self.lCyl[:,0])
		xdiff = Y-X
		
		X,Y = pl.meshgrid(self.lCyl[:,1],self.lCyl[:,1])
		ydiff = Y-X

		idxtri = pl.tril_indices(Y.shape[0],0)
		if onlyUpperTriangle:
			xdiff[idxtri] = 0
			ydiff[idxtri] = 0
		return xdiff**2+ydiff**2
	
	##Return a matrix where the distance between the two centres of the cyl is calcutaled. Useful when the radius is not constant
	def radMatSq(self,onlyUpperTriangle=False):	
		rad=self.lCyl[:,2]
		Drad=pl.zeros((rad.size,rad.size))
		
		for i in range(rad.size):
			for j in range(rad.size):
				Drad[i,j]=(rad[i]+rad[j])**2
	
		return Drad
		

	## Calculates the structure factor and plots it. Parameters identical to polarPlotS
	#
	#@see polarPlotS
	def polarPlot(self,show=True,clim=(0,1),smooth=False,colorbar=True):
		S = self.calcS()
		S[:,0:3]=0
		if clim == "auto":
			#clim = (0.007,S.max()-0.01)
			clim = (0.001,S.max())
		self.polarPlotS(S,clim,smooth,colorbar)
		if show:
			pl.show()
		pl.close

	## Plot a given structure factor (S) as a polar plot
	#
	#@param S			The structure factor to be plotted (dimensions of k and phi should fit what is returned
	#							for the particular instance of the object, so most likely calculated by calcS)
	#@param clim		Limits of the coloraxis for the plot given as a tupple. If set to "auto", then it is
	#							scaled from 0 to max(S)
	#@param smooth		Smooth the structure factor before plotting. This is mostly for testing, and is not
	#							really recommended to change from default value
	#@param colorbar	(bool) Show colorbar in plot
	#@param useAxis	Specify an axis on which to plot the 
	#  plot (defaults to creating a new axis)
	#
	#@return nothing
	#							
	#@see polarPlot
	def polarPlotS(self,S,clim=(0,1),smooth=False,colorbar=True,useAxis=False):
		if smooth:
			print("Smoothing is work in progress (since it is not on a square grid...)")
			#And no reason to import the filter every time...
			S = gaussian_filter(S,smooth)

		if useAxis:
			ax = useAxis
		else:
			ax = pl.subplot(1,1,1,projection='polar',aspect=1.)
		cax = ax.pcolormesh(self.lphi.ravel(),self.lk,S.T,cmap='gist_yarg',vmin=clim[0],vmax=clim[1])
		#cax = ax.pcolormesh(self.lphi.ravel(),self.lk,S.T,cmap='RdBu',vmin=clim[0],vmax=clim[1])
		pl.ylim(0,self.lk[-1])
		if colorbar:
			pl.gcf().colorbar(cax)


	## Method to visualise the cylinders
	#
	# Creates a figure (on the current plot axis) that displays all cylinders with their respective
	# radii. 
	#
	#@param	sameColour		If True, then all cylinders are plotted in the same colour, which is more
	#								visually appealing, but it makes it harder to track their movement.
	#@param	show				If False, then plt.show() will not be called from within the method
	#
	#@return nothing
	def visualise(self,sameColour=False,show=True):
		#By gca() it is easier to integrate in subplots
		pl.gca().set_xlim(self.bbox[0],self.bbox[1])
		pl.gca().set_ylim(self.bbox[0],self.bbox[1])
		pl.gca().set(adjustable='box',aspect='equal')
		if sameColour:
			colourwheel = list("k"*len(self.lCyl))
		else:
			colourwheel = list("bgrcmyk"*(len(self.lCyl)//7+1))
		for cyl in self.lCyl:
			circle = pl.Circle( (cyl[0],cyl[1]), cyl[2], color=[0,0,0,0.3])
			pl.gca().add_artist(circle)
		if show:
			pl.show()
			


	## integration over k-range
	#
	# Integrates S (calculated by calcS) over the whole range of k specified by klim.
	#
	#@param S	The structure factor calculated by calcS
	#
	#@return		S integrated over k (but not over phi)
	#
	#@see calcS, integratePhi, integrateKPhi
	def integrateK(self,S):
		return pl.dot(S,self.wk)

	## integration over phi-range
	#
	# Integrates S (calculated by calcS) over phi from 0 to 360 degrees.
	#
	#@param S	The structure factor calculated by calcS
	#
	#@return		S integrated over phi (but not over k)
	#
	#@see calcS, integrateK, integrateKPhi
	def integratePhi(self,S):
		return 1/(2*pl.pi) * pl.dot(S.T,self.wphi)

	## integration over k- and phi-range
	#
	# Integrates S (calculated by calcS) over all of phi and k.
	#
	#@param S	The structure factor calculated by calcS
	#
	#@return		S integrated over phi and k (scalar value)
	#
	#@see calcS, integrateK, integratePhi
	def integrateKPhi(self,S):
		return self.integratePhi(self.integrateK(S).T)


	## Calculate structure factor over all specified value of phi and k (unless deltaK is specified)
	#
	#@param deltaK		If this 2D vector is specified, then S is calculated only for that vector
	#
	#@return Structure factor (S) for all specified k and phi values
	#
	#@see calcdS
	def calcS(self,deltaK=None):
		#Single deltaK
		#S = 0
		#for pos in self.lCyl:
		#	pos = pl.array(pos[:2])
		#	S += pl.exp(-1.j *pl.dot(deltaK,pos))
		#S = S*S.conj()
		#S /= len(self.lCyl)**2
		#x1/x2-based deltaK
		#def _calcS(lCyl,deltaK):
		#	S = pl.exp(-1.j * pl.dot(deltaK,lCyl[:,:2].T))
		#	S = abs(S.sum())**2/len(lCyl)**2
		#	return S
		#S = 0.*x1
		#for i in range(S.shape[0]):
		#	for j in range(S.shape[1]):
		#		S[i,j] = _calcS(self.lCyl,pl.array([x1[i,j],x2[i,j]]).T)
		if deltaK == None:
			x1 = pl.multiply(self.lk,pl.cos(self.lphi))
			x2 = pl.multiply(self.lk,pl.sin(self.lphi))
			deltaK = pl.zeros((x1.size,2))
			deltaK[:,0] = x1[:].ravel()
			deltaK[:,1] = x2[:].ravel()


		S = pl.exp(-1.j * pl.dot(deltaK,self.lCyl[:,:2].T))
		#If deltaK is just one point, we are summing over axis 0
		try:
			S = abs(S.sum(1))**2/len(self.lCyl)**2
		except ValueError:
			S = abs(S.sum())**2/len(self.lCyl)**2
		#If deltaK is calculated from x1,x2 where we know the
		# original shape, then return S in that shape. Otherwise
		# deltaK is just assumed to be a list
		try:
			S = S.reshape(x1.shape)
		except:
			pass
		return S

	## Calculates the gradient of S for all k and phi values with respect to all cylinder positions
	#
	#@param deltaK		This parameter works the same way as for calcS
	#
	#@return Gradients of the structure factor (dS) for all cylinder positions. First for all x-values
	#			and then for all y-values. This ordering is the same as can be retrieved from getx.
	#
	#@see calcS,getx
	def calcdS(self,deltaK=None):
		if not (deltaK == None):
			N = len(self.lCyl)
			dS = pl.zeros(2*N,dtype='complex')	#For both x and y
			sumpart = 0
			for i in range(N):
				pos = pl.array(self.lCyl[i][:2])
				sumpart += pl.exp(1.j *pl.dot(deltaK,pos))

				dS[i] = -1.j*deltaK[0]*pl.exp(-1.j*pl.dot(deltaK,pos))
				dS[i+N] = -1.j*deltaK[1]*pl.exp(-1.j*pl.dot(deltaK,pos))

			dS = 2*(sumpart*dS).real
			dS /= N**2
			return dS

		else:
			x1 = pl.multiply(self.lk,pl.cos(self.lphi))
			x2 = pl.multiply(self.lk,pl.sin(self.lphi))
			deltaK = pl.zeros((x1.size,2))
			deltaK[:,0] = x1[:].ravel()
			deltaK[:,1] = x2[:].ravel()

			dS = pl.zeros((x1.shape[0],x1.shape[1],self.Nx),dtype='complex')	#For both x and y
			deltaK = deltaK.reshape((x1.shape[0],x1.shape[1],2))
			#--Non-vectorised implementation for reference:
			#for m in range(x1.shape[0]):
			#	for n in range(x1.shape[1]):
			#		sumpart = 0
			#		for i in range(self.Ncyl):
			#			pos = pl.array(self.lCyl[i][:2])
			#			sumpart += pl.exp(1.j *pl.dot(deltaK[m,n,:],pos))
			#
			#			dS[m,n,i] = -1.j*deltaK[m,n,0]*pl.exp(-1.j*pl.dot(deltaK[m,n,:],pos))
			#			dS[m,n,i+self.Ncyl] = -1.j*deltaK[m,n,1]*pl.exp(-1.j*pl.dot(deltaK[m,n,:],pos))
			#
			#		dS[m,n,:] = 2*(sumpart*dS[m,n,:]).real
			#dS /= self.Ncyl**2
			tmp = pl.dot(deltaK,self.lCyl[:,:2].T)
			tmp1 = pl.exp(-1.j*tmp).T
			tmp2 = pl.exp(1.j*tmp)
			sumpart = tmp2.sum(2)

			dS[:,:,:self.Ncyl] = (-1.j*deltaK[:,:,0].T*tmp1).T
			dS[:,:,self.Ncyl:] = (-1.j*deltaK[:,:,1].T*tmp1).T
			dS = 2*pl.multiply((sumpart.T[pl.newaxis]).T,dS).real
			dS /= self.Ncyl**2
			return dS

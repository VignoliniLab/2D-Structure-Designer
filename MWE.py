from __future__ import division
import pylab as pl
from matplotlib.ticker import MaxNLocator
from scipy.optimize import fmin_l_bfgs_b

from cylinderEnsemble import cylinderEnsemble
#from mmaClass import mmaClass	#This library is not publicly available
from powertools import outputFolder
from termcolor import colored, cprint
from scipy.stats import norm
import os


_FD = 0
_collisionConstraint = 0
_show = 0
_useMMA = False	#True: MMA, False: BFGS
movelim = 0.04	#40 nm should be enough..
NcollisionRuns = 2 #2 suggested value
saveSteps=1
showFinal=1
showDist=0
showInitial=1
randomNoOpt=0
Optimisation=1


def main(FOUT):
	seed = int(pl.rand()*1e8)
	OUT = outputFolder(FOUT)
	#OUT = ""
	print("Seed: {:}\t\t(now this is actually correctly implemented)".format(seed))
	pl.seed(seed)
	maxiter =51 #51 suggested value
	Nk = 101
	bbox = (-2.5,2.5)
	Nphi = 91
	itshift = 0
	fillFactor = 1.001 #defines the clearance used to then fill up the voids
	k0 = 2*pl.pi/(2*r)
	klim = (1.,6/4*k0)
	ce = cylinderEnsemble(klim,Nk,Nphi,bbox)

	sigmaK =   0.2
	sigmaPhi = 6/180*pl.pi

	#Decide which SF you want to target

	#Hex
	#Sref = sum([gaussPeak(ce,k0,phi0/180*pl.pi,sigmaK,sigmaPhi) for phi0 in [0,60,120,180,240,300]])
	#Sref = Sref/ce.integrateKPhi(Sref)
	
	#Rect (excluding some corners)
	#Sref = sum([gaussPeak(ce,k0,phi0/180*pl.pi,sigmaK,sigmaPhi) for phi0 in [0,90,180,270]])

	#Elliptical
	Sref = elliptical(ce,k0,k0,1)
	Sref = Sref/ce.integrateKPhi(Sref)

	#Layered
	#aF=2.17
	#Sref = gaussPeak(ce,k0,0/180*pl.pi,sigmaK,sigmaPhi)+gaussPeak(ce,k0,180/180*pl.pi,sigmaK,sigmaPhi)+gaussPeak(ce,k0/aF,90/180*pl.pi,sigmaK,sigmaPhi)+gaussPeak(ce,k0/aF,270/180*pl.pi,sigmaK,sigmaPhi)
	#Sref = Sref/ce.integrateKPhi(Sref)


	if 1: #put 1 if you want to visualise the targeted SF
		Sref /= Sref.max()
		ce.polarPlotS(Sref)
		pl.savefig(OUT+"TargetSF.pdf")
		pl.close('all')
		#pl.show()
		#exit()

	fnum = int((ff*((-bbox[0]+bbox[1])**2))/(pl.pi*((r)**2))) #number of cylinders for given ff
	print("Desired Cyl")
	print(fnum)

	def normDist(r,rsigma): #try to have an always positive normal dist
		x=pl.normal(r,rsigma)
		while x<0:
			x=pl.normal(r,rsigma)
		return x
	lCyl = pl.array([ ( (pl.random()-0.5)*2*bbox[0], (pl.random()-0.5)*2*bbox[1], normDist(r,rsigma) ) for i in range(fnum)]) #smaller than ff

	ce.updateCylinders(lCyl)
	ce.radMatSq(lCyl)

	ce.removeCollisions()
	if 0: #put 1 if you want to check that the Collisions Removal worked well
		print(ce.lCyl.shape[0])
		ce.visualise(sameColour=True,show=False)
		pl.savefig(OUT+"Collision.pdf")
		pl.close
	


	#Decide how to FillUp: Mono or Poly. The Poly one uses ff as input (it has to be calculated at every step!)
	
	ce.fillUpVoidsPoly(r,rsigma,fillFactor,ff) 
	#ce.fillUpVoidsMono(r,fillFactor,fnum)


	if showInitial or randomNoOpt: #put 1 if you want to check that the FillingUp worked well
		if 0:
			a=ce.Ncyl
			ce.removeCollisions()
			b=ce.Ncyl
			cprint(a-b,"red")
			ce.visualise(sameColour=False,show=False)
			cprint('Desired Cyl','cyan')
			print(fnum)
			cprint('Final Cyl','cyan')
			print(ce.lCyl.shape[0])	

		ce.saveToTxt("/Users/macbookpro/anaconda3/Script/2D generator/structureFactorDesigner/"+"2D"+str(index)+".txt")
		if showInitial:  
			ce.polarPlot(clim="auto",show=False)
			pl.savefig(OUT+"InitialSF.pdf")
			pl.close()

		if showDist:
			cprint("Plotting radius distribution","cyan")
			(mu, sigma) = norm.fit(ce.lCyl[:,2])
			xmin, xmax = pl.xlim()
			x = pl.linspace(xmin, xmax, 1000)
			p = norm.pdf(x, mu, sigma)
			pl.figure()
			pl.hist(ce.lCyl[:,2],bins=15,normed=True)
			pl.plot(x, p, 'k', linewidth=2)
			pl.xlabel('Radius (µm)')
			pl.ylabel('Probability')
			pl.title(r'$\mathrm{Radii\ distribution\:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
            #pl.xlim ((0,0.6))
			pl.savefig(OUT+"RadiusDist.pdf")
			#pl.show()
			print(mu)
			print(sigma)


	if Optimisation:
		if saveSteps: #do you want to save the starting guess?
			ce.visualise(show=False)
			pl.savefig(OUT+"start.pdf") #saves the starting configuration in the output folder
			pl.close()

		for ii in range(NcollisionRuns):
			if ii > 0:
				
				ce.removeCollisions()

				ce.fillUpVoidsPoly(r,rsigma,fillFactor,ff) 

				itshift += 100

			x0 = ce.getx()
			N = x0.size
			M = 1

			iterHist = []
			if _useMMA:	#we go for the else
				mma = mmaClass(N,M,xmin=bbox[0],xmax=bbox[1],movelim=movelim,volfrac=1.)
				xmma,f0,df0dx,f,dfdx = mma.initializeInputs()
				xmma[:] = x0
				if _FD:
					mma.FDcheck(elm=[0,1,5,6,3]);maxiter=50

				for it in range(1,maxiter):
						f0[:],df0dx[:] = interfaceCalcSAnddS(xmma,Sref,ce,plot=False,saveIter=itshift+it,iterHist=iterHist,outFolder=OUT)
						if _collisionConstraint:
							f[0],dfdx[:,0] = interfaceCollisionConstraint(xmma,ce)
						else:
							f[0],dfdx[:,0] = mma.volfrac(xmma)
						print ("{:}\t{:.3f}\t\t{:.3f}".format(str(it+itshift).zfill(4),f0[0],f[0]))

						xmma[:] = mma.mmasub(it,xmma,f0,df0dx,f,dfdx)
						if abs(mma.xold2-xmma).max() < 1e-3:
							break
				xres = xmma
				f = f0[0]

			else:
				xbfgs = x0+0. #bfgs is the name of the algorithm
				it = 1
				def _wrapperFunc(xonly):
					return interfaceCalcSAnddS(xonly,Sref,ce,plot=False,saveIter=itshift+it,iterHist=iterHist,outFolder=OUT)

				if 0: #FD check
					leps = pl.linspace(0.00001,0.1,20)
					abserror=pl.zeros(len(leps))
					idxabs=0
					idxelm = 7 #do both for x and y(second half)
					fref,dfref = _wrapperFunc(xbfgs)
					for eps in leps:
						xperturb = xbfgs+0.
						xperturb[idxelm] += eps
						fnew,tmp = _wrapperFunc(xperturb)
						error = abs((fnew-fref)/eps-dfref[idxelm])
						#print(error)
						abserror[idxabs]=error
						idxabs+=1
					#print(abserror)
					pl.scatter(leps,abserror)
					pl.plot(leps,abserror)
					pl.ylabel('Absolute Error')
					pl.xlabel('Epsilon')
					pl.annotate('idx=7 (x point), cSep=0.1', xy=(0.05, 0.95), xycoords='axes fraction')
					#pl.legend('idx=1 (x point), cSep=2')
					pl.savefig(OUT+"FD_x7_cSep=0.1.pdf")
					pl.show()
					
					#exit()
				xres,f0,d = fmin_l_bfgs_b(_wrapperFunc,xbfgs,bounds=[bbox]*xbfgs.size,maxiter=maxiter,iprint=1)
				#xres = xbfgs #why?
				f = f0
					

		#print xres
		pl.save(OUT+"xres",xres)
		print ("Final val:", interfaceCalcSAnddS(xres,Sref,ce,True)[0])

		if saveSteps: 
			ce.saveToTxt(OUT+"2D"+str(index)+".txt")  

		ce.removeCollisions()
		ce.saveToTxt("/Users/macbookpro/anaconda3/Script/2D generator/structureFactorDesigner/"+"2D"+str(index)+".txt")  

		if showFinal: # put 1 if you want to plot the final results
			print("Plotting final structure factor")
			ce.polarPlot(clim="auto",show=False)
			pl.savefig(OUT+"finishSF.pdf")
			pl.close() 

			print("Plotting final design after removal of collisions")
			pl.close
			ce.visualise(sameColour=False,show=False)
			pl.savefig(OUT+"finish.pdf")
			pl.close() 

		cprint("Desired Cyl",'cyan')
		print(ff)
		aOccupied=0
		for elm in ce.lCyl:
			aOccupied+=pl.pi*elm[2]**2
		cprint('Final ff','cyan')
		print(aOccupied/(((-ce.bbox[0]+ce.bbox[1])**2)))



def gaussPeak(ce,k0,phi0,sigmaK,sigmaPhi):
	tmp = ce.lphi-phi0
	tmp[tmp>pl.pi] = (tmp[tmp>pl.pi] - 2*pl.pi)
	tmp[tmp<-pl.pi] = (tmp[tmp<-pl.pi] + 2*pl.pi)
	Ktrans,Phitrans = pl.meshgrid(ce.lk-k0,tmp)

	return 1/(sigmaK*sigmaPhi*pl.sqrt(2*pl.pi))*pl.exp(-0.5 * ((Ktrans/sigmaK)**2+(Phitrans/sigmaPhi)**2))

##This function is used to generate the reference structure factor by calculating a 2D ellipse
def elliptical(ce,a,b,sigmaE):
    K,Phi = pl.meshgrid(ce.lk,ce.lphi)
    result = pl.zeros(K.shape)
    idx1 = K > a*b/pl.sqrt( (b*pl.cos(Phi))**2+(a*pl.sin(Phi))**2) #> if you want a disk, < if you want a filled ellipse
    idx2 = K-(sigmaE) < a*b/pl.sqrt( (b*pl.cos(Phi))**2+(a*pl.sin(Phi))**2)
    result[idx1*idx2] = 1
    return result

def randomSF (ce):
	K,Phi = pl.meshgrid(ce.lk,ce.lphi)
	result = pl.random(K.shape)
	return result

def interfaceCalcSAnddS(x,Sref,ce,plot=False,saveIter=False,iterHist=[],maxIter=False,outFolder=""):
	r = ce.lCyl[:,2] 
	ce.UpdateCylinderFromX(x,r)

	lk,wk = ce.lk,ce.wk #calcGauss(Nk,klim[0],klim[1])
	phi,wphi = ce.lphi,ce.wphi

	#Sref = ce.Sref

	S = ce.calcS()
	dS = ce.calcdS()

	#integration over k and phi for elements in dS
	def dintegrateKPhi(dS):
		tmp = (ce.wk[pl.newaxis].T*dS).sum(1)
		return 1/(2*pl.pi)*(ce.wphi*tmp).sum(0)

	SA = ce.integrateKPhi(S)
	dSA = dintegrateKPhi(dS)
	#return SA,dSA

	SN = S/SA
	dSN = (dS*SA- pl.swapaxes(S[pl.newaxis].T*dSA[pl.newaxis],0,1))/SA**2
	#return SN[20,30],dSN[20,30,:]

	#Miminize average difference
	cost = ce.integrateKPhi(abs(SN-Sref)**2)
	dcost = 2*dintegrateKPhi(pl.swapaxes((SN-Sref)[pl.newaxis].T,0,1)*dSN)
	#return cost,dcost

	#Minimize standard deviation
	c1 = 100
	A = ce.integrateKPhi(S*0+1)
	mu = cost/A
	dmu = dcost/A
	SD = SN-Sref
	dSD = dSN

	STD = c1*ce.integrateKPhi((SD-mu)**2)/A 
	dSTD = 2*c1*dintegrateKPhi((dSD-dmu)*pl.swapaxes((SD-mu)[pl.newaxis].T,0,1))/A
	#print "\t\t",cost

	iterHist += [STD]
	#cprint(iterHist,'red')
	if saveSteps:
		pl.close()
		#Subplots explained:
		#https://stackoverflow.com/questions/3584805/in-matplotlib-what-does-the-argument-mean-in-fig-add-subplot111
		fig = pl.figure()
		fig.add_subplot(131)
		ce.visualise(show=False)
		pl.xticks([])
		pl.yticks([])
		fig.add_subplot(132,projection='polar')
		ce.polarPlotS(S,clim=(0,0.5),colorbar=False,useAxis=pl.gca())
		pl.xticks([])
		pl.yticks([])
		fig.add_subplot(133)
		pl.plot(range(len(iterHist)),iterHist)
		pl.plot(len(iterHist)-1,iterHist[-1],'*')
		pl.ylim(0,max(iterHist))
		if maxIter:
			pl.xlim(0,maxIter)
		else:
			pl.xlim(0,len(iterHist))
		x0,x1 = pl.gca().get_xlim()
		y0,y1 = pl.gca().get_ylim()
		pl.gca().set_aspect(abs(x1-x0)/abs(y1-y0))
		pl.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
		if len(iterHist) > 2:
			pl.xticks([0,len(iterHist)-1])
		else:
			pl.xticks([0])
		#Save the figure in a parallel thread while continuing the optimization
		#(and pray that we finish before pl.close())
		#... this doesn't work. Maybe because the thread doesn't have a window manager?
		#OUT2 = "/Users/vej22/Documents/Bacteria/Modelling/TopOptStructureFactor/"+OUT
		#thread = Thread(target=fig.savefig,kwargs={'fname':OUT2+"iter"+str(saveIter).zfill(4)+".png",
		#	'bbox_inches':'tight','dpi':150})
		#thread.start()
		fig.savefig(outFolder+"iter"+str(saveIter).zfill(4)+".png",bbox_inches='tight',dpi=150)
		pl.close()
	#needed for polydisperse
	while True:
		pen,dpen = interfaceCollisionConstraint(x,ce)
		cSep =0.002 #how strong we penalise (IT NEEDS TO BE OPTIMISED FOR EVERY STRUCTURE YOU HAVE IN MIND)
		return STD+cSep*pen,dSTD+cSep*dpen
		break
	return STD,dSTD
	

def interfaceCollisionConstraint(x,ce): 
	r = ce.lCyl[:,2]
	ce.UpdateCylinderFromX(x,r)
	N = ce.Ncyl
	f = 0
	df = pl.zeros(x.size)
	#d = 0.9/2	#distance between objects
				# (I have given a 0.1 slack, since heavi(0) = 0.5)
	beta = 4
	def heavi(x):
		return -0.5*(pl.tanh(beta*x)/pl.tanh(beta)-1)
	def dHeavi(x):
		return -0.5*beta*(1-pl.tanh(beta*x)**2)/pl.tanh(beta)

	if 0:
		x = pl.linspace(-3,3)
		pl.plot(x,heavi(x))
		pl.show()
		exit()

	#c = 1e-3	#Set low for FD check
	c = 100
	fillFactor=1.01  #maybe try smaller, play with it
	lCylNew = ce.lCyl
	for i in range(N):
		for j in range(i+1,N):
			#dist = (lCylNew[i][0]-lCylNew[j][0])**2 + (lCylNew[i][1]-lCylNew[j][1])**2
			#radiiDist=((lCylNew[i][2]+lCylNew[j][2])*fillFactor)**2
			dist = pl.sqrt((lCylNew[i][0]-lCylNew[j][0])**2 + (lCylNew[i][1]-lCylNew[j][1])**2)
			radiiDist=((lCylNew[i][2]+lCylNew[j][2])*fillFactor)
			ftmp = c*(dist-radiiDist)
			#print(dist)
			#print(radiiDist)
			#print(ftmp)
			f += heavi(ftmp)	
			df[i]   += dHeavi(ftmp)*c*2*(lCylNew[i][0]-lCylNew[j][0])
			df[j]   -= dHeavi(ftmp)*c*2*(lCylNew[i][0]-lCylNew[j][0])
			df[i+N] += dHeavi(ftmp)*c*2*(lCylNew[i][1]-lCylNew[j][1])
			df[j+N] -= dHeavi(ftmp)*c*2*(lCylNew[i][1]-lCylNew[j][1])
	
	#exit()
	return f ,df




if __name__ == "__main__":
	indexMin=1
	indext=1
	indexMax=2
	for i in range(indext,indexMax):
		index=i
		step=4
		if index<(indexMin+step):
			r = 0.15
			rsigma=r*5/100
			ff=0.4
			main("Rad150nm_Poly5_IsotropicSF_Sk02Sp6_Ff04_5µm"+str(index))

	os.system('say "your program has finished"')

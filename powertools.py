from __future__ import division
import subprocess
import shutil
import time
import sys
import os


def outputFolder(RESDIR,force=False,logfile='log.out'):
	if RESDIR[-1] != "/":
		RESDIR += "/"
	if os.path.exists(RESDIR):
		if not force:
			try:
				#Python 2
				ans = raw_input("The dir " + RESDIR + " already exists, continue? (y/n)")
			except NameError:
				#Python 3
				ans = input("The dir " + RESDIR + " already exists, continue? (y/n)")
			if not ans=='y':
				exit("Quitting instead of overwriting")
	else:
		os.makedirs(RESDIR)
	print("Resultdir is: " + RESDIR)
	shutil.copyfile(os.path.abspath(sys.argv[0]),RESDIR+sys.argv[0])

	tee = subprocess.Popen(['tee', '%s/%s' % (RESDIR, logfile)],
			stdin=subprocess.PIPE) 
	sys.stdout.flush()
	os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
	sys.stderr.flush()
	os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

	return RESDIR

#_ittime = pl.zeros(2)	#First entry is updater per printIterTime, second per printSubTime
def printIterTime():
	tic = time.time()
	if not _ittime[0] == 0.:
		print("Iteration time: {:.2f}s".format(tic-_ittime[0]))
	_ittime[0]=tic
	_ittime[1]=tic

#To be used "within" printIterTime
def printSubTime(description):
	tic = time.time()
	if not _ittime[1] == 0.:
			print(description + " time: {:.2f}s".format(tic-_ittime[1]))
	_ittime[1]=tic
	
#_subtimestring = [""]
def gatherSubTime(description):
	tic = time.time()
	if not _ittime[1] == 0.:
			_subtimestring[0] += description + " time: {:.2f}s\n".format(tic-_ittime[1])
	_ittime[1]=tic

def printGatheredTimes():
	print(_subtimestring[0])
	_subtimestring[0] = ""
	

##
#
# A dictionary class with a couple of extra features:
#  * If an unindexed entry is accessed an empty list is returned,
#		which means that commands like
#			DH["myVar"] += [2]
#		does not require initialisation of DH["myVar"]. This is convenient
#		for looping
#
#	* entries created using const() are accessed as an attribute. For example
#		DH.const("myVar",2)
#		print DH.myVar
class dataHandler(dict):
	def __init__(self, *args, **kwargs):
		super(dataHandler,self).__init__(*args,**kwargs)
		self._const = []

	def __getitem__(self,key):
		if self.has_key(key):
			return super(dataHandler,self).__getitem__(key)
		else:
			self[key] = []
			return self[key]
	
	def const(self,key,val):
		if not key in dir(self):
			self._const += [key]
			setattr(dataHandler,key,val)
		else:
			print("dataHandler: name '{:s}' is already in use".format(key))

	def consts(self,dic):
		for key in dic:
			self.const(key,dic[key])

	def dump(self,PATH):
		f = open(PATH,'w')
		for key in self._const:
			data = getattr(dataHandler,key)
			f.write("{:s} !{:s} !{:s} !{:s}\n".format(key.ljust(15),data.__str__().ljust(15),type(key),type(data)))
		f.close()


if __name__ == "__main__":
	a = dataHandler()
	a.const("it",5)
	a.consts({"it"  :5,
				"it2" :3})
	print(a.it)
	print(a.it2)
	a["ab"] += [2]
	print(a["ab"])

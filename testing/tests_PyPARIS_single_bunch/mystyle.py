from matplotlib import rc, rcdefaults
import pylab as pl



def mystyle(fontsz=10):
	
	font = {#'family' : 'normal',
			#'weight' : 'bold',
			'size'   : fontsz}
	print(fontsz)
	rcdefaults()
	rc('font', **font)
	
def mystyle_arial(fontsz=10, dist_tick_lab=7):
	
	rcdefaults()
	rc('font',**{'family':'sans-serif','sans-serif':['arial'], 'size':fontsz})
	rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=dist_tick_lab)
	


def sciy():
	pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
	
def scix():
	pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='x') 

def extr(filename):
	with open(filename) as f:
		words = [line.split() for line in f]
		entropy = open( 'entropy.dat', 'w')
		entropy.write( str(words[1][8]) )
		entropy.close()
	


extr('entropycal.dat')	

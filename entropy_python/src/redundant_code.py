    
    
    # Don't need to reshape the covariance matrix following lines are redundant (to be removed)
    #
    # reshape the covariance matrix to 3 x (3n)^2/3
    covar.matar = np.loadtxt('covarmat.dat')
    if len(covar.matar[0]) > 3:
	    natoms = len(covar.matar) / 3
	    reshar = covar.matar.reshape(1, (3*natoms)**2)
	    cleng = ((3*natoms)**2) / 3
	    covar.matar = reshar.reshape(cleng, 3) * 0.01
	    np.savetxt('covarmat_undiag.dat', covar.matar)


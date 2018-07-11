import numpy as np
H=np.array([[2, 0], [0, 2]])

[w,x0]=np.linalg.eig(H)

Gamma0Tdiag = 1. / (1 + np.exp(np.diag(w)))
Gamma0Tdiag = np.diag(Gamma0Tdiag)

x=np.array([x0[0]])
Gamma0Trealspace = x * Gamma0Tdiag * x.T

#quenched Hamiltonian -> van Hemmen
#H = H0 + deltalambda * diag(quenchlineuniform)
#define t
expiHt = expm(1j * H * t)

expminusiHt = expiHt.conjugate().T

Gammat = expminusiHt * Gamma0Trealspace * expiHt
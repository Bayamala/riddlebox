%% The following consists of excerpts of a code I used to time - evolve free fermions, described by hoppings and chemical potential without pairing:

%% H0 is the initial hopping Hamiltonian written as a matrix and H is the same with a slightly shifted chemical potential (which enter as diagonal terms):

%% Diagonalize Hamiltonian Matrix

import numpy
import scipy

[V, D] = eig(H0)

%% Compute correlation functions at temperature T

Gamma0Tdiag = 1. / (1 + exp(beta * diag(D)))
Gamma0Tdiag = diag(Gamma0Tdiag)

%% This gives the correlations in the diagonalized space(i.e., contains only occupations), so I rotate back into real space:

Gamma0Trealspace = V * Gamma0Tdiag * V.conjugate().T

% Time evolution

H = H0 + deltalambda * diag(quenchlineuniform)

expiHt = expm(1j * H * t)

expminusiHt = expiHt.conjugate().T

% one time step:

Gammat = expminusiHt * Gamma0Trealspace * expiHt

% % With pairing term, one has to be a bit more careful.

% starting Hamiltonian (BCS Theory and chemical potential)

% H = sum_{i, j} Jij c_i ^ + c_j
% + 1 / 2 sum_{i, j} Kij(c_i ^ + c_j ^ + c_i c_j)      Kij = Jij i < j, -Jij i > j
% +mu sum_i c_i ^ + c_i

A = -J * (Jij + Jij.conjugate().T) + mu*eye(N)

B=-J * (Jij - Jij.conjugate().T) %has to be adjusted if pairing term different strength than hopping, see van Hemmens paper for definition of A and B

% % the following should achieve the same as the diagonalization using van Hemmens code

% % singular value decomposition to find eigenstates of H0:

[PhiDag0, Lambda0, PsiTDag0] = svd(A0-B0)
% PhiDag * Lambda * PsiTDag.conjugate().T=A-B

E0modes = diag(A0 - Lambda0)

E0gs = trace(A0 - Lambda0) / 2

    Phi = PhiDag.conjugate().T

    Psi = PsiT.conjugate().T

% % % transformation between real space and eigenstates

    g = (Phi + Psi) / 2

    h = (Phi - Psi) / 2

%how to understand this?
    T = kron([[1, 0],[0, 0]], g)+kron([[0, 0],[0, 1]], g)+kron([[0, 1],[0, 0]], h)+kron([[0, 0],[1, 0]], h)

% % % % % I just realize I did not finish this part where I extract the covariance matrix from g and h, but a thermal state should look like:

    eigenmodeoccupation = diag(1. / (exp(diag(Lambda) / temperature) + 1))
    Gammathermal = T.conjugate().T*(kron([[1,0],[0,0]],eye(chainlength,chainlength)-eigenmodeoccupation)+kron([[0,0],[0,1]],eigenmodeoccupation))*T
    %how does this formular come about?

% % % % % % % % % Lambda=Lambda0?

% Time evolution of covariance matrix in the diagonalized Hamiltonian

Uteohelp = expm(1j * deltat * Lambda)

Uteo = kron([[1, 0],[0, 0]], Uteohelp)+kron([[0, 0],[0, 1]], Uteohelp.conjugate().T)

% % % time evolution two times Gammat?

Gammat = T.conjugate().T * Uteo.conjugate().T * T * Gammat * T.conjugate().T * Uteo * T

% Observables from covariance matrix

    cidagcj = Gammat(chainlength + 1:2 * chainlength, chainlength + 1: 2 * chainlength)

    cicjdag = Gammat(1:chainlength, 1: chainlength)

    cicj = Gammat(1:chainlength, chainlength + 1: 2 * chainlength)

    cidagciMatrix = kron(diag(real(cidagcj)), ones(1, chainlength))

    cidagcicjdagcjMatrix = cidagciMatrix. * cidagciMatrix.'

    occ = real(diag(cidagcj))
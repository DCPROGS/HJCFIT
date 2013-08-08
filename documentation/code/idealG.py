from numpy import dot, identity, abs, all
from dcprogs.likelihood import QMatrix, IdealG, expm

qmatrix = QMatrix([ [-3050,        50,  3000,      0,    0], 
                    [2./3., -1502./3.,     0,    500,    0], 
                    [   15,         0, -2065,     50, 2000], 
                    [    0,     15000,  4000, -19000,    0], 
                    [    0,         0,    10,      0,  -10] ], 2)

idealG = IdealG(qmatrix)
print idealG

idealG_fa = dot(expm(2e-4 * qmatrix.ff), qmatrix.fa)
assert all( abs(idealG_fa - idealG.fa(2e-4)) < 1e-8 )

inversion = -0.5 * identity(2) - qmatrix.aa
assert all( abs( dot(inversion, idealG.laplace_af(-0.5)) - qmatrix.af ) < 1e-8 )

import time
import math
import sys
import numpy as np
from scipy.optimize import minimize

from dcpyps import dcio
from dcpyps import dataset
from dcpyps import mechanism
from dcprogs.likelihood import Log10Likelihood

# LOAD DATA: Burzomato 2004 example set.
scnfiles = [["../../DCPYPS/dcpyps/samples/glydemo/keA.scn"], 
            ["../../DCPYPS/dcpyps/samples/glydemo/keB.scn"],
            ["../../DCPYPS/dcpyps/samples/glydemo/keC.scn"], 
            ["../../DCPYPS/dcpyps/samples/glydemo/keD.scn"]]
tres = [0.000040, 0.000040, 0.000040, 0.000040]
tcrit = [0.0006, -0.070, -0.020, -0.020]
conc = [300e-6, 3000e-6, 30000e-6, 100000e-6]

recs = []
bursts = []
for i in range(len(scnfiles)):
    rec = dataset.SCRecord(scnfiles[i], conc[i], tres[i], tcrit[i])
    rec.record_type = 'recorded'
    recs.append(rec)
    bursts.append(rec.bursts.intervals())
    rec.printout()

# LOAD PRIMED MECHANISM USED in Lape et al 2012
mecfn = "../../DCPYPS/dcpyps/samples/mec/ke05.mec"
version, meclist, max_mecnum = dcio.mec_get_list(mecfn)
mec = dcio.mec_load(mecfn, meclist[1][0])

# PREPARE RATE CONSTANTS.
rates = mec.unit_rates()
mec.set_rateconstants(rates)

# Fixed rates.
for i in range(len(mec.Rates)):
    mec.Rates[i].fixed = False

# Constrained rates.
mec.Rates[20].is_constrained = True
mec.Rates[20].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[20].constrain_args = [24, 1.5]
mec.Rates[25].is_constrained = True
mec.Rates[25].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[25].constrain_args = [21, 2.0]
mec.Rates[30].is_constrained = True
mec.Rates[30].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[30].constrain_args = [34, 3.0]
mec.Rates[32].is_constrained = True
mec.Rates[32].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[32].constrain_args = [34, 2.0]
mec.Rates[33].is_constrained = True
mec.Rates[33].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[33].constrain_args = [31, 2.0]
mec.Rates[35].is_constrained = True
mec.Rates[35].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[35].constrain_args = [31, 3.0]
mec.update_constrains()

mec.set_mr(True, 15, 0)
mec.set_mr(True, 23, 1)
mec.set_mr(True, 29, 2)
mec.update_mr()

mec.printout(sys.stdout)
theta = np.log(mec.theta())

kwargs = {'nmax': 2, 'xtol': 1e-12, 'rtol': 1e-12, 'itermax': 100,
    'lower_bound': -1e6, 'upper_bound': 0}
likelihood = []

for i in range(len(recs)):
    likelihood.append(Log10Likelihood(bursts[i], mec.kA,
        recs[i].tres, recs[i].tcrit, **kwargs))

def dcprogslik(x, args=None):
    mec.theta_unsqueeze(np.exp(x))
    lik = 0
    for i in range(len(conc)):
        mec.set_eff('c', conc[i])
        lik += -likelihood[i](mec.Q) * math.log(10)
    return lik

iternum = 0
def printiter(theta):
    global iternum
    iternum += 1
    if iternum % 100 == 0:
        lik = dcprogslik(theta)
        print("iteration # {0:d}; log-lik = {1:.6f}".format(iternum, -lik))
        print(np.exp(theta))

lik = dcprogslik(theta)
print ("\nStarting likelihood = {0:.6f}".format(-lik))
start = time.clock()
wallclock_start = time.time()
success = False
result = None
while not success:
    #res = minimize(dcprogslik, np.log(theta), method='Powell', callback=printit, options={'maxiter': 5000, 'disp': True})
    result = minimize(dcprogslik, theta, method='Nelder-Mead', callback=printiter,
        options={'xtol':1e-4, 'ftol':1e-4, 'maxiter': 5000, 'maxfev': 10000,
        'disp': True})
    if result.success:
        success = True
    else:
        theta = result.x
        
end = time.clock()
wallclock_end = time.time()
print ("\nHJCFIT Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
        %time.localtime()[0:6])
print ('CPU time in simplex=', end - start)
print ('Wallclock time in simplex=', wallclock_end - wallclock_start)
print ('\n\nresult=')
print (result)

print ('\n Final log-likelihood = {0:.6f}'.format(-result.fun))
print ('\n Number of iterations = {0:d}'.format(result.nit))
print ('\n Number of evaluations = {0:d}'.format(result.nfev))
mec.theta_unsqueeze(np.exp(result.x))
print ("\n Final rate constants:")
mec.printout(sys.stdout)
print ('\n\n')

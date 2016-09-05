from dcpyps import mechanism
from dcprogs.mpihelpers import MPILikelihoodSolver

# LOAD DATA: Burzomato 2004 example set.
scnfiles = [["../samples/glydemo/A-10.scn"],
            ["../samples/glydemo/B-30.scn"],
            ["../samples/glydemo/C-100.scn"],
            ["../samples/glydemo/D-1000.scn"]]
tres = [0.000030, 0.000030, 0.000030, 0.000030]
tcrit = [0.004, -1, -0.06, -0.02]
conc = [10e-6, 30e-6, 100e-6, 1000e-6]

mysolver = MPILikelihoodSolver()
mysolver.load_data(scnfiles, tres, tcrit, conc)

recs = mysolver.recs
bursts = mysolver.bursts

# LOAD FLIP MECHANISM USED Burzomato et al 2004
mecfn = "../samples/mec/demomec.mec"
# PREPARE RATE CONSTANTS.
rates = [4500.0, 700.0, 2500.0, 1800.0, 900.0, 18000.0, 200.0,
         0.1100E+06, 4900.0, 0.4000E+09, 1850.0, 10000.0, 5000.0,
         0.7500E+09, 8500.0, 1050.0, 3500.0, 0.5000E+07, 2300.0,
         0.9500E+07, 1950, 0.130000E+08]
mysolver.load_mec(mecfn, rates)

for i in range(len(mysolver.mec.Rates)):
    mysolver.mec.Rates[i].fixed = False

# Constrained rates.
mysolver.mec.Rates[21].is_constrained = True
mysolver.mec.Rates[21].constrain_func = mechanism.constrain_rate_multiple
mysolver.mec.Rates[21].constrain_args = [17, 3]
mysolver.mec.Rates[19].is_constrained = True
mysolver.mec.Rates[19].constrain_func = mechanism.constrain_rate_multiple
mysolver.mec.Rates[19].constrain_args = [17, 2]
mysolver.mec.Rates[16].is_constrained = True
mysolver.mec.Rates[16].constrain_func = mechanism.constrain_rate_multiple
mysolver.mec.Rates[16].constrain_args = [20, 3]
mysolver.mec.Rates[18].is_constrained = True
mysolver.mec.Rates[18].constrain_func = mechanism.constrain_rate_multiple
mysolver.mec.Rates[18].constrain_args = [20, 2]
mysolver.mec.Rates[8].is_constrained = True
mysolver.mec.Rates[8].constrain_func = mechanism.constrain_rate_multiple
mysolver.mec.Rates[8].constrain_args = [12, 1.5]
mysolver.mec.Rates[13].is_constrained = True
mysolver.mec.Rates[13].constrain_func = mechanism.constrain_rate_multiple
mysolver.mec.Rates[13].constrain_args = [9, 2]
mysolver.mec.update_constrains()

mysolver.mec.set_mr(True, 7, 0)
mysolver.mec.set_mr(True, 15, 1)

mysolver.mec_printout()

mysolver.set_likelihood_func()

mysolver.run_optimizer()

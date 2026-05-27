from HJCFIT.likelihood import QMatrix, IdealG, MissedEventsG

# Define parameters.
qmatrix = QMatrix([ [-3050,        50,  3000,      0,    0], 
                    [2./3., -1502./3.,     0,    500,    0], 
                    [   15,         0, -2065,     50, 2000], 
                    [    0,     15000,  4000, -19000,    0], 
                    [    0,         0,    10,      0,  -10] ], 2)
tau = 1e-4
 
eG = MissedEventsG(qmatrix, tau)
idealG = IdealG(qmatrix)

tcritical = 5e-3

print("Equilibrium Vectors\n"            \
      "=======================\n\n"          \
      "Ideal Likelihood\n"                   \
      "----------------\n\n"                 \
      "  * initial: {ideal_initial!r}\n"     \
      "  * final: {ideal_final!r}\n\n\n"     \
      "Missed-events Likelihood\n"           \
      "------------------------\n\n"         \
      "  * initial: {equi_initial!r}\n"      \
      "  * final: {equi_final!r}\n\n\n\n"    \
      "CHS Vectors\n"                    \
      "===============\n\n"                  \
      "Missed-events Likelihood\n"           \
      "------------------------\n\n"         \
      "  * tcritical: {tcritical}\n"         \
      "  * initial: {chs_initial!r}\n"       \
      "  * final: {chs_final!r}"             \
      .format(
        ideal_initial = idealG.initial_vectors,
        ideal_final   = idealG.final_vectors,
        equi_initial  = eG.initial_vectors,
        equi_final    = eG.final_vectors,
        chs_initial   = eG.initial_CHS_vectors(tcritical),
        chs_final     = eG.final_CHS_vectors(tcritical),
        tcritical     = tcritical
      )
)

from behave import when, then
from test_setup import register_type
register_type()


@given('a list of {n:Integer} random ideal likelihoods')
def step(context, n):
  from dcprogs.random import qmatrix as random_qmatrix
  from dcprogs.likelihood import IdealG
  qmatrices, Gs, i = [], [], 10*n
  while len(Gs) != n:
    i -= 1
    if i == 0: raise AssertionError('Could not instanciate enough likelihoods.')
    qmatrix = random_qmatrix()
    try: G = IdealG(qmatrix)
    except: continue
    else:
      Gs.append(G)
      qmatrices.append(qmatrix)
  if not hasattr(context, 'qmatrices'): context.qmatrices = []
  if not hasattr(context, 'likelihoods'): context.likelihoods = []
  context.qmatrices.extend(qmatrices)
  context.likelihoods.extend(Gs)




@when('IdealG objects are instantiated with the q-matrices')
def step(context):
  from dcprogs.likelihood import IdealG
  context.idealgs = [IdealG(u) for u in context.qmatrices]

@when('IdealG objects are instantiated with the matrices and nopens')
def step(context):
  from dcprogs.likelihood import IdealG
  context.idealgs = [IdealG(u.matrix, u.nopen) for u in context.qmatrices]

@when('the {name} equilibrium occupancies are computed')
def step(context, name): 
  if not hasattr(context, 'occupancies'): context.occupancies = []
  equname = '{0}_occupancies'.format(name)
  for G in context.likelihoods:
    context.occupancies.append( getattr(G, equname) )




@then('computing af for each time yields exp(t Q_AA) Q_AF')
def step(context):
  from numpy import abs, all, dot
  from scipy.linalg import expm
  for idealg, matrix in zip(context.idealgs, context.qmatrices):
    for t in context.times:
      value = dot(expm(t * matrix.aa), matrix.af)
      try: assert all(abs(idealg.af(t) - value) < context.tolerance)
      except:
        print(matrix)
        print(t)
        raise

@then('computing fa for each time yields exp(t Q_FF) Q_FA')
def step(context):
  from numpy import abs, all, dot
  from scipy.linalg import expm
  for idealg, matrix in zip(context.idealgs, context.qmatrices):
    for t in context.times:
      value = dot(expm(t * matrix.ff), matrix.fa)
      try: assert all(abs(idealg.fa(t) - value) < context.tolerance)
      except:
        print(matrix)
        print(t)
        raise

@then('computing laplace_af for each scale yields (sI - Q_AA)^-1 Q_AF')
def step(context):
  from numpy import abs, all, dot, identity
  from numpy.linalg import inv
  for idealg, matrix in zip(context.idealgs, context.qmatrices):
    for scale in context.scales:
      value = dot(inv(scale * identity(matrix.aa.shape[0]) - matrix.aa), matrix.af)
      try: assert all(abs(idealg.laplace_af(scale) - value) < context.tolerance)
      except:
        print(matrix)
        print(scale)
        raise

@then('computing laplace_fa for each scale yields (sI - Q_FF)^-1 Q_FA')
def step(context):
  from numpy import abs, all, dot, identity
  from numpy.linalg import inv
  for idealg, matrix in zip(context.idealgs, context.qmatrices):
    for scale in context.scales:
      value = dot(inv(scale * identity(matrix.ff.shape[0]) - matrix.ff), matrix.fa)
      try: assert all(abs(idealg.laplace_fa(scale) - value) < context.tolerance)
      except:
        print(matrix)
        raise


@then('the initial occupancies exists and is the kernel of I - laplace_af * laplace_fa')
def step(context):
  from numpy.linalg import inv, svd
  from numpy import abs, all, dot, identity
  for matrix, idealg in zip(context.qmatrices, context.idealgs):
    occupancies = idealg.initial_occupancies
    kernel = dot( dot(inv(matrix.aa), matrix.af), dot(inv(matrix.ff), matrix.fa) )
    kernel = identity(kernel.shape[0]) - kernel
    U, singvals, V = svd(kernel)

    try:
      assert sum(abs(singvals) < context.tolerance) == 1
      assert all(dot(occupancies, kernel) < context.tolerance)
    except:
      print(matrix)
      print("Equilibrium: {0}".format(occupancies))
      print("Kernel Application: {0}".format(dot(occupancies, kernel)))
      raise
    
@then('the final occupancies exists and is the kernel of I - laplace_fa * laplace_af')
def step(context):
  from numpy.linalg import inv, svd
  from numpy import abs, all, dot, identity
  for matrix, idealg in zip(context.qmatrices, context.idealgs):
    occupancies = idealg.final_occupancies
    kernel = dot( dot(inv(matrix.ff), matrix.fa), dot(inv(matrix.aa), matrix.af) )
    kernel = identity(kernel.shape[0]) - kernel
    U, singvals, V = svd(kernel)

    try:
      assert sum(abs(singvals) < context.tolerance) == 1
      assert all(dot(occupancies, kernel) < context.tolerance)
    except:
      print(matrix)
      print("Equilibrium: {0}".format(occupancies))
      print("Kernel Application: {0}".format(dot(occupancies, kernel)))
      raise

@then('the {name} equilibrium occupancies are the only solution to the equilibrium equations')
def step(context, name):
  from numpy.linalg import svd
  from numpy import dot, identity, abs, all
  for qmatrix, G, occ in zip(context.qmatrices, context.likelihoods, context.occupancies):
    eqmatrix = dot(G.laplace_af(0), G.laplace_fa(0)) if name == "initial"                         \
               else dot(G.laplace_fa(0), G.laplace_af(0))
    eqmatrix -= identity(eqmatrix.shape[0])

    left, sings, right = svd(eqmatrix)
    sings = sorted(sings)

    try:
      assert abs(sings[0]) < context.tolerance
      assert sings[0] < sings[1] + context.tolerance
      assert all(abs(dot(occ,  eqmatrix)) < context.tolerance)
    except: 
      print(G)
      print(" * occupancies: {0}".format(occ))
      print(" * matrix:\n{0}".format(eqmatrix))
      print(" * error: {0}".format(dot(occ, eqmatrix)))
      print(" * singular values: {0}".format(sings))
      raise

@then('the components of the {name} equilibrium occupancies sum to one')
def step(context, name):
  from numpy import sum
  for qmatrix, G, occ in zip(context.qmatrices, context.likelihoods, context.occupancies):
    assert (sum(occ) - 1e0) < context.tolerance


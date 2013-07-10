from behave import when, then
from test_setup import register_type
register_type()


@when('IdealG objects are instantiated with the q-matrices')
def step(context):
  from dcprogs.likelihood import IdealG
  context.idealgs = [IdealG(u) for u in context.matrices]
@when('IdealG objects are instantiated with the matrices and nopens')
def step(context):
  from dcprogs.likelihood import IdealG
  context.idealgs = [IdealG(u.matrix, u.nopen) for u in context.matrices]

@then('computing af for each time yields exp(t Q_AA) Q_AF')
def step(context):
  from numpy import abs, all, dot
  from scipy.linalg import expm
  for idealg, matrix in zip(context.idealgs, context.matrices):
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
  for idealg, matrix in zip(context.idealgs, context.matrices):
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
  for idealg, matrix in zip(context.idealgs, context.matrices):
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
  for idealg, matrix in zip(context.idealgs, context.matrices):
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
  for matrix, idealg in zip(context.matrices, context.idealgs):
    occupancies = idealg.occupancies_initial
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
  for matrix, idealg in zip(context.matrices, context.idealgs):
    occupancies = idealg.occupancies_final
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

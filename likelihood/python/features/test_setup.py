""" Holds some general setup for behave. """

def Matrix(string): 
  """ Creates matrices from specific strings """
  from numpy import array, identity
  if string == "Qmatrix":
    return array([[ -3050,        50,  3000,      0,    0 ], 
                  [ 2./3., -1502./3.,     0,    500,    0 ],  
                  [    15,         0, -2065,     50, 2000 ],  
                  [     0,     15000,  4000, -19000,    0 ],  
                  [     0,         0,    10,      0,  -10 ] ])
  if string == "empty": return array([])
  if string == "spam": return identity(5).tolist()
  if string == "numpy_spam": return array([['a', 'b', 'c']*3])
  if string == "rectangle": return identity(6).reshape(3, 12)

def register_type(): 
  from behave import matchers
  matchers.register_type(Integer=lambda x: int(x))
  matchers.register_type(Float=lambda x: float(x))
  matchers.register_type(Eval=lambda x: eval(x))
  matchers.register_type(Bool=lambda x: bool(x))
  matchers.register_type(Matrix=Matrix)

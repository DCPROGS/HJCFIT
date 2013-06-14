from numpy import array

class DeterminantEq(object):
    def __init__(self, matrix, nopen, tau, doopen=True):
        import ctypes
        from numpy import array
        libloc = "/Users/mdavezac/usr/src/dcprogs/build/likelihood/liblikelihood.dylib"
        self._likelihood = ctypes.CDLL(libloc)
        self._likelihood.create_determinant_eq.argtypes = [
            ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_double),  ctypes.c_int,
            ctypes.c_double, ctypes.c_bool ]
        self._likelihood.create_determinant_eq.restype = ctypes.c_void_p
        self._likelihood.delete_determinant_eq.argtypes = [ctypes.c_void_p]
        self._likelihood.call_determinant_eq.argtypes = [ctypes.c_void_p, ctypes.c_double]
        self._likelihood.call_determinant_eq.restype = ctypes.c_double
        self._likelihood.str_determinant_eq.argtypes = [ctypes.c_void_p]
        self._likelihood.str_determinant_eq.restype = ctypes.c_char_p
        
        matrix = array(matrix, dtype="float64")
        self._c_handle = self._likelihood.create_determinant_eq(
             matrix.shape[0], matrix.shape[1],
             matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
             int(nopen), float(tau), bool(doopen))
    def __del__(self): 
        if hasattr(self, '_c_handle'):
          c_handle = self._c_handle
          del self._c_handle
          self._likelihood.delete_determinant_eq(c_handle)
    def __call__(self, s):
        return self._likelihood.call_determinant_eq(self._c_handle, s)
    def __str__(self): 
        return self._likelihood.str_determinant_eq(self._c_handle)
            
matrix = array([[ -3050,        50,  3000,      0,    0 ], 
                [ 2./3., -1502./3.,     0,    500,    0 ],  
                [    15,         0, -2065,     50, 2000 ],  
                [     0,     15000,  4000, -19000,    0 ],  
                [     0,         0,    10,      0,  -10 ] ])
a = DeterminantEq(matrix, 2, 0.15, False)
print a
print a(-17090.192769)


%{
  PyObject* _dcprogs_dtype() {
    return PyArray_TypeObjectFromType(DCProgs::numpy::type<DCProgs::t_real>::value);
  }
  PyObject* eig(PyObject* _matrix) {

    using namespace DCProgs;
    DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
    Eigen::EigenSolver<t_rmatrix> eigsolver(matrix);
    if(eigsolver.info() != Eigen::Success)  {
      std::ostringstream sstr("Could not solve eigenvalue problem.");
      sstr << numpy_io(matrix) << "\n";
      throw errors::Mass(sstr.str());
    }
    t_cvector const eigenvalues = eigsolver.eigenvalues();
    t_cmatrix const eigenvectors = eigsolver.eigenvectors();
    if( (eigenvalues.imag().array().abs() > 1e-12).all() )  {
      Object<> py_eigs = steal_ref(numpy::wrap_to_numpy(eigenvalues));
      if(not py_eigs) throw errors::PythonErrorAlreadyThrown();
      Object<> py_vecs = steal_ref(numpy::wrap_to_numpy(eigenvectors));
      if(not py_vecs) throw errors::PythonErrorAlreadyThrown();
      return PyTuple_Pack(2, py_eigs.release(), py_vecs.release());
    }
    Object<> py_eigs = steal_ref(numpy::wrap_to_numpy(eigenvalues.real().eval()));
    if(not py_eigs) throw errors::PythonErrorAlreadyThrown();
    Object<> py_vecs = steal_ref(numpy::wrap_to_numpy(eigenvectors.real().eval()));
    if(not py_vecs) throw errors::PythonErrorAlreadyThrown();
    return PyTuple_Pack(2, py_eigs.release(), py_vecs.release());
  }
  PyObject* inv(PyObject* _matrix) {
    using namespace DCProgs;
    DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
    return numpy::wrap_to_numpy(matrix.inverse().eval()); 
  }
  double det(PyObject* _matrix) {
    using namespace DCProgs;
    DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
    return static_cast<t_real>(matrix.determinant());
  }
  PyObject* svd(PyObject* _matrix) {

    using namespace DCProgs;
    DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
    Eigen::JacobiSVD<t_rmatrix> svd(matrix, Eigen::ComputeThinU|Eigen::ComputeThinV);
    Object<> py_eigs = steal_ref(numpy::wrap_to_numpy(svd.singularValues().eval()));
    if(not py_eigs) throw errors::PythonErrorAlreadyThrown();
    Object<> pyU = steal_ref(numpy::wrap_to_numpy(svd.matrixU().eval()));
    if(not pyU) throw errors::PythonErrorAlreadyThrown();
    Object<> pyV = steal_ref(numpy::wrap_to_numpy(svd.matrixV().transpose().eval()));
    if(not pyV) throw errors::PythonErrorAlreadyThrown();
    return PyTuple_Pack(3, pyU.release(), py_eigs.release(), pyV.release());
  }
  PyObject* expm(PyObject* _matrix) {
    using namespace DCProgs;
    DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
    return numpy::wrap_to_numpy(matrix.exp().eval());
  }
%}
PyObject* _dcprogs_dtype();
%feature("autodoc", "Computes eigenvalues\n\n"
         "Only have this so we can compute eigenvalues of 128-bit reals, if dcprogs is so "
         "compiled.") eig;
PyObject* eig(PyObject* _matrix);
%feature("autodoc", "Computes inverse of a matrix\n\n"
         "Only have this so we can compute inverse of 128-bit reals, if dcprogs is so "
         "compiled.") inv;
PyObject* inv(PyObject* _matrix);
DCProgs::t_real det(PyObject* _matrix);
PyObject* svd(PyObject* _matrix);
PyObject* expm(PyObject* _matrix);


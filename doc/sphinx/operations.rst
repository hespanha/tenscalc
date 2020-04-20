.. include:: defs.rst

====================================
 Operations on symbolic expressions
====================================

Arithmetic operations
=====================

|tenscalc| supports most of the basic arithmetic operations in
|matlab|, in most cases using the same or a very similar syntax.

.. list-table:: Basic arithmetic operations
   :header-rows: 1
   :stub-columns: 1
   :widths: 2 4 4
   :align: left
      
   * - Usage
     - Description
     - Notes

   * - .. function:: X + Y
       .. function:: plus(X,Y)
     - Entry-by-entry addition of tensors of the same size or of a
       scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``plus()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.
     
   * - .. function:: X - Y
       .. function:: minus(X,Y)
     - Entry-by-entry subtraction of tensors of the same size or
       of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``minus()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

   * - .. function:: X .* Y
       .. function:: times(X,Y)
     - Entry-by-entry multiplication of tensors of the same size or
       of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``times()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

   * - .. function:: X * Y
       .. function:: mtimes(X,Y)
     - Matrix multiplication of tensors that adapts to the size of the
       operands as follows:
       
       * regular matrix multiplication when ``X`` and ``Y`` are both
         matrices (tensors with 2 dimensions) with
         ``size(X,2)==size(Y,1)``
       * matrix by column-vector multiplication when ``X`` is a matrix
         (tensor with 2 dimensions) and ``Y`` a vector (tensor with 1
         dimension) with ``size(X,2)==size(Y,1)``
       * row vector by matrix multiplication when ``X`` is a vector
         (tensor with 1 dimension) and ``Y`` a matrix (tensors with 2
         dimensions) with ``size(X,1)==size(Y,1)``
       * inner product when ``X`` and ``Y`` are both vectors (tensors
         with 1 dimension) with ``size(X,1)==size(Y,1)``
       * entry-by-entry multiplication with either ``X`` or ``Y``
         are scalars (tensor with 0 dimensions).
     - Depending on the sizes of the parameters, this operation may
       behave quite differently from |matlab|'s matrix multiplication
       ``mtimes()``.

   * - .. function:: X ./ Y
       .. function:: rdivide(X,Y)
     - Entry-by-entry right division of tensors of the same size
       or of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``rdivide()``, expansion upon
       singleton dimensions is not performed automatically to match
       matrix sizes.

   * - .. function:: X .\\ Y
       .. function:: ldivide(X,Y)
     - Entry-by-entry left division of tensors of the same size or
       of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``ldivide()``, expansion upon
       singleton dimensions is not performed automatically to match
       matrix sizes.

   * - .. function:: sum(X,vecdim)
       .. function:: sum(X,'all')
     - Sum of entries of the tensor ``X`` along the directions
       specified by the vector ``vecdim``, or over all
       dimensions. Resulting in a vector with the same size as ``X``,
       but with the dimensions in ``vecdim`` removed.
     - Similar syntax to |matlab|

   * - .. function:: max(X)
       .. function:: max(X,[],vecdim)
       .. function:: max(X,[],'all')
     - Maximum entry of the tensor ``X`` along its 1st dimension (1st
       form), the dimensions specified in ``vecdim`` (2nd form), or
       along all dimensions (3rd form).
     - Similar syntax to |matlab|, with the exception that it does not
       return a second output with indices.

   * - .. function:: max(X,Y)
     - For two tensors X and Y with the same size, returns a tensor M
       also with the same size, but with entries taken from X or Y,
       depending on which entry is largest. If X is a scalar, then M
       has the same size as Y and its entries are the largest of the
       corresponding entry of Y or the (only) entry of X. Similarly if
       Y is a scalar.
     - Similar syntax to |matlab|.

   * - .. function:: min(X)
       .. function:: min(X,[],vecdim)
       .. function:: min(X,[],'all')
     - Minimum entry of the tensor ``X`` along its 1st dimension (1st
       form), the dimensions specified in ``vecdim`` (2nd form), or
       along all dimensions (3rd form).
     - Similar syntax to |matlab|, with the exception that it does not
       return a second output with indices.

   * - .. function:: min(X,Y)
     - For two tensors X and Y with the same size, returns a tensor M
       also with the same size, but with entries taken from X or Y,
       depending on which entry is smallest. If X is a scalar, then M
       has the same size as Y and its entries are the smallest of the
       corresponding entry of Y or the (only) entry of X. Similarly if
       Y is a scalar.
     - Similar syntax to |matlab|.

   * - .. function:: full(X)
     - Converts sparse tensor to full in the sense that subsequent
       computations will not take advantage of the information that
       some entries are known to be zero.

     - The main use of this function is to force the results of a
       computation to be returned in a linear memory structure with
       all the zeros filled in appropriately. This is generally
       required to return ``n``-dimensional tensors to |matlab| with
       ``n>2``, since |matlab| does not support sparse arrays with
       more than 2 dimensions.

   * - .. function:: tprod(A1,index1,A2,index2,...)
     - Tensor product
     - See :ref:`tprod-section`

     
.. _tprod-section:

Tensor product
--------------

The function |tprod| provides a very general and flexible
multiplication operation between tensors, which includes the usual
matrix multiplication ``*`` as a special case, the summation over rows
and/or columns |sum|, computing the trace of a matrix |trace|,
extracting a matrix main diagonal |diag|, computing the Euclidean norm
of a vectors |norm|; as well as generalizations of all these
operations to ``n``-dimensional tensors:

.. function:: tprod(A,a,B,b,C,c,...)

   :param A: 1st tensor
   :type A: |tcalculus| symbolic tensor
   :param a: indices for ``A`` with length ``ndims(A)``
   :type a: vector of integers 

   :param B: 2nd tensor
   :type B: |tcalculus| symbolic tensor
   :param b: indices for ``B`` with length ``ndims(B)``
   :type b: vector of integers

   :param C: 3rd tensor
   :type C: |tcalculus| symbolic tensor
   :param c: indices for ``C`` with length ``ndims(C)``
   :type c: vector of integers

|tprod| returns a tensor ``Y`` obtained using a summation-product operation of the form

.. math::

   Y(y_1, y_2, ...) = \sum_{s_1} \sum_{s_2} \cdots A(a_1,a_2,...) B(b_1,b_2,...) C(c_1,c_2,...) \cdots

with the matching between the result tensor indices :math:`y_1,
y_2,...`, the summation indices :math:`s_1, s_2,...`, and the input
tensors indices :math:`a_1, a_2,..., b_1,b_2,...,c_1,c_2,...` is
determines as follows:


* A negative values of ``-1`` in one or several of the input tensor
  indices :math:`a_1, a_2,..., b_1,b_2,...,c_1,c_2,...` means those
  particular indices should be summed under the 1st summation
  operation.

* A negative values of ``-2`` in one or several of the input tensor
  indices :math:`a_1, a_2,..., b_1,b_2,...,c_1,c_2,...` means those
  particular indices should be summed under the 2nd summation
  operation.

  ...

* The absence of a negative index means that there are no summations.

* A positive value of ``+1`` in one or several of the input tensor
  indices means those indices should match the 1st index :math:`y_1`
  of the result tensor :math:`Y`.

* A positive value of ``+2`` in one or several of the input tensor
  indices means those indices should match the 2nd index :math:`y_1`
  of the result tensor :math:`Y`.

  ...

* The absence of a positive index means that the result has an empty
  size (i.e., it is a scalar).

  

A few examples are helpful to clarify the syntax and highlight the
flexibility of |tprod|::

  Tvariable A [m,n];
  Tvariable x [n]
  Tvariable B [n,k]

  tprod(A,[1,-1],x,[-1]);    % same as the matrix-vector product A*x
  tprod(A,[1,-1],B,[-1,2]);  % same as the matrix-matrix product A*B

  tprod(A,[2,1]);            % same as the transpose A'

  tprod(A,[1,-1]);           % same as sum(A,2)

  tprod(x,[-1],x,[-1])       % same as the square of the Euclidean norm x'*x
  tprod(A,[-1,-2],A,[-1,-2]) % same as the square of the Frobenius norm
  
  Tvariable C [n,n]
  tprod(C,[1,1])             % vector with the main diagonal of A
  tprod(C,[-1,-1])           % same as trace(A)
  
Logical-valued operations
=========================

|tenscalc| uses logical-valued operations mostly to specify
optimization constraints.

.. list-table:: Logical-valued operations
   :header-rows: 1
   :stub-columns: 1
   :widths: 3 4 4
   :align: left
      
   * - Usage
     - Description
     - Notes

   * - .. function:: X==Y
       .. function:: eq(X,Y)
     - Entry-by-entry equality comparison of tensors of the same size
       or of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``eq()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

   * - .. function:: X>=Y
       .. function:: ge(X,Y)
     - Entry-by-entry greater than or equal to comparison of tensors
       of the same size or of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``ge()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

       From the perspective of a constrained optimization numerical
       solver, due to finite numerical precision, X>=Y and X>Y
       represent the same constraint.
       
   * - .. function:: X>Y
       .. function:: gt(X,Y)
     - Entry-by-entry greater than comparison of tensors of the same
       size or of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``gt()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

   * - .. function:: X<=Y
       .. function:: le(X,Y)
     - Entry-by-entry smaller than or equal to comparison of tensors
       of the same size or of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``le()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

       From the perspective of a constrained optimization numerical
       solver, due to finite numerical precision, X<=Y and X<Y
       represent the same constraint.

   * - .. function:: X<Y
       .. function:: lt(X,Y)
     - Entry-by-entry smaller than comparison of tensors of the same
       size or of a scalar with a tensor of arbitrary size.
     - Unlike |matlab|'s regular ``lt()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

   * - .. function:: all(X,dim)
       .. function:: all(X,'all')
     - Checks if the entries of the tensor ``X`` are nonzero and
       performs the Boolean operation ``and`` along the dimensions
       specified in ``vecdim`` (1st form) or along every dimension
       (2nd form), producing the logical value ``true`` if all entries
       are nonzero. The result is a tensor with the same size as X,
       but with the dimensions in ``vecdim`` removed.
     - Similar syntax to |matlab|.

   * - .. function:: any(X,dim)
       .. function:: any(X,'all')
     - Checks if the entries of the tensor ``X`` are nonzero and
       performs the Boolean operation ``or`` along the dimensions
       specified in ``vecdim`` (1st form) or along every dimension
       (2nd form), producing the logical value ``true`` if at least
       one entry is nonzero. The result is a tensor with the same size
       as X, but with the dimensions in ``vecdim`` removed.
     - Similar syntax to |matlab|.

Entry-wise operations
=====================

The following functions are applied to every entry of a tensor.

.. list-table:: Entry-wise operations
   :header-rows: 1
   :stub-columns: 1
   :widths: 3 4 4
   :align: left
      
   * - Usage
     - Description
     - Notes

   * - .. function:: exp(X)
     - Exponential of tensor entries.
     - Similar syntax to |matlab|.

   * - .. function:: log(X)
     - Natural logarithm of tensor entries.
     - Similar syntax to |matlab|.

   * - .. function:: sin(X)
     - Sine of tensor entries in radians.
     - Similar syntax to |matlab|.

   * - .. function:: cos(X)
     - Cosine of tensor entries in radians.
     - Similar syntax to |matlab|.

   * - .. function:: tan(X)
     - Tangent of tensor entries in radians.
     - Similar syntax to |matlab|.

   * - .. function:: atan(X)
     - Inverse tangent in radians of tensor entries.
     - Similar syntax to |matlab|.

   * - .. function:: sqr(X)
     - Square of tensor entries.
     - Similar to ``X.^2`` or ``X .* X``.

   * - .. function:: cube(X)
     - Cube of tensor entries.
     - Similar to ``X.^3`` or ``X .* X .* X``.

   * - .. function:: X.^Y
       .. function:: power(X,Y)
     - Element-wise ``X`` raised to the power ``Y``.
     - Similar syntax to |matlab|, but |tenscalc| requires the power
       ``Y`` to be a regular numeric scalar, *not* a matrix/vector nor
       Tcalculus symbolic expression.
       
   * - .. function:: sqrt(X)
     - Square root of tensor entries.
     - Similar syntax to |matlab|.

   * - .. function:: round(X)
     - Round to nearest integer
     - Similar syntax to |matlab|, except that it does **not** support
       a second argument specifying a desired number of digits for
       rounding to a decimal.
   * - .. function:: ceil(X)
     - Round to nearest integer towards +infinity.
     - Similar syntax to |matlab|.
   * - .. function:: floor(X)
     - Round to nearest integer towards -infinity.
     - Similar syntax to |matlab|.

   * - .. function:: sign(X)
     - Signum function applied to the entries ``X``, equal to 1, 0, or
       -1, depending on whether the corresponding extry of ``X`` is
       positive, zero, or negative.
     - Similar syntax to |matlab|,
   * - .. function:: heaviside(X)
     - Step or heaviside function applied to the entries ``X``, equal
       to 1, 0.5, or 0, depending on whether the corresponding entry
       of ``X`` is positive, zero, or negative.
     - Similar syntax to |matlab|.

   * - .. function:: abs(X)
     - Absolute value of tensor entries.
     - Similar syntax to |matlab|.

   * - .. function:: relu(X)
     - Rectified linear unit activation function applied to the entries ``X``.
     - Similar to ``max(X,0)``.

   * - .. function:: srelu(X)
     - Smooth rectified linear unit activation function applied to the entries ``X``.
     - Similar to ``log(1+exp(X))`` or ``x+log(1+exp(-x))``.


   * - .. function:: normpdf(X)
     - normpdf(X) returns the pdf of the standard normal distribution
       evaluated at the entries of ``X``.
     - Similar syntax to |matlab|, except that it does **not** support
       second and third arguments specifying mean and standard
       deviation different than 0 and 1, respectively.


Linear algebra
==============

|tenscalc| supports several basic linear algebra operations, but for
operations that require some form of factorization to be performed
efficiently, such as::

     det()     logdet()     inv()     mldivide()     \     traceinv()  

|tenscalc| requires the user to explicitly select the factorization
desired::

     ldl()     lu()

This is accomplished by passing to the function the factorized matrix,
as in::

     det(lu(A))            det(ldl(A))
     logdet(lu(A))         logdet(ldl(A))
     inv(lu(A))            inv(ldl(A))
     traceinv(lu(A))       traceinv(ldl(A))
     mldivide(lu(A),Y)     mldivide(ldl(A),Y)
     lu(A)\B               ldl(A)\B

This syntax works because in |tenscalc| the factorization functions
|lu| and |ldl| return the *whole factorization as a single
entity*, which can then be passed to any function that take a
factorization as input (such as the functions listed above). This
behavior is distinct from regular |matlab| for which |lu| and
|ldl| return factorizations through multiple outputs.

|tenscalc|'s code generation makes sure that redundant computations
are not executed, therefore the following two snippets of code result
in the same computation::
  
  Tvariable A [10,10]
  Tvariable b [10]
  facA=lu(A);
  y=det(facA);
  x=facA\b;

or::

  Tvariable A [10,10]
  Tvariable b [10]
  y=det(lu(A));
  x=lu(A)\b;

Specifically note that, even though ``lu(A)`` appears twice in the
bottom snippet, the matrix ``A`` is only factored once.

Which factorization to use?
---------------------------

* The LDL factorization is faster and requires less
  memory. However, it has some limitations:

  + LDL should only be used for symmetric matrices. When used on a
    matrix that is not symmetric, all the entries above the main
    diagonal are ignored.

  + |tenscalc|'s LDL factorization only works for matrices that do not
    have zeros in the main diagonal. Structural zeros in the main
    diagonal will result in an error at code generation
    time. Non-structural zeros (i.e., zeros that cannot be determined
    at code generation time) will lead to divisions by zero at run
    time.

* The LU factorization is a little slower and requires twice as
  much memory to store both the L and U factors. However, it can be
  applied to non-symmetric matrices and matrices with zeros in the main
  diagonal.

Both the LU and the LDL factorizations, use "pychologically
lower/upper triangular matrices", i.e., matrices that are triangular
up to a permutation, with permutations selected to minimize the
fill-in for sparse matrices and reduce computation time (see
|matlab|'s documentation for ``lu()`` and ``ldl()`` with sparse
matrices).

.. list-table:: Linear Algebra operations
   :header-rows: 1
   :stub-columns: 1
   :widths: 3 4 4
   :align: left
      
   * - Usage
     - Description
     - Notes

   * - .. function:: diag(v,k)
       .. function:: diag(v)
       .. function:: diag(A)
     - When ``v`` is an ``n``-vector, returns a square matrix
       with ``n+abs(k)`` rows/columns, with the ``k``-th diagonal equal
       to ``v``. ``k=0`` corresponds to the main diagonal, ``k>0`` above the
       main diagonal and ``k<0`` below.
       
       When ``k`` ommited, it is assumed equal to (main diagonal).
       
       When ``A`` is an ``n``-by-``n`` a matrix, returns an
       ``n``-vector with the main diagonal of ``A``.
     - Similar syntax to |matlab|, except that it does **not** support
       a second argument specifying a diagonal other than the main
       diagonal, when ``A`` is a matrix.

   * - .. function:: trace(A)
     - Trace of a matrix, i.e., sum of the diagonal elements of ``A``,
       which is also the sum of the eigenvalues of ``A``.
     - Similar syntax to |matlab|.

   * - .. function:: A.'
       .. function:: A'
       .. function:: transpose(A)
       .. function:: ctranspose(A)
     - Transpose of a real-valued matrix.
     - Similar syntax to |matlab|, except that |tenscalc| does not
       support complex-valued variables and therefore ``transpose``
       and ``ctranspose`` return the same values.

   * - .. function:: lu(A)
     - LU factorization using "pychologically lower/upper triangular
       matrices", i.e., matrices that are triangular up to a
       permutation, with permutations selected to minimize the fill-in
       for sparse matrices and reduce computation time (see |matlab|'s
       documentation for lu with sparse matrices).
     - The output of this function includes the whole factorization as
       a single entity, in a format that can be passed to functions
       that require factorizations (such as mldivide, inv, det,
       logdet, traceinv), but should *not* be used by functions that
       are not expecting a factorized matrix as an input.
       
   * - .. function:: ldl(A)
     - LDL factorization using "pychologically lower triangular
       matrices", i.e., matrices that are triangular up to a
       permutation, with permutations selected to minimize the fill-in
       for sparse matrices and reduce computation time (see |matlab|'s
       documentation for lu with sparse matrices).

       All entries of ``A`` above the main diagonal are ignored and
       assumed to be equal to the one below the main diagonal,
       *without performing any test regarding of whether or not this
       is true*.
     
     - The output of this function includes the whole factorization as
       a single entity, in a format that can be passed to functions
       that require factorizations (such as mldivide, inv, det,
       logdet, traceinv), but should *not* be used by functions that
       are not expecting a factorized matrix as an input.

       
   * - .. function:: lu(A) \\ B
       .. function:: mldivide(lu(A),B)
       .. function:: ldl(A) \\ B
       .. function:: mldivide(ldl(A),B)
     - Left matrix division, which is the solution to the system of
       equations ``A*X=B`` where ``A`` must be a nonsingular square
       matrix (tensor with 2 dimensions) and ``B`` may either be a
       vector (tensor with 1 dimension) or a matrix (tensor with 2
       dimensions) with ``size(A,1)==size(A,2)==size(B,1)``.

     - |matlab| allows for non-square and possibly singular matrices
       ``A``, in which case the least-squares solution to ``A*X=B`` is
       returned. Currently, |tenscalc| requires ``A`` to be square
       and nonsingular.

       When an the LDL factorization is used, all entries of ``A``
       above the main diagonal are ignored and assumed to be equal to
       the one below the main diagonal, *without performing any test
       regarding of whether or not this is true*.
       
   * - .. function:: inv(lu(A))
       .. function:: inv(ldl(A))
     - Inverse of a square matrix.
     - Similar syntax to |matlab|.

       When an the LDL factorization is used, all entries of ``A``
       above the main diagonal are ignored and assumed to be equal to
       the one below the main diagonal, *without performing any test
       regarding of whether or not this is true*.

   * - .. function:: det(lu(A))
       .. function:: det(ldl(A))
     - Determinant of a square matrix.
     - Similar syntax to |matlab|.

       When an the LDL factorization is used, all entries of ``A``
       above the main diagonal are ignored and assumed to be equal to
       the one below the main diagonal, *without performing any test
       regarding of whether or not this is true*.

   * - .. function:: logdet(lu(A))
       .. function:: logdet(ldl(A))
     - Natural logarithm of the determinant of a square matrix.
     - Results in the same value ``log(det(A))``, but in the
       context of optimizations, it is *more efficient to use
       ``logdet(A)`` because the latter simplifies the computation of
       the derivative.*
       
       When an the LDL factorization is used, all entries of ``A``
       above the main diagonal are ignored and assumed to be equal to
       the one below the main diagonal, *without performing any test
       regarding of whether or not this is true*.

   * - .. function:: traceinv(lu(A))
       .. function:: traceinv(ldl(A))
     - Natural logarithm of the determinant of a square matrix.
     - It results in the same value ``trace(inv(A))``, but in the
       context of optimizations, it is *more efficient to use
       ``traceinv(A)`` because the latter simplifies the computation of
       the derivative.*
       
       When an the LDL factorization is used, all entries of ``A``
       above the main diagonal are ignored and assumed to be equal to
       the one below the main diagonal, *without performing any test
       regarding of whether or not this is true*.

   * - .. function:: ldl_d(A)
     - Returns the diagonal matrix in the LDL factorization.
     - All entries of ``A`` above the main diagonal are ignored and
       assumed to be equal to the one below the main diagonal,
       *without performing any test regarding of whether or not this
       is true*.

   * - .. function:: lu_d(A)
     - return the main diagonal of the U matrix in an LU
       factorization. The main diagonal of the L matrix is equal to 1.
     - All entries of ``A`` above the main diagonal are ignored and
       assumed to be equal to the one below the main diagonal,
       *without performing any test regarding of whether or not this
       is true*.

   * - .. function:: norm(x)
       .. function:: norm(x,2)
     - 2-norm of a scalar or a vector
     - Similar syntax to |matlab|, but restricted to vectors.

       The 2-norm should be avoided in optimization criteria because
       it is not differentiable at the origin. See :ref:`smoothness`
       on how to overcome this issue.

   * - .. function:: norm(x,1)
     - 1-norm of a scalar, a vector, or matrix
     - The 1-norm should be avoided in optimization criteria because
       it is not differentiable at points where the optimum often
       lies.  See :ref:`smoothness` on how to overcome this issue.

   * - .. function:: norm(x,inf)
     - infinity norm of a scalar, a vector, or matrix
     - The infinity norm should be avoided in optimization criteria
       because it is not differentiable at points where the optimum
       often lies. See :ref:`smoothness` on how to overcome this
       issue.

   * - .. function:: norm2(x)
       .. function:: norm2(x,S)
     - ``norm2(x)`` returns the sum of the square of all entries of the
       tensor ``x``, which for vectors and matrices corresponds to the square of the
       Euclidean and Frobenius norm of ``x``.

       ``norm2(x,S)`` returns the value of the quadratic form ``<x,Sx>``. This
       form is only applicable when ``x`` is a vector (tensor with 1
       dimension) and ``S`` a square matrix (tensor with 2 dimensions).
     - ``norm2(x)`` is similar to ``norm(x,2)^2`` and also to
       ``sum(x.^2,'all')``, but computes derivatives more efficiently.

   * - .. function:: norm1(x)
     - returns the sum of the absolute value of all entries of the
       tensor ``x``, which for vectors corresponds to the 1-norm of ``x``.
     - ``norm1(x)`` is similar to ``sum(abs(x),'all')``, but the former
       computes derivatives more efficiently.

       |norm1| should be avoided in optimization criteria because it
       is not differentiable at points where the optimum often
       lies. See :ref:`smoothness` on how to overcome this issue.
       
   * - .. function:: norminf(x)
     - returns the largest absolute value of all entries of the tensor
       ``x``, which for vectors corresponds to the infinity-norm of ``x``.
     - ``norminf(x)`` is similar to ``max(abs(x),[],'all')``, but the
        former computes derivatives more efficiently.

       |norminf| should be avoided in optimization criteria because it
       is not differentiable at points where the optimum often
       lies. See :ref:`smoothness` on how to overcome this issue.
       
.. warning::

   Calling `\\`, |mldivide|, |inv|, |det|, |logdet|, |traceinv|,
   |ldl_d|, or |lu_d| with a matrix that has not been factorize will
   lead to syntax errors, often not easy to directly relate to the
   missing factorization.

.. _smoothness:

Avoiding lack of smoothness
---------------------------

In general norm are not smooth:

* The 2-norm is not smooth at the origin because of the |sqrt|
* The 1-norm is not smooth at any point where one coordinate is equal to zero because of the |abs|
* The infinity-norm is not smooth at any point where two or more
  entries have the equal largest absolute values because of the |max|

For optimizations, this is particular problematic if the optimum
occurs precisely at points where the norm is not differentiable, *which
is almost always the case for the 1-norm and the infinity norm.*

However, it is generally possible to avoid this problem by introducing
auxiliary slack variables and constraints that make the optimization
smooth, without losing convexity.

* A minimization involving a 2-norm of the form:

     :math:`\min\Big\{ \text{norm}(x,2)+f(x) : x\in\mathbb{R}^n, F(x)\ge 0, G(x)= 0 \Big\}`

  can be reformulated as

     :math:`\min\Big\{ v+f(x) : v\in\mathbb{R}, x\in\mathbb{R}^n, v>0, v^2 \ge \text{norm2}(x), F(x)\ge 0, G(x)= 0 \Big\}`

* A minimization involving a 1-norm of the form:
  
     :math:`\min\Big\{ \text{norm}(x,1)+f(x) : x\in\mathbb{R}^n, F(x)\ge 0, G(x)= 0 \Big\}`

  can be reformulated as

     :math:`\min\Big\{ \text{sum}(v)+f(x) : x,v\in\mathbb{R}^n, -v\le x \le v, F(x)\ge 0, G(x)= 0 \Big\}`

* A minimization involving an infinity-norm of the form:
  
     :math:`\small\min\Big\{ \text{norm}(x,\text{inf})+f(x) : x\in\mathbb{R}^n, F(x)\ge 0, G(x)= 0 \Big\}`

  can be reformulated as

     :math:`\small\min\Big\{ v+f(x) : v\in\mathbb{R},\; x\in\mathbb{R}^n, -v\le x \le v, F(x)\ge 0, G(x)= 0 \Big\}`

   
Calculus - differentiation
==========================

The following function computes the partial derivatives of a symbolic
expression with respect to the entries of a tensor-valued symbolic
variable:

.. function:: gradient(f,x)

   :param f: tensor-valued expression to be differentiated
   :type f: |tcalculus| symbolic tensor
   :param x: variable (created with |tvariable|) with respect to which
             the derivatives will be taken
   :type x: |tcalculus| tensor-values symbolic variable

When

* ``f`` is a tensor with size ``[n1,n2,...,nN]``
* ``x`` is a tensor-valued variable (created with |tvariable|) with size
  ``[m1,m2,...,mM]``

then ``g=gradient(f,x)`` results in a tensor with size
``[n1,n2,...,nN,m1,m2,...,mM]`` with

.. math::

  g(i1,i2,...,iN,j1,j2,...,jM)=\frac{d\,f(i1,i2,...,iN)}{d x(j1,j2,...,jM)}

For example, if ``f`` is a scalar (with size ``[]``) and ``x`` an
``n``-vector (with size ``[n]``), then ``g=gradient(f,x)`` results in
an ``n``-vector (with size ``[n]``) with the usual gradient:

.. math::

  g(i)=\frac{d\,f}{d x(i)}

If we then compute ``h=gradient(g,x)``, we obtain an ``n``-by-``n``
matrix (with size ``[n,n]``) with the Hessian matrix:

.. math::

  h(i,j)=\frac{d\,g(i)}{d x(j)} = \frac{d^2\,f}{d x(i)\,d x(j)}

The computation of first and second derivatives can be
streamlined using the following function:

.. function:: hessian(f,x[,y])

   :param f: tensor-valued expression to be differentiated
   :type f: |tcalculus| symbolic tensor
   :param x: variable (created with |tvariable|) with respect to which
             the 1st derivatives will be taken
   :type x: |tcalculus| tensor-values symbolic variable
   :param y: variable (created with |tvariable|) with respect to which
             the 2nd derivatives will be taken (optinal, when omitted
             the 2nd derivatives will also be taken with respect to x)
   :type y: |tcalculus| tensor-values symbolic variable

When

* ``f`` is a tensor with size ``[n1,n2,...,nN]``
* ``x`` a tensor-valued variable (created with |tvariable|) with size
  ``[m1,m2,...,mM]``
* ``y`` a tensor-valued variable (created with |tvariable|) with size
  ``[l1,l2,...,lL]``

then ``[h,g]=hessian(f,x,y)`` results in

* a tensor ``g`` with size
  ``[n1,n2,...,nN,m1,m2,...,mM]`` with

  .. math::

     g(i1,i2,...,iN,j1,j2,...,jM)=\frac{d\,f(i1,i2,...,iN)}{d x(j1,j2,...,jM)}

* a tensor ``h`` with size
  ``[n1,n2,...,nN,m1,m2,...,mM,l1,l2,...,lL]`` with

  .. math::

     h(i1,i2,...,iN,j1,j2,...,jM,k1,k2,...,kL)=\frac{d\,f(i1,i2,...,iN)}{d x(j1,j2,...,jM)d y(j1,j2,...,jM)}

In practice, ``[h,g]=hessian(f,x,y)`` is equivalent to::

  g=gradient(f,x);
  h=gradient(g,y);
  

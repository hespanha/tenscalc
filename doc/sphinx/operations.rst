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
  
Boolean operation
=================

|tenscalc| uses boolean operations mostly to specify optimization
constraints.

.. list-table:: Boolean operations
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
       or of a scalar with a tensor of arbitrary size
     - Unlike |matlab|'s regular ``eq()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

   * - .. function:: X>=Y
       .. function:: ge(X,Y)
     - Entry-by-entry greater than or equal to comparison of tensors
       of the same size or of a scalar with a tensor of arbitrary size
     - Unlike |matlab|'s regular ``ge()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

       From the perspective of a constrained optimization numerical
       solver, due to finite numerical precision, X>=Y and X>Y
       represent the same constraint.
       
   * - .. function:: X>Y
       .. function:: gt(X,Y)
     - Entry-by-entry greater than comparison of tensors of the same
       size or of a scalar with a tensor of arbitrary size
     - Unlike |matlab|'s regular ``gt()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

   * - .. function:: X<=Y
       .. function:: le(X,Y)
     - Entry-by-entry smaller than or equal to comparison of tensors
       of the same size or of a scalar with a tensor of arbitrary size
     - Unlike |matlab|'s regular ``le()``, expansion upon singleton
       dimensions is not performed automatically to match the tensors'
       sizes.

       From the perspective of a constrained optimization numerical
       solver, due to finite numerical precision, X<=Y and X<Y
       represent the same constraint.

   * - .. function:: X<Y
       .. function:: lt(X,Y)
     - Entry-by-entry smaller than comparison of tensors of the same
       size or of a scalar with a tensor of arbitrary size
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
     - Similar syntax to |matlab|

   * - .. function:: any(X,dim)
       .. function:: any(X,'all')
     - Checks if the entries of the tensor ``X`` are nonzero and
       performs the Boolean operation ``or`` along the dimensions
       specified in ``vecdim`` (1st form) or along every dimension
       (2nd form), producing the logical value ``true`` if at least
       one entry is nonzero. The result is a tensor with the same size
       as X, but with the dimensions in ``vecdim`` removed.
     - Similar syntax to |matlab|

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
     - Similar syntax to |matlab|

   * - .. function:: log(X)
     - Natural logarithm of tensor entries.
     - Similar syntax to |matlab|

   * - .. function:: sin(X)
     - Sine of tensor entries in radians.
     - Similar syntax to |matlab|

   * - .. function:: cos(X)
     - Cosine of tensor entries in radians.
     - Similar syntax to |matlab|

   * - .. function:: tan(X)
     - Tangent of tensor entries in radians.
     - Similar syntax to |matlab|

   * - .. function:: atan(X)
     - Inverse tangent in radians of tensor entries.
     - Similar syntax to |matlab|

   * - .. function:: reciprocal(X)
     - Reciprocal of tensor entries.
     - Similar to ``1./X``

   * - .. function:: sqr(X)
     - Square of tensor entries
     - Similar to ``X.^2`` or ``X .* X``

   * - .. function:: cube(X)
     - Cube of tensor entries
     - Similar to ``X.^3`` or ``X .* X .* X``

   * - .. function:: sqrt(X)
     - Square root of tensor entries
     - Similar syntax to |matlab|

   * - .. function:: round(X)
     - Round to nearest integer
     - Similar syntax to |matlab|, except that it does support a
       second argument specifying a desired number of digits for
       rounding to a decimal.
   * - .. function:: ceil(X)
     - Round to nearest integer towards +infinity
     - Similar syntax to |matlab|
   * - .. function:: floor(X)
     - Round to nearest integer towards -infinity
     - Similar syntax to |matlab|

   * - .. function:: abs(X)
     - INCOMPLETE
     - Similar syntax to |matlab|
   * - .. function:: relu(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: srelu(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: heaviside(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: normpdf(X)
     - INCOMPLETE
     - INCOMPLETE


Linear algebra
==============

.. list-table:: Linear Algebra operations
   :header-rows: 1
   :stub-columns: 1
   :widths: 3 4 4
   :align: left
      
   * - Usage
     - Description
     - Notes

   * - .. function:: A \\ B
       .. function:: mldivide(A,B)
     - Left matrix division, which is the solution to the system of
       equations ``A*X=B`` where ``A`` must be a nonsingular square
       matrix (tensor with 2 dimensions) and ``B`` may either be a
       vector (tensor with 1 dimension) or a matrix (tensor with 2
       dimensions) with ``size(A,1)==size(A,2)==size(B,1)``.

       ALSO TALK ABOUT FACTORIZATION
       
     - |matlab| allows for non-square and possibly singular matrices
       ``A``, in which case the least-squares solution to ``A*X=B`` is
       returned. Currently, |tenscalc| requires ``A`` to be square
       and nonsingular.

   * - .. function:: diag(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: trace(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: transpose(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: ctranspose(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: inv(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: det(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: chol(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: ldl(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: ldl_l(X)
     - INCOMPLETE
     - INCOMPLETE
     
   * - .. function:: ldl_d(X)
     - INCOMPLETE
     - INCOMPLETE
     
   * - .. function:: lu(X)
     - INCOMPLETE
     - INCOMPLETE
     
   * - .. function:: logdet(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: traceinv(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: norm(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: norm1(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: norm2(X)
     - INCOMPLETE
     - INCOMPLETE

   * - .. function:: norminf(X)
     - INCOMPLETE
     - INCOMPLETE

     
Calculus - differentiation
==========================


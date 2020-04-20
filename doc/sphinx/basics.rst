.. include:: defs.rst

========
 Basics
========

|tenscalc| computations and optimizations are defined using symbolic
expressions supported by the class |tcalculus|, named after *tensor
calculus*\ .

.. _tensors-section:

What are Tensors?
=================

In |tenscalc|, tensors are multi-dimensional real-valued arrays over
which one performs numerical operations. The |size| of a tensor is a
vector that specifies how many entries the tensor has in each
dimension. The following table shows common tensors and their sizes:

.. list-table:: Common tensor types
   :header-rows: 1
		     
   * - name
     - number of dimensions
     - |tenscalc|'s size
     - mathematical set
       
   * - scalar
     - 0
     - ``[]``
     - :math:`\mathbb{R}`

   * - :math:`n`-vector
     - 1
     - ``[n]``
     - :math:`\mathbb{R^n}`

   * - :math:`n`-by-:math:`m` matrix
     - 2
     - ``[n,m]``
     - :math:`\mathbb{R^{n\times m}}`

   * - tensor with :math:`d` of dimensions
     - ``d``
     - ``[n1,n2,...,nd]``
     - :math:`\mathbb{R^{n_1\times n_2\times \cdots\times n_d}}`

.. warning::
       
   It is common (at least within |matlab|) to make no distinction
   between a scalar, a 1-vector, and a 1-by-1 matrix. In |matlab| it
   is also common not to distinguish between an :math:`n`-vector and
   an :math:`n`-by-1 matrix. However, |tenscalc| does make a
   distinction between all these and does not allow, e.g., adding a
   2-vector with a 2-by-1 matrix, without an explicit |reshape|
   operation.

|tenscalc| is slightly more permissive with operations between scalars
and tensors of larger sizes in that it does allow summations between a
scalar and a tensor of any size, with the understanding that the
scalar is added to *every entry* of the tensor, which is consistent
with |matlab|'s syntax.

See :ref:`size-conversion-section` regarding the rules for automatic
conversion of sizes between |matlab| matrices and |tenscalc| tensors.


Building |tcalculus| symbolic expressions
===========================================

|tenscalc| symbolic expressions are supported by the class

.. class:: Tcalculus

but you never really need to call this class directly. Instead,
symbolic expressions are built like regular |matlab| expressions using
many of the usual |matlab| operators and functions, which have been
overloaded to operate on |tenscalc| symbolic expressions.

.. warning:: To speed up operations, the class |tcalculus| keeps some
   information in a few auxiliary global variables. Before defining a
   new computation or optimization it is a good idea to clear these
   variables with

     clear all

   If this is not possible, one can use the following method to just
   clear |tcalculus|'s auxiliary global variables
     
   .. method:: Tcalculus.clear

   However, after using this method, any instances of the |tcalculus|
   class that remains in memory becomes invalid and will lead to
   errors if used in a subsequent computations or optimizations
   **without necessarily resulting in a syntax error**. It is thus
   strongly advised to use ``clear all`` rather than
   ``Tcalculus.clear``.
   

	       
.. _Tvariable-section:

Symbolic variables
------------------

Most symbolic expressions start with symbolic variables that represent
tensors that will only be assigned numeric values at code-execution
time. Symbolic variables are created using

.. function:: Tvariable name size

   :param name: Name of the variable to be created
   :type name: character string		     
   :param size: Size of the tensor as a vector of integers.
		When omitted, the empty size ``[]`` is assumed, which
                corresponds to a scalar.
   :type size: vector of integers
		       
The variable is created in the caller's workspace::
  
   >> clear all
   >> Tvariable a [3,4];
   >> whos
      Name      Size            Bytes  Class        Attributes
      a         3x4                39  Tcalculus
      ans       3x4                39  Tcalculus

|tvariable| can also be used as a regular |matlab| function using the syntax:

.. function:: Tvariable(name,size)

   :returns: |tcalculus| symbolic expression holding the symbolic
             variable
	    
In this case, |tvariable| returns a symbolic expression that contains
the variable. Note that the name of the symbolic variable is still
given by ``name`` regardless of the name of the variable that you use
to store it in::
      
  >> clear all
  >> y=Tvariable('a',[3,4]);
  >> whos
     Name      Size            Bytes  Class        Attributes
     y         3x4                39  Tcalculus
  >> disp(y)
  variable1   = variable('a',) : [3,4]


.. _size-conversion-section:

Using |matlab| matrices in |tenscalc| symbolic expressions
----------------------------------------------------------

Often symbolic expressions involve regular |matlab| matrices, as in adding a symbolic variable created with |tvariable| with a regular |matlab| matrix::

  >> Tvariable a [2,2]
  >> b=a+[2,1;3,4]

The following rules are used to convert |matlab| matrices to |tenscalc| symbolic expressions. In essence, the conversion is completely transparent *unless the matrix is 1-by-1 or n-by-1.*

.. _size-conversion-table:

.. list-table:: Rules for automatic convertion between |matlab| matrices and |tenscalc| tensors
   :header-rows: 1
		     
   * - |matlab| variable
     - |matlab|'s size
     - |tenscalc| variable
     - |tenscalc|'s size

   * - 1-by-1 matrix
     - ``[1,1]``
     - scalar
     - ``[]``

   * - ``n``-by-1 column matrix
     - ``[n,1]``
     - ``n``-vector
     - ``[n]``

   * - 1-by-``n`` row matrix
     - ``[1,n]``
     - 1-by-``n`` tensor
     - ``[1,n]``

   * - ``m``-by-``n`` matrix
     - ``[m,n]``
     - ``m``-by-``n`` tensor
     - ``[m,n]``

   * - ``n1``-by-``n2``- ... -by-``nd`` matrix
     - ``[n1,n2,...,nd]``
     - ``n1``-by-``n2``- ... -by-``nd`` tensor
     - ``[n1,n2,...,nd]``

|tenscalc| also provides a few functions that can be used to directly
create symbolic tensors that can be useful when either the rules above
do not have the desired effect (e.g., in the somewhat unlikely case
wants to create a ``1-by-1`` tensor, rather than a scalar) or when it
is more efficient to directly create the symbolic tensor.

.. list-table:: Functions to create (constant) tensors
   :header-rows: 1
   :stub-columns: 1
   :widths: 3 4 4
   :align: left
      
   * - Usage
     - Description
     - Notes

   * - .. function:: Tzeros([n1,...,nd])
       .. function:: Tzeros(n1,...,nd)
     - Returns a tensor with all entries equal to 0. The tensor size
       can be provided as a single vector with the number of entries
       in all directions, or as multiple input parameters, one per
       dimension.

     - This function can be used interchangeably with the regular
       |matlab| function ``zeros``, as long as the the size conversion
       rules in :ref:`size-conversion-table` are appropriate.
       
   * - .. function:: Tones([n1,...,nd])
       .. function:: Tones(n1,...,nd)
     - Returns a tensor with all entries equal to 1. The tensor size
       is specified as in |tzeros|.
     
     - This function can be used interchangeably with the regular
       |matlab| function ``ones``, as long as the the size conversion
       rules in :ref:`size-conversion-table` are appropriate.
     
   * - .. function:: Teye([n1,...,nk,n1,...,nk])
       .. function:: Teye(n1,...,nk,n1,...,nk)
     - Returns an "identity" tensor with all entries equal to zero,
       except for the entries with the indice 1 equal to the indice
       k+1, the indice 2 equal to the indice k+2, ..., as in

       ``(i1,i2,...,ik,i1,i2,...,ik)``

       These entries are all equal to 1.


     - For ``k=1``, this function can be used interchangeably with the
       regular |matlab| function ``eye``, but is convenient to
       generate "identity" tensors with a larger number of dimensions.
 
   * - .. function:: Tconstant(mat,size)
     - Returns a tensor with entries given by the matrix ``mat`` and
       size specified by ``size`` as in |tzeros|.

     - This function is typically used to override the size conversion
       rules in :ref:`size-conversion-table`.
     
       
.. warning::

   When called with a single argument |tzeros|, |tones|, and |teye|
   differ from their |matlab| counterparts: |tzeros| and |tones|
   return vectors (i.e., tensors with 1 dimension), rather than square
   matrices, and |teye| generates an error since it has no
   counter-part for vectors.
     
     
Accessing information about the size of symbolic variables
==========================================================

.. list-table:: Functions to access information about symbolic expressions
   :header-rows: 1
   :stub-columns: 1
   :widths: 1 3 3
   :align: left

   * - Usage
     - Description
     - Notes

   * - .. function:: size(X)
       .. function:: [n1,n2,...,nd]=size(X)
       .. function:: size(X,dim)
     - Size of a tensor
     - Similar syntax to |matlab|
       
   * - .. function:: ndims(X)
     - Number of dimensions in a tensor
     - Similar syntax to |matlab|.  Note that a scalar (which has empty
       ``size = []``) always has 0 dimensions.
       
   * - .. function:: numel(X)
     - Number of elements in a tensor
     - Similar syntax to |matlab|.  Note that a scalar (which has empty
       ``size = []``) always has 1 elements.
       
   * - .. function:: length(X)
     - Length of a tensor
     - Similar syntax to |matlab|.  Note that a scalar (which has empty
       ``size = []``) always has 1 elements.
       
   * - .. function:: isempty(X)
     - True for an empty tensor
     - Similar syntax to |matlab|.  Note that a scalar (which has empty
       ``size = []``) always has 1 elements so it is never empty.
       
Indexing and resizing symbolic variables
========================================

|tenscalc| uses the same syntax as |matlab| to index tensors:

* ``X(i,j,k,...)`` returns a subtensor formed by the elements of ``X``
  with subscripts vectors ``i,j,k,...``

  The resulting tensor has the same number of dimensions as ``X``,
  with lengths along each dimension given by ``length(i), length(j),
  length(k),...``

* A colon ``:`` can be used as a subscript to indicate all subscripts
  on that particular dimension.

* The keyword ``end`` can be used within an indexing expression to
  denote the last index. Specifically, ``end = size(X,k)`` when used
  as part of the ``k``\ th index.

* ``subsref(X,S)`` with a subscript reference structure ``S`` behaves
  as in regular |matlab|, with he caveat that entries of the tensor
  must always be indexed using subscripts the type ``()``.  However,
  this does not preclude the construction of cells of |tcalculus|
  objects.


.. list-table:: Functions reshape and resize symbolic variables
   :header-rows: 1
   :stub-columns: 1
   :widths: 2 2 4
   :align: left

   * - Usage
     - Description
     - Notes

   * - .. function:: reshape(X,size)
       .. function:: reshape(X,n1,...,nd)
     - Reshape array 
     - Similar syntax to |matlab|, except that it does not support using
       ``[]`` as one of the dimensions.


   * - .. function:: repmat(X,size)
       .. function:: repmat(X,n1,...,nd)
     - Replicate and tile tensor
     - Similar syntax to |matlab|, with the understanding that the
       number of dimensions of ``X`` must be strictly preserved. Often
       |repmat| needs to be combined with |reshape| to obtain the
       desired effect. For example to replicate twice a 3-vector ``X``
       to create a 3-by-2 matrix by placing the copies of ``X`` side
       by side, one needs::
	 
	 Y=repmat(reshape(X,3,1),1,2);

       or, to create a 2-by-3 matrix by placing the copies of ``X``
       one on top of the other, one needs::
	 
         Y=repmat(reshape(X,1,3),2,1);

   * - .. function:: cat(dim,A,B,...)
     - Concatenate tensors along the dimension ``dim``
     - Similar syntax to |matlab|, with the understanding that the
       number of dimensions of the input tensors ``A,B,...`` must
       match for all but dimensions but ``dim``. |tenscalc| will try
       to add singleton dimensions to make |cat| succeed as best as it
       can. For example in the following code, the scalar variable
       ``a`` with size ``[]`` is augmented to the size ``[1,1]`` so
       that concatenation to a 2-by-3 matrix is possible::

         >> Tvariable a [];
         >> b=[a,a,a;a,a,a];
         >> size(b)
         ans =
	     2     3
	 
       and in the following code, the vector variable ``b`` with size
       ``[3]`` is augmented to the size ``[3,1]`` so that
       concatenation to a 3-by-2 matrix is possible::
       
         >> Tvariable b [3];
	 >> c=[b,b];
	 >> size(c)
         ans =
	     3     2
	 
   * - .. function:: vertcat(A,B,...)
     - Concatenate tensors along the 1st dimension, which is equivalent to
       ``cat(1,A,B,...)` and also to ``[A;B;...]``
     - Similar syntax to |matlab|, with the same caveats as |cat|.
       
   * - .. function:: horzcat(dim,A,B,...)
     - Concatenate tensors along the 2nd dimension, which is equivalent to
       ``cat(2,A,B,...)`` and also to ``[A,B,...]``
     - Similar syntax to |matlab|, with the same caveats as |cat|.


   * - .. function:: vec2tensor(X,sz,subs)
       .. function:: vec2tensor(X,sz,subs,dim)
     - Expands a vector to a sparse tensor
     - In the 1st form
         :param X: |tcalculus| ``n``-vector (tensor with size ``[n]``)
	 :param sz: vector with ``d`` integers
	 :param subs: ``n``-by-``d`` matrix of subscripts

       and returns a |tcalculus| tensor ``Y`` with size ``sz``, with
       the nonzero entries taken from ``X``, with
       ``Y(subs(i,:))=X(i)`` for ``i=1:n``

       In the 2nd form
         :param X: |tcalculus| tensor ``X`` with ``size(X,dim)=n``
	 :param sz: vector with ``d`` integers
	 :param subs: ``n``-by-``d`` matrix subscripts

       and returns a |tcalculus| tensor ``Y`` with size similar to
       that of ``X``, but the ``dim`` dimension expanded to the sizes
       in ``sz``, and the nonzero entries taken from ``X``, with
       ``Y(...,subs(i,:),...)=X(...,i,...)`` for ``i=1:n``, where the
       ``...`` denote indices of the dimensions before and after
       ``dim``

       This function is typically used to create sparse tensors from
       the entries of a (typically full) vector. See :ref:`structured-matrices`
     
.. warning::
       
   Unlike |matlab|, |tenscalc| does *not* allow for subscripted assignments, such as

     ``A(1,2)=5``

   This functionality needs to be achieved with the concatenation
   operations |cat|, |vertcat|, |horzcat|.

.. _structured-matrices:

Creating structured matrices
----------------------------

The function |vec2tensor| is very useful to create structured matrices from vectors, to be used as optimization variables

* Diagonal matrix::

    % Creates an NxN diagonal matrix
    Tvariable v [N];
    A=vec2tensor(v,[N,N],[1:N;1:N]');

* Lower triangular matrix::
  
    % Creates an NxN lower triangular matrix
    [i,j]=find(ones(N));
    k=find(i>=j);
    Tvariable v length(k);
    A=vec2tensor(v,[N,N],[i(k),j(k)]);

* Symmetric matrix::

    % Creates an NxN lower triangular matrix
    [i,j]=find(ones(N));
    kl=find(i>j);
    k0=find(i==j);
    Tvariable v length(k0)+length(kl);
    A=vec2tensor([v;v(1:length(kl))],[N,N],[i(kl),j(kl);i(k0),j(k0);j(kl),i(kl)]);

In all these examples, one would set ``v`` to be an optimization
variable, that implicitly represents the structured matrix.

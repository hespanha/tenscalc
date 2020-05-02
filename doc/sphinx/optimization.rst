.. include:: defs.rst

========================
Constrained optimization
========================

|tenscalc| can be used to generate code to solve optimizations or Nash
equilibria. All solvers use a primal-dual interior-point method and
share several common parameters.

Types of code
=============

|C| code

  |C| code is self-contained and library free (aside from the |C|
  standard library). The |C| code is encapsulated into a dynamic library
  that can be linked to |matlab| or used independently.

  This type of code is extremely efficient and fast for small to
  medium-size problems; typically up to a few thousands of variables
  and constraints. In this domain |C| code is typically **10-100 times
  faster.**

  In addition to the |C| code, a |matlab| class is also generated to
  set parameter values and call the solver from within |matlab|.

  The scripts that generate |C| code start with the prefix ``cmex2``.

|matlab|

  |matlab| code is integrated into a |matlab| class that is used to set
  parameter values and call the solver.

  The scripts that generate |matlab| code start with the prefix ``class2``.

.. note::

   The |matlab| classes that call the |C| and the |matlab| solvers
   appear indistinguishable to the user, aside from the speed of
   execution. 
  

Minimization
============

|tenscalc| can generate optimized code to solve constrained
optimizations of the form

.. math::
   x^* \in \arg\min_x\Big\{ f(x) : F(x)\ge 0, G(x)=0 \Big\}

where the variable ``x`` can include multiple tensors and the equality
and inequality constraints can be expressed by equalities and
inequalities involving multiple tensors.

The following two scripts are used to generate code to solve this type
of constrained minimization

.. function:: cmex2optimizeCS(parameter1,value1,parameter2,values2, ...)
.. function:: class2optimizeCS(parameter1,value1,parameter2,values2, ...)

   :param parameter1: parameter to set
   :type parameter1: string
   :param value1: parameter to set
   :type value1: type depends on the parameter
   :param parameter2: parameter to set
   :type parameter2: string
   :param value2: parameter to set, ...
   :type value2: type depends on the parameter
   :returns: name of the |matlab| class created
   :rtype: string
		       
Both scripts take several parameters (many of them optional) in the
form of pairs, consisting of a parameter name and its value, as in::

  classname=cmex2optimizeCS('classname','myoptimization', ...
                            'objective',norm2(x-y), ...
                            'optimizationVariables', { x, y }, ...
			    'constraints', { x>=-1, x<=1, y>=10 }, ...
			    'outputExpressions', { x, y });

Alternative, parameters to the scripts can be passed using a structure, as in::

   opt.classname='myoptimization';
   opt.objective=norm2(x-y);
   opt.optimizationVariables={ x, y };
   opt.constraints={ x>=-1, x<=1, y>=10 };
   opt.outputExpressions={ x, y };
   classname=cmex2optimizeCS(opt);

.. note:: This type of parameter passing and validation is enabled by
   the |funpartools| toolbox, which also enables several other
   advanced features. See :ref:`pedigrees`.

The function |cmex2optimizeCS| generates |C| code, whereas
|class2optimizeCS| generates |matlab| code, but both functions take
the same set of parameters and generate |matlab| classes that are
indistinguishable to the user.

The following table list the most commonly used parameters for
|cmex2optimizeCS| and |class2optimizeCS|. For the full set of
parameters, use::

   cmex2optimizeCS help
   class2optimizeCS help



.. list-table:: Selected parameters for |cmex2optimizeCS| and
                |class2optimizeCS|. For the full set of parameters use
                ``cmex2optimizeCS help`` or ``class2optimizeCS help``
   :header-rows: 1
   :stub-columns: 1
   :widths: 1 1 3
   :align: left
      
   * - Parameter
     - Allowed values
     - Description

   * - 'optimizationVariables'
     - cell-array of |tcalculus| tensor variables created using |tvariable|
     - Variables to be optimized.
       
   * - 'objective'
     - scalar |tcalculus| tensor
     - Criterion to optimize.
   
   * - 'constraints'
     - cell-array of |tcalculus| tensors, each involving one of the
       following operations ``==``, ``>=``, ``<=``, ``>``, ``<``.
     - Optimization constraints.

       .. warning::

	  For constraint satisfaction, there is no difference between
	  ``<=`` and ``<`` or between ``>=`` and ``>``.

   * - 'outputExpressions'
     - cell-array of |tcalculus| tensors
     - Expressions (typically involving the optimization variables)
       that the solver should return upon termination.

       .. warning::

	  |matlab| does not support sparse arrays with more than 2
          dimensions, therefore to include in ``outputExpressions``
          tensors with more than 2 dimensions one must make them |full|, as in::

	    cmex2optimizeCS('classname','myoptimization', ...
                            'objective',norm2(x-x0), ...
                            'optimizationVariables', { x }, ...
			    'parameters', {x0, ...
	 	  	    'constraints', { x>=-1, x<=1 }, ...
    			    'outputExpressions', { full(x) });
	  
       See :ref:`special-variables-optimization`.

   * - 'parameters'
     - cell-array of |tcalculus| tensor variables created using |tvariable|
     - Parameters that the optimization criterion and constraints may
       depends upon; to be provided to the solver, but not optimized.
       
   * - 'sensitivityVariables'
     - cell-array of |tcalculus| tensor variables created using |tvariable|
     - Optimization variables with respect to which we want to compute cost
       sensitivity.

       See :ref:`sensitivity-variables`.

   * - 'classname'
     - string
     - Name of the class to be created.  A |matlab| class will be
       created with this name plus a ``.m`` extension.

       See :ref:`solver-class`.

       .. warning:: This parameter should not be used with pedigrees
                    are enabled using ``'pedigreeClass'``.  See
                    :ref:`pedigrees`

       
   * - 'folder'
     - string
     - Path to the folder where the class and cmex files will be
       created.  The folder will be created if it does not exist and
       it will be added to the begining of the path if not there
       already.
       
   * - 'simulinkLibrary'
     - string
     - Name of a simulink library to be created with Simulink blocks that can be used
       to call the different solver functions.
       The blocks are created with direct feedthrough.
       
       No library and no simulink blocks are created if simulinkLibrary is an empty string.
       
   * - 'gradTolerance'
     - positive number, default 1e-4
     - Maximum norm for the gradient below which the first order
       optimality conditions assumed to by met.
       
   * - 'equalTolerance'
     - positive number, default 1e-4
     - Maximum norm for the vector of equality constraints below which
       the equalities are assumed to hold.
       
   * - 'desiredDualityGap'
     - positive number, default 1e-5
     - Value for the duality gap that triggers the end of the
       constrained optimization. The overall optimization terminates
       at the end of the first Newton step for which the duality gap
       becomes smaller than this value.
       
   * - 'solverVerboseLevel'
     - integer
     - Level of verbose for the solver outputs:
       
       * 0 - the solver does not produce any output and does not report timing information
       * 1 - the solver does not produce any output but reports timing information
       * 2 - the solver only report a summary of the solution status when the optimization terminates
       * 3 - the solver reports a summary of the solution status at each iteration
       * >3 - the solver produces several (somewhat unreadable) outputs at each iteration step
       
   * - 'addEye2Hessian'
     - nonnegative real, default 1e-9
     - Add to the Hessian matrix appropriate identity matrices scaled by this constant.
       
       A larger value for addEye2Hessian has two main effects:
       
       1) Improves the numerical conditioning of the system of equations that
	  finds the Newton search direction.
       2) Moves the search direction towards the gradient descent of
	  the Lagragian (and away from the Newton direction).
	  
       Both effects improve the robustness of the solver, but this is typically
       achieved at the expense of slower convergence.
       
       For convex problems, one typically chooses addEye2Hessian equal to the
       square root of the machine precision.
       
       For non-convex problems, one can try to increase this parameter when
       the Newton direction actually causes and increase of the Lagragian.
       
   * - 'muFactorAggressive'
     - real in the interval (0,1), default .3333
     - Multiplicative factor used to update the barrier parameter mu.
       This value is used when there is good progress along the Newton
       direction.  Nice convex problems can take as low as 1/100, but
       poorly conditioned problems may require as high as 1/3.  This
       parameter is only used when ``skipAffine=true``.

   * - 'muFactorConservative'
     - real in the interval (0,1), default .75
     - Multiplicative factor used to update the barrier parameter mu
       (must be smaller than 1).  This value is used when there is
       poor or no progress along the Newton direction. A value not
       much smaller than one is preferable.

   * - 'skipAffine'
     - [``false``, ``true``], default ``true``
     - When false the barrier parameter ``mu`` is updated based on how
       much progress can be achieved in the search direction obtained
       with ``mu=0``.  This is known as the 'affine search direction
       step'.  This step (obtained by setting to ``skipAffine=false``)
       can significantly speed up convergence by rapidly decreasing
       the barrier parameter.  However, it can be fragile for tough
       non-convex problems.

   * - 'delta'
     - [2,3], default 3
     - Delta parameter used to determine mu based on the affine search
       direction step: Set ``delta=3`` for well behaved problems (for
       an aggressive convergence) and ``delta=2`` in poorly
       conditioned problems (for a more robust behavior).  This
       parameter is only used when ``skipAffine=false``.

   * - 'alphaMin'
     - real in the interval (0,1], default 1e-7
     - Minimum value for the scalar gain in Newton's method line
       search, below which a search direction is declared to have
       failed.
       
   * - 'alphaMax'
     - real in the interval (0,1], default 1
     - Maximum value for the scalar gain in Newton's method line
       search.  Should only be set lower than 1 for very poorly scaled
       problems.
       
   * - 'useLDL'
     - [``false``, ``true``], default ``true``
     - When ``true`` the search directions are computed using an
       LDL instead of an LU factorization.

       In general, the LDL factorization leads to faster code.
       However, the current implementation is restricted to a pure
       diagonal matrix (no 2x2 blocks in the D factor) so it may fail
       with the message 'ldl needs pivoting'. If this happens either
       set ``useLDL=false`` or use a nonzero value for
       ``addEye2Hessian``.

   * - 'smallerNewtonMatrix'
     - [``false``, ``true``], default ``false``
     - When ``true`` the matrix that needs to be inverted to compute a
       Newton step is reduced by first eliminating the dual variables
       associated with inequality constraints.  However, often the
       smaller matrix is not as sparse so the computation may actually
       increase.
       
   * - 'compilerOptimization'
     - [``'-O0'``, ``'-O1'``, ``'-O2'``, ``'-O3'``, ``'-Ofast``' ], default ``'-O0'``
     - Optimization parameters passed to the |C| compiler.

       :-O1: often generates the fastest code, whereas
       :-O0: compiles the fastest

       This parameter is only used for |C|-code solvers.
       
   * - 'minInstructions4loop'
     - integer, default 50
     - Minimum number of similar instruction that will be execute as
       part of a ``for(;;)`` loop, rather than being executed as
       independent |C| commands.  When equal to ``inf``,
       instructions will never be grouped into foor loops.

       This parameter is only used for |C|-code solvers.
       
   * - 'maxInstructionsPerFunction'
     - integer, default 100
     - Maximum number of instructions to be included into a single
       function.  When equal to ``inf``, there is no limit on the size
       of a single function.

       Large values of ``maxInstructionsPerFunction`` and therefore
       large functions give more opportunities for compiler
       optimization, but can resul in very slow compilation
       (especially with compiler optimization turned on).

       This parameter is only used for |C|-code solvers.
       
   * - 'absolutePath'
     - [``false``, ``true``], default ``true``
     - When ``true`` the the cmex functions use an absolute path to
       open the dynamic library, which means that the dynamic library
       cannot be moved away from the folder where it was created.

       When ``false`` no path information about the dynamic library is
       included in the cmex function, which must then rely on the
       OS-specific method used to find dynamic libraries. See
       documentation of 'dlopen' for linux and OSX or 'LoadLibrary'
       for Microsoft Windows.

       This parameter is only used for |C|-code solvers.
       
   * - 'pedigreeClass'
     - string
     - When nonempty, the function outputs are saved to a set of files.
       All files in the set will be characterized by a 'pedigree',
       which decribes all the input parameters that were used in the script.
       This variable contains the name of the file class and may include a path.
     
       This is supported by the |funpartools| toolbox. See ``help
       createPedigree``.

       See :ref:`pedigrees`

   * - 'executeScript'
     - [``'yes'``, ``'no'``, ``'asneeded'``], default ``'yes'``
     - Determines whether or not the body of the function should be executed:
       
       :yes: the function body should always be executed.
       :no: the function body should never be executed and therefore
            the function returns after processing all the input
            parameters.
       :asneeded: if a pedigree file exists that match all the input
		  parameters (as well as all the parameters of all
		  'upstream' functions) the function body is not
		  executed, otherwise it is execute.
       
       See :ref:`pedigrees`.

.. warning::

   |cmex2optimizeCS| and |class2optimizeCS| only solve for first-order
   optimality conditions so they may produce either a minimium or a
   maximum (or a saddle-point). It is up to the user to restrict the
   search domain to make sure that the desired optimum was found, or
   test if that was the case after the fact. The ``Hess_`` output can
   help in that regard.

.. _`special-variables-optimization`:

Special variables to include in ``'outputExpressions'``
-------------------------------------------------------

The following |tcalculus| variables are assigned special values and
can be using in ``outputExpressions``:
       
:``lambda1_``, ``lambda2_``, ...: Lagrangian multipliers associated
                                  with the inequalities constraints
                                  (in the order that they appear and
                                  with the same size as the
                                  corresponding constraints)
:``nu1_``, ``nu2_``, ...: Lagrangian multipliers associated with the
                          equality constraints (in the order that they
                          appear and with the same size as the
                          corresponding constraints)
:``Hess_``: Hessian matrix used by the (last) Newton step to update
            the primal variables (not including addEye2Hessian).
:``dHess_``: D factor in the LDL factorization of the Hessian matrix
             used by the (last) Newton step to update the primal
             variables (including addEye2Hessian, unlike ``Hess_``).
:``Grad_``: gradient of Lagrangian at the (last) Newton step.
:``mu_``: barrier parameter at the (last) Newton step.
:``u_``: vector stacked with all primal variables at the (last) Newton step.
:``F_``: vector stacked with all equalities at the (last) Newton step.
:``G_``: vector stacked with all inequalities at the (last) Newton step.
:``nu_``: vector stacked with all dual equality variables at the (last) Newton step.
:``lambda_``: vector stacked with all dual inequality variables at the
              (last) Newton step.

.. comment .. warning:: INCOMPLETE - MISSING SENSITIVITY VARIABLES (perhaps not needed if using |Tvars2optimizeCS|)

.. warning::
   
   To be able to include these variables as input parameters, they
   have to be previously created using |tvariable| *with the
   appropriate sizes*.  Eventually, their values will be overridden by
   the solver to reflect the values listed above.
   
.. _`sensitivity-variables`:

Sensitivity variables
---------------------

.. warning:: This section of the documentation is still **incomplete**.

.. function:: Tvars2optimizeCS(parameter1,value1,parameter2,values2, ...)

.. list-table:: Selected parameters for |Tvars2optimizeCS|. For the
                full set of parameters use ``Tvars2optimizeCS help``
   :header-rows: 1
   :stub-columns: 1
   :widths: 1 1 3
   :align: left
      
   * - Parameter
     - Allowed values
     - Description

   * - 'optimizationVariables'
     - cell-array of |tcalculus| tensor variables created using |tvariable|
     - Variables to be optimized.
       
   * - 'sensitivityVariables'
     - cell-array of |tcalculus| tensor variables created using |tvariable|
     - Optimization variables with respect to which we want to compute cost
       sensitivity.
	     
   * - 'objective'
     - scalar |tcalculus| tensor
     - Criterion to optimize.
   
   * - 'constraints'
     - cell-array of |tcalculus| tensors, each involving one of the
       following operations ``==``, ``>=``, ``<=``, ``>``, ``<``
     - Optimization constraints.

   * - 'addEye2Hessian'
     - nonnegative real, default 1e-9
     - Add to the Hessian matrix appropriate identity matrices scaled by this constant.
       
       A larger value for addEye2Hessian has two main effects:
       
   * - 'smallerNewtonMatrix'
     - [``false``, ``true``], default ``false``
     - When ``true`` the matrix that needs to be inverted to compute a
       Newton step is reduced by first eliminating the dual variables
       associated with inequality constraints.  However, often the
       smaller matrix is not as sparse so the computation may actually
       increase.

.. _`solver-class`:

Solver Class
------------

The |matlab| class created to set parameter values and call the solver
from within |matlab| has the following methods:

.. function:: obj=classname()

   creates class and, for |C| code, loads the dynamic library containing the |C| code
		  
.. function:: delete(obj)

   deletes the class and, for |C| code, unloads the dynamic library
   
.. function:: setP_{parameter}(obj,value)

   sets the value of one of the parameters
   
.. function:: setV_{variable}(obj,value)

   sets the value of one of the optimization variables
   
.. function:: [y1,y2, ...]=getOutputs(obj)

   gets the values of the ``outputExpressions``
   
.. function:: [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter))

   :param mu0: initial value for the barrier variable
   :param maxIter: maximum number of Newton iterations
   :param saveIter: iteration # when to save the "hessian" matrix (for
                    subsequent pivoting/permutations/scaling
                    optimization) only saves when allowSave is true.

		    * When ``saveIter=0``, the hessian matrix is saved at
		      the last iteration; and
	       
		    * When ``saveIter<0``, the hessian matrix is not saved.

		    The "hessian" matrix will be saved regardless of
		    the value of ``saveIter``, when the solver exists with
		    ``status=4``.
	       
   :returns status: solver exist status
		  
                  *  0  = success
                  *  >0 = solver terminated unexpectedly
		     
                  * nonzero status indicates the reason for termination
                    in a binary format:
		    
		    - bit  0 = 1 - (primal) variables violate constraints
                    - bit  1 = 1 - dual variables are negative
                    - bit  2 = 1 - failed to invert hessian
                    - bit  3 = 1 - maximum # of iterations reached
		      
		    when the solver exists because the maximum # of iterations
                    was reached (bit 3 = 1), the remaining bits provide
                    information about the solution returned
		   
                    - bit  4 = 1 - gradient larger then ``gradTolerance``
                    - bit  5 = 1 - equality constraints violate ``equalTolerance``
                    - bit  6 = 1 - duality gap larger than ``desiredDualityGap``
                    - bit  7 = 1 - barrier variable larger than minimum value
                    - bit  8 = 1 - scalar gain alpha in Newton direction smaller than ``alphaMin``
                    - bit  9 = 1 - scalar gain alpha in Newton direction smaller than .1
                    - bit 10 = 1 - scalar gain alpha in Newton direction smaller than .5
   :returns iter: number of iterations
   :returns time: solver's compute time (in secs).

Pedigrees
---------

Pedigrees can save a lot of time for functions that take some time to
execute, like the |tenscalc|'s functions that generate code::

    |cmex2optimizeCS|      |cmex2equilibriumLatentCS|      |cmex2compute|
    |class2optimizeCS|     |class2equilibriumLatentCS|     |class2compute|

Essentially, every time a function is executed with pedigrees enables,
all its inputs are saved in a *pedigree* file and the function's
outputs are also saved. In case the function produces output files
(e.g., |C| or |matlab| code), the files are stored with unique names
for possible subsequent reuse. When the function is called again, it
checks whether it has been previously called with *the same exact
inputs*:

* If it has, then the previously saved outputs can be retrieved from
  the appropriate files and the function does not need to be
  recomputed.

* Otherwise, the function is executed and an additional pedigree and
  the associated outputs are saved for potential subsequent use.

Pedigrees are enables by specifying the input parameter
``'pedigreeClass'``, which is the common prefix used for the names of
all the files used to save the pedigree and outputs of the
function. Each time the function is called, this prefix is augmented
with a unique suffix that reflects the date and time the function was
called. To be precise, the *first* date/time the function was called
with that particular set of inputs. Any subsequent calls with the same
exact inputs reuse the previously saved pedigree.

The following call to |cmex2optimizeCS| enables the used of pedigrees::

  classname=cmex2optimizeCS('pedigreeClass','tmp_myopt', ...
                            'executeScript','asneeded',... 
                            'objective',norm2(x-y), ...
                            'optimizationVariables', { x, y }, ...
	  		    'constraints', { x>=-1, x<=1, y>=10 }, ...
			    'outputExpressions', { x, y });

and an instance of the solver generated by this call to
|cmex2optimizeCS| can be created using::
  
  obj=feval(classname);
	
In the call to |cmex2optimizeCS|, the parameter pair:
    ``'pedigreeClass','tmp_myopt'``

specifies that pedigrees should be enabled and that all pedigree files
should start with the prefix ``tmp_myopt``. The |cmex2optimizeCS|
parameter pair:

    ``'executeScript','asneeded'``

further specifies that the |cmex2optimizeCS| function should only be
executed if "it is needed". Specifically, if |cmex2optimizeCS| has
previously been called with the same exact set of inputs, then
|cmex2optimizeCS| should *not* be executed and, instead, the
previously generated code should be reused. Alternatively,

    ``'executeScript','yes'``

would specify that the code should be regenerated again *even it*
|cmex2optimizeCS| *has been called before with the same exact inputs.*

.. warning::

   When pedigrees are enabled through ``'pedigreeClass'`` one should
   **not** specify the class name using the parameter ``'classname'``.

   The actual classname will be chosen so that it is unique for a
   specific set of inputs and returned to the user as the output. In
   this way, regenerating code with a new set of inputs will not
   overwrite existing code.

.. note::

   Every time a function that uses pedigrees is called with a
   different set of inputs, a new pedigree is created for potential
   subsequent reuse. Because of this, one can easily get 100s of
   pedigree files. The good news is that pedigree files can be safely
   removed, because they are simply used to save time.

   It is a good practice to name all your pedigree files with a unique
   prefix that marks the file as "safe-to-remove". In all |tenscalc|
   examples, we use the prefix ``'tmp'`` to mark these files, which
   means that all files started with the 3 letters ``'tmp'`` are safe
   to remove.


.. note::
   
   Pedigrees are stored both as a ``.mat`` file (for fast retrieval)
   as well as a human-readable ``.html`` file. In general, there is
   little reason to dig into pedigree files, but all the inputs are
   there in case one is curious about previous calls to
   code-generation functions.

   Pedigrees are enabled by the |funpartools| toolbox and are
   available to any function that uses this toolbox to process input
   and output parameters.


Nash equilibrium
================

|tenscalc| can generate optimized code to compute Nash equilibrium

.. math::
   u^*\in\arg \min_u\Big\{ f(u,d^*,x) : F(u,d^*,x)\ge 0, G(u,d^*,x)=0 \Big\}

   d^*\in\arg \min_d\Big\{ g(u^*,d,x) : F(u^*,d,x)\ge 0, G(u^*,d,x)=0 \Big\}

where the variables ``u``, ``d``, ``x`` can include multiple tensors
and the equality and inequality constraints can be expressed by
equalities and inequalities involving multiple tensors.

The *latent* variable ``x`` that appears in both minimizations, must
be fully determined by the equality constraints and is thus not really
a free optimization variable.

The following two scripts are used to generate code to compute this
type of Nash equilibrium

.. function:: cmex2equilibriumLatentCS(parameter1,value1,parameter2,values2, ...)
.. function:: class2equilibriumLatentCS(parameter1,value1,parameter2,values2, ...)

   :param parameter1: parameter to set
   :type parameter1: string
   :param value1: parameter to set
   :type value1: type depends on the parameter
   :param parameter2: parameter to set
   :type parameter2: string
   :param value2: parameter to set, ...
   :type value2: type depends on the parameter
   :returns: name of the |matlab| class created
   :rtype: string
		       
The function |cmex2optimizeCS| generates |C| code, whereas
|class2optimizeCS| generates |matlab| code, but both functions take
the same set of parameters and generate |matlab| classes that are
indistinguishable to the user.

.. list-table:: Selected parameters for |c| and
                |class2equilibriumLatentCS|. For the full set of
                parameters use ``class2equilibriumLatentCS help`` or
                ``class2equilibriumLatentCS help``
   :header-rows: 1
   :stub-columns: 1
   :widths: 2 3 3
   :align: left
      
   * - Parameter
     - Allowed values
     - Description

   * - 'P1optimizationVariables'
     
       'P2optimizationVariables'
     - cell-array of |tcalculus| tensor variables created using |tvariable|
     - Variables to be optimized by player 1 and player 2
       
   * - 'P1objective'
     
       'P2objective'
     - scalar |tcalculus| tensor
     - Criteria to optimize for player 1 and player 2
   
   * - 'P1constraints'
     
       'P2constraints'
     - cell-array of |tcalculus| tensors, each involving one of the
       following operations ``==``, ``>=``, ``<=``, ``>``, ``<``
     - Optimization constraints for player 1 and player 2

       .. warning::

	  For constraint satisfaction, there is no difference between
	  ``<=`` and ``<`` or between ``>=`` and ``>``.

   * - 'latentVariables'
     - cell-array of |tcalculus| tensor variables created using |tvariable|
     - Latent optimization variables common to both players

   * - 'latentConstraints'
     - cell-array of |tcalculus| tensors, each involving one of the
       following operations ``==``, ``>=``, ``<=``, ``>``, ``<``
     - Optimizations constraints common to both players that
       implicitely define the latext variables.

   * - 'outputExpressions'
     - cell-array of |tcalculus| tensors
     - Expressions (typically involving the optimization variables)
       that the solver should return upon termination.

       See :ref:`special-variables-nash`:

.. _`special-variables-nash`:

Special variables to include in ``'outputExpressions'`` for Nash solver
-----------------------------------------------------------------------

The following |tcalculus| variables are assigned special values and
can be using in ``outputExpressions:``
       
:``P1lambda1_``, ``P1lambda2_``, ... , ``P2lambda1_``, ``P2lambda2_``, ...:
   Lagrangian multipliers associated with the inequalities constraints
   for player 1 and 2 (in the order that they appear and with
   the same size as the corresponding constraints)
:``P1nu1_``, ``P1nu2_``, ... , ``P2nu1_``, ``P2nu2_``, ... ``P1xnu1_``, ``P1xnu2_``, ... , ``P2xnu1_``, ``P2xnu2_``, ...:
   Lagrangian multipliers associated with the equality constraints
   player 1 and 2 (in the order that they appear and with the same
   size as the corresponding constraints). The ``P1x`` and ``P2x``
   variables correspond to the ``latentConstraints``.
:``Hess_``: Hessian matrix used by the (last) Newton step to update
            the primal variables (not including addEye2Hessian).

.. warning::
   
   To be able to include these variables as input parameters, they
   have to be previously created using |tvariable| *with the
   appropriate sizes*.  Eventually, their values will be overridden by
   the solver to reflect the values listed above.


.. include:: defs.rst

==================================
 Code generation for computations
==================================

.. warning:: This section of the documentation is still **incomplete**.
   

**Computation trees**

  Sets of computations organized as dependency graphs, with multiple
  roots (input variables) and leaves (output expressions).

Defining a computation
======================

.. class:: csparse

.. method:: cs=csparse()

   Creates an empty computation tree
	    
   :returns: empty |csparse| object
		      
.. method:: declareSet(cs,destination,functionName)

   Adds an input variable to the tree
	    
   :param cs: |csparse| object
   :param destination: |tcalculus| tensor variable created using |tvariable|
   :param functionname: name of the function to be created
   :type functionname: string
	    
.. method:: declareGet(cs,sources,functionName)
	    
   Adds output expressions to the tree
	    
   :param cs: |csparse| object
   :param sources: cell array of |tcalculus| tensor-valued expressions
   :param functionname: name of the function to be created
   :type functionname: string

.. method:: declareCopy(cs,destinations,sources,functionName)

   Sets the value of input variables, based on the value of output
   expressions

   :param cs: |csparse| object
   :param sources: cell array of |tcalculus| tensor-valued expressions
   :param destinations: cell array of |tcalculus| tensor variables created using |tvariable|
   :param functionname: name of the function to be created
   :type functionname: string

		       
Code generation
===============

.. function:: cmex2compute(parameter1,value1,parameter2,values2, ...)
.. function:: class2compute(parameter1,value1,parameter2,values2, ...)

   Generate code to set values of input variables, get values of
   output expressions, and copy output expressions to input variables.
	      
   :param parameter1: parameter to set
   :type parameter1: string
   :param parameter2: parameter to set, ...
   :type parameter2: string
   :returns: name of the |matlab| class created
   :rtype: string
		       
   




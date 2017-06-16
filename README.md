# TensCalc

## Table of Contents

* [Description](#description)
* [Installation](#installation)
* [Usage](#usage)
* [Issues](#issues)
* [Contact Information](#contact-information)
* [License Information](#license-information)

## Description

The *TensCalc* *Matlab* toolbox provides an environments for
performing nonlinear constrained optimization. 

The variables to be optimized can be multi-dimensional arrays of any
dimension (tensors) and the cost functions and inequality constraints
are specified using *Matlab*-like formulas. 

Interior point methods are used for the numerical optimization, which
uses formulas for the gradient and the hessian matrix that are
computed symbolically in an automated fashion. 

The package can either produce optimized *Matlab* code or C code. The
former is preferable for very large problems, whereas the latter for
small to mid- size problems that need to be solved in just a few
milliseconds. The C code can be used from inside *Matlab* using an
(automatically generated) cmex interface or in standalone
applications. No libraries are required for the standalone code.

## Installation

The *TensCalc* toolbox supports *Matlab* running under:

- OSX (tested extensivly)
- linux (tested lightly)
- Microsoft Windows (very little testing)

To install

1. Install the [FunParTools](https://github.ucsb.edu/hespanha/funpartools) toolbox.

2. Install the [CmexTools](https://github.ucsb.edu/hespanha/cmextools) toolbox.

3. Download *TensCalc* using one of the following options:

	1. downloading it as a zip file from
		https://github.ucsb.edu/hespanha/tenscalc/archive/master.zip
	   and unziping to an appropriate location

	2. cloning this repository with Git, e.g., using the shell command
	   ```sh
	   svn checkout https://github.ucsb.edu/hespanha/tenscalc.git
	   ```
	  
	3. checking out this repository with svn, e.g., using the shell command
	   ```sh
       git clone https://github.ucsb.edu/hespanha/tenscalc.git
       ```

	The latter two options are recommended because you can
    subsequently use `svn update` or `git pull` to upgrade *TensCalc*
    to the latest version.

	After this, you should have at least the following folders:

	* `tenscalc/lib`
	* `tenscalc/examples`
	* `tenscalc/doc`

4. Add `tenscalc/lib` and `tenscalc/lib/csparse` to your matlab path. 
   From inside the folder `tenscalc/lib`, this can be done with

	```matlab
	addpath(fileparts(which(`class2optimizeCS`)));
	addpath([fileparts(which(`class2optimizeCS`)),'/csparse']);
	savepath
	```

5. Compile a set of `cmex` files needed by *TensCals* using the
   *Matlab*, from inside the folder `tenscalc/lib`:

	```matlab
	compileInstructionsTable
	```
	
	This requires the *CmexTools* toolbox.
	
6. To test if all is well, go to `tenscalc/examples` and try a few example, such as

	```matlab
	mls
	sls
	l1l2estimationCS
	```
	
## Usage

### What are tensors?

Tensors are essentially multi-dimensional arrays, but one needs to
keep in mind that in *Matlab* every variable is an array of dimension 2
or larger. Unfortunately, this is not suitable for *TensCalc*, which
also needs arrays of dimension 0 (i.e., scalars) and 1 (i.e.,
vectors). This can create con- fusion because *Matlab* automatically
“upgrades” scalars and vectors to matrices (by adding singleton
dimensions), but this is not done for *TensCalc* expressions. More on
this later, but for now let us keep going with the quick start.

### STVEs?

The basic objects in *TensCalc* are symbolic tensor-valued expressions
(STVEs). These ex- pressions typically involve symbolic variables that
can be manipulated symbolically, evaluated for specific values of its
variables, and optimized.

### Need speed?

Prior to numerical optimization, STVEs must be “compiled” for
efficient computation. This compilation can take a few seconds or even
minutes but results in highly efficient *Matlab* or C code. Big payoffs
arise when you need to evaluate or optimize an expression multiple
time, for different values of input variables. *TensCalc*’s compilation
functions thus always ask you to specify input parameters. Much more
on *TensCalc*’s compilations tools can be found in CSparse’s
documentation.

### Creating STVEs. 

The following sequence of *TensCalc* command can be used to declare an STVE to be
used in a simple least-squares optimization problem

``` matlab
N=100; n=8;
Tvariable A [N,n];
Tvariable b N;
Tvariable x n;
y=A*x-b; 
J=norm2(y)
```

### Optimizing STVEs

To perform an optimization also need to create an appropriate
specialized *Matlab* class, say called minslsu, using the following
command:
```matlab
cmex2optimizeCS('classname','Cminslsc',...
                'method','primalDual',...
                'objective',J,...
                'optimizationVariables',{x},...
                'constraints',{x>=0,x<=.05},...
                'outputExpressions',{J,x},...
                'parameters',{A,b},...
                'solverVerboseLevel',2);
```

The goal of this class is to minimize the symbolic expression `J` with
respect to the variable `x`, subject to the constrains 'x>=0' and
'x<=.05`. The symbolic variables `A` and `b` are declared as
parameters that can be changed from optimization to
optimization. Setting the `solverVerboseLevel` to 3, asks for a
moderate amount of debugging information to be printed while the
solver is executed (one line per iteration of the solver).

Once can see the methods available for the class `Cminslsc`
generated by `cmex2optimizeCS` using the usual
`help` command, which produces:

```matlab
>> help minslsu
% Create object
obj=minslsu();
% Set parameters
setP_A(obj,{[100,8] matrix});
setP_b(obj,{[100,1] matrix});
% Initialize primal variables
setV_x(obj,{[8,1] matrix});
% Solve optimization
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[y1,y2]=getOutputs(obj);
```

The following commands creates an instance of the class and preforms
the optimization for specific parameter values:

```matlab
thisA=rand(N,n);
thisb=rand(N,1);
x0=.02*rand(n,1);

obj=Cminslsc();
setP_A(obj,thisA);
setP_b(obj,thisb);
setV_x(obj,x0);
mu0=1;
maxIter=20;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(-1));
[Jcstar,xcstar]=getOutputs(obj);
```

The parameters `mu0` and `maxIter` passed to the solver are the
initial value of the barrier variable and the maximum number of Newton
iterations, respectively.

### Need more examples?

The example above and many others can be found in `tenscalc\examples`.

### Full documentaion

Full documentation for this toolbox can be found in

* `doc/tenscalc.pdf`

additiona tecnhnical information can be found at

* `doc/ipm.pdf`
* `doc/csparse.pdf`
* `doc/computationgraphs.pdf`
* `doc/timeseries.pdf`


## Issues

* While most *Matlab* scripts are agnostic to the underlying operating
  systems (OSs), the use of `mex` functions depends heavily on the
  operating systems.

  Our goal is to build a toolbox that works across multiple OSs; at
  least under OSX, linux, and Microsft Windows. However, most of our
  testing was done under OSX so one should expect some bugs under the
  other OSs. Sorry about that.

* Getting the solver to converge can be difficult for problems that
  are numerically ill conditioned.

* The next biggest issue is the use of fairly obscure error messages
  when things go wrong.
  

## Contact Information

Joao Hespanha (hespanha@ucsb.edu)

http://www.ece.ucsb.edu/~hespanha

University of California, Santa Barbara
	
## License Information

Copyright 2010-2017 Joao Hespanha

This file is part of Tencalc.

TensCalc is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

TensCalc is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.

function types=instructionTypes()
%
% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.
    
types={
    'I_set'                % no operation needed (value determined by set)
                           %   no parameters
    'I_load'               % load memory location with constant value
                           %   operands = none
                           %   parameters = scalar double with value
                           %                to be loaded
%% C/ASM instructions
    'I_sum'                % sum(parameters.*operands)
                           %   parameters = array with +1/-1 specifying
                           %                addition/subtraction for each
                           %                element of vector
    'I_sumprod'            % sum(prod(operands,1))
                           %      a1*a2*...*am+b1*b2*...*bm+...
                           %   operands = a1,a2,...,am,b1,b2,...,bm,c1,...
                           % In case operands are one long vector:
                           %   parameters = [m, # sums]
    'I_inv'                % inverse: 1/denominator
                           %   operands = denominator
    'I_minus_inv_sqr'      % inverse: -1/denominator^2
                           %   operands = denominator
    'I_2_inv_cube'         % inverse: 2/denominator^3
                           %   operands = denominator
    'I_div'                % division: numerator/denominator
                           %   operands = numerator, denominator
    'I_minus_dot'          % symmetric of dot-product: -a1*b1-a2*b2-...
                           %   operands = a1,b1,a2,b2,....
    'I_minus_dot_div'      % symmetric of dot-product & divide: (-a1*b1-a2*b2-...)/a0
                           %   operands = a0,a1,b1,a2,b2,....
    'I_plus_minus_dot'     % subtract dot-product: a0-a1*b1-a2*b2-...
                           %   operands = a0,a1,b1,a2,b2,....
    'I_plus_sqr'           % norm-2 square: a1^2+a2^2+...
                           %   operands = a1,a2,....
    'I_plus_minus_dot_div' % subtract dot-product & divide: (a0-a1*b1-a2*b2-...)/b0
                           %   operands = a0,b0,a1,b1,a2,b2,....
    'I_min'                % computes the min: min(a1,a2,...);
                           %   operands = a1,a2,....
    'I_min0'               % computes the min: min(0,a1,a2,...);
                           %   operands = a1,a2,....
    'I_max'                % computes the max: max(a1,a2,...);
                           %   operands = a1,a2,....
    'I_max0'               % computes the max: max(0,a1,a2,...);
                           %   operands = a1,a2,....
    'I_max_abs'            % infinity norm: max{ |a1|,|a2|,... }
                           %   operands = a1,a2,....
    'I_clp'                % computes the largest c>=0 such that: ai + c bi >=0
                           %   operands = a1,b1,a2,b2,....an,bn
    'I_exp'                % computes exp(a1)
                           %   operand = a1
    'I_log'                % computes log(a1)
                           %   operand = a1
    'I_cos'                % computes cos(a1)
                           %   operand = a1
    'I_minus_cos'          % computes -cos(a1)
                           %   operand = a1
    'I_sin'                % computes sin(a1)
                           %   operand = a1
    'I_minus_sin'          % computes -sin(a1)
                           %   operand = a1
    'I_abs'                % computes abs(a1)
                           %   operand = a1
    'I_sqrt'               % computes sqrt(a1)
                           %   operand = a1
    'I_Dsqrt'              % computes .5/sqrt(a1)
                           %   operand = a1
    'I_DDsqrt'             % computes -.25/sqrt(a1)^3
                           %   operand = a1
    'I_sqr'                % computes a1^2
                           %   operand = a1
    'I_2times'             % computes 2*a1
                           %   operand = a1
    'I_2'                  % computes 2
                           %   operand = a1 (ignored)
    'I_cube'               % computes a1^3
                           %   operand = a1
    'I_3sqr'               % computes 3*a1^2
                           %   operand = a1
    'I_6times'             % computes 6*a1
                           %   operand = a1
    'I_atan'               % computes atan(a1)
                           %   operand = a1
    'I_Datan'              % computes 1/(1+a1^2)
                           %   operand = a1
    'I_DDatan'             % computes -2*a1/(1+a1^2)^2
                           %   operand = a1;
    
%% Atomic C instructions;
    'I_luS2A'              % performs an atomic LU factorization of a
                           % matrix with entries in the scratchbook
                           %   parameters = [ LU atomic variable ID,
                           %                  #rows (=#cols)]
                           %   operands = nnz elements (in order of ir)
    'I_luS2Asym'           % performs an atomic LU factorization of a
                           % matrix with entries in the scratchbook
                           %   parameters = [ LU atomic variable ID,
                           %                  #rows (=#cols)]
                           %   operands = nnz elements (in order of ir)
    'I_mldivideA2F1'       % given an atomic LU factorization for a
                           % matrix A and a vector b, solves the
                           % system of equations A x = b;
                           % and returns the 1st entry of x
                           %   parameters = []
                           %   operands [ I_luS2A instruction, b entries ]
    'I_mldivideA2Fn'       % copy one entry of solution x from
                           % I_mldivideA2F1 to scratchbook
                           %   parameters = [ entry # (1-based)]
%                          %   operands = [I_luS2A instruction,I_mldivideA2F1 intruction ]
    
%% Matlab instructions
    'I_Mnorm2'
    'I_Mnorminf'
    'I_Mplus'
    'I_Mmtimes'
    'I_Mtimes'
    'I_Msum'
    'I_Mmin'
    'I_Mones'
    'I_Mzeros'
    'I_Meye'
    'I_Mdiag'
    'I_Mclp'
    'I_Mctranspose'
    'I_Msubsref'
    'I_Mtprod'
    'I_Mtprod_matlab'
    'I_Mfull'         
    'I_Mreshape'         
    'I_Mrepmat'         
    'I_Mcat'         
    'I_Mlu'
    'I_Mldl'
    'I_Mchol'
    'I_Mmldivide_l1'
    'I_Mmldivide_u'
    'I_Mmldivide_u1'
    'I_Mmldivide_d'
    'I_Mrdivide'
    'I_Mcompose'
};

m=mat2cell(int32(1:length(types)),1,ones(1,length(types)));
types=cell2struct(m,types,2);

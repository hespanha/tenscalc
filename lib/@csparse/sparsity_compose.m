function [subsY,instrY,instructions]=sparsity_compose(obj,thisExp,instruction)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'compose' that preserves sparsity.
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

osize=getOne(obj.vectorizedOperations,'osize',thisExp);
operands=getOne(obj.vectorizedOperations,'operands',thisExp);

subsY=getOne(obj.vectorizedOperations,'subscripts',operands(1));
instr1=getOne(obj.vectorizedOperations,'instructions',operands(1));

instrY=newInstructions(obj,instruction,{[]},num2cell(instr1,2),thisExp);

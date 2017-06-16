function [subsY,instrY,instructions]=sparsity_max(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'max_nz', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
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
dimension=getOne(obj.vectorizedOperations,'parameters',thisExp);
osize1=getOne(obj.vectorizedOperations,'osize',operands(1));

subsX=getOne(obj.vectorizedOperations,'subscripts',operands(1));
instrX=getOne(obj.vectorizedOperations,'instructions',operands(1));

subsY=subsX;
subsY(dimension,:)=[];
[subsY,iSmall,iLarge]=unique(subsY','rows');
subsY=subsY';

instrY=nan(length(iSmall),1);
for i=1:length(iSmall)
    k=(iLarge==iSmall(i));
    operands=instrX(k);
    if length(operands)==prod(osize1(dimension))
        if length(operands)>1
            instrY(i)=newInstruction(obj,obj.Itypes.I_max,[],operands,thisExp);
        else
            instrY(i)=operands;
        end
    else
        instrY(i)=newInstruction(obj,obj.Itypes.I_max0,[],operands,thisExp);
    end
end



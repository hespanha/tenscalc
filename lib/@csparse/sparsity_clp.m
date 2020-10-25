function [subsY,instrY,instructions]=sparsity_clp(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'clp', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements. For x >= 0,
%    clp(x,dx) = max { alpha >0 : x + alpha dx >= 0 }
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

subsY=zeros(0,1);

%% Get instructions for operands
subsX=getOne(obj.vectorizedOperations,'subscripts',operands(1));
instrX=getOne(obj.vectorizedOperations,'instructions',operands(1));

subsdX=getOne(obj.vectorizedOperations,'subscripts',operands(2));
instrdX=getOne(obj.vectorizedOperations,'instructions',operands(2));

[kX,kdX]=ismember(subsX',subsdX','rows');

if any(~ismember(subsdX',subsX','rows'))
    error(['sparsity_clp: structurally zero values for x not allowed ' ...
           'when dx is not structurally zero'])
end

%% Determine instructions for Y
operands=[instrX(kX)';instrdX(kdX(kX))'];
instrY=newInstruction(obj,obj.Itypes.I_clp,[],operands(:),thisExp);

%disp(obj)

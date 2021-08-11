function [subsY,instrY,instructions]=sparsity_clp(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'clp', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements. For x >= 0,
%    clp(x,dx) = max { alpha >0 : x + alpha dx >= 0 }
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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
end

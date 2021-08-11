function [subsY,instrY,instructions]=sparsity_norm1(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'norm2', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    osize=getOne(obj.vectorizedOperations,'osize',thisExp);
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);

    subsY=zeros(0,1);

    %% Get instructions for operand
    subsX=getOne(obj.vectorizedOperations,'subscripts',operands(1))';
    instrX=getOne(obj.vectorizedOperations,'instructions',operands(1));


    %% Determine instructions for Y
    instrY=newInstruction(obj,obj.Itypes.I_plus_abs,[],instrX(:)',thisExp);
end

function [subsY,instrY,instructions]=sparsity_compose(obj,thisExp,instruction)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'compose' that preserves sparsity.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    osize=getOne(obj.vectorizedOperations,'osize',thisExp);
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);

    subsY=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instr1=getOne(obj.vectorizedOperations,'instructions',operands(1));

    instrY=newInstructions(obj,instruction,{[]},num2cell(instr1,2),thisExp);
end
function [subsY,instrY,instructions]=sparsity_compose_full(obj,thisExp,instruction)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'compose' that completely fills the matrix.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    osize=getOne(obj.vectorizedOperations,'osize',thisExp);
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);

    subs1=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instr1=getOne(obj.vectorizedOperations,'instructions',operands(1));

    subsY=memory2subscript(osize,1:prod(osize));
    [lia,locb]=ismember(subsY',subs1','rows');

    % fill matrix
    instrFull=nan(size(subsY,2),1);
    % non-zero
    instrFull(lia)=instr1(locb(lia));
    if ~all(lia)
        % zero
        instrFull(~lia)=newInstructions(obj,obj.Itypes.I_load,num2cell(zeros(sum(~lia),1)),{[]},thisExp);
    end
    instrY=newInstructions(obj,instruction,{[]},num2cell(instrFull,2),thisExp);
end
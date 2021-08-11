function [subsY,instrY]=sparsity_logdet_ldl(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'logdet_ldl', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% logdet_ldl(LDL):
% 1) computes det(D) by
%    . computing the product of the main diagonal entries of LDL
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    verboseLevel=0;

    operands=getOne(obj.vectorizedOperations,'operands',thisExp);

    osizeLDL=getOne(obj.vectorizedOperations,'osize',operands(1));
    subsLDL=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrLDL=getOne(obj.vectorizedOperations,'instructions',operands(1));
    if length(osizeLDL)~=2 || osizeLDL(1)~=osizeLDL(2)
        osizeLDL
        error('sparsity_logdet_ldl: in logdet_ldl(LDL), LDL must be a square (2D) matrix\n');
    else
        n=double(osizeLDL(1));
    end

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsity_logdet_ldl (%3d): LDL      size=%-10s, nnz=%d\n',...
                thisExp,['[',index2str(osizeLDL),']'],length(instrLDL));
    end

    % find instructions for diagonal
    k=find(subsLDL(1,:)==subsLDL(2,:)); % find diagonal elements
    instrX=instrLDL(k);

    subsY=zeros(0,1);
    if length(instrX)==n
        % log(prod)
        %instrY=newInstruction(obj,obj.Itypes.I_sumprod,[length(instrX),1],instrX(:)',thisExp);
        %instrY=newInstruction(obj,obj.Itypes.I_log,[],instrY,thisExp);
        % sum(log)
        instrY=newInstructions(obj,obj.Itypes.I_abs,{[]},num2cell(instrX,2),thisExp);
        instrY=newInstructions(obj,obj.Itypes.I_log,{[]},num2cell(instrY,2),thisExp);
        parameters=ones(1,numel(instrY));
        instrY=newInstruction(obj,obj.Itypes.I_sum,parameters,instrY,thisExp);
    else
        error('sparsity_logdet_ldl: LDL structurally singular in logdet_ldl(LDL)\n');
    end
end

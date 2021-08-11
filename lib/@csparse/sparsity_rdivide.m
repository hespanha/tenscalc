function [subsY,instrY,instructions]=sparsity_rdivide(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'rdivide', allocates the required memory, and
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

    %% Get instructions for operands
    osize1=getOne(obj.vectorizedOperations,'osize',operands(1));
    subsX1=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrX1=getOne(obj.vectorizedOperations,'instructions',operands(1));

    osize2=getOne(obj.vectorizedOperations,'osize',operands(2));
    subsX2=getOne(obj.vectorizedOperations,'subscripts',operands(2));
    instrX2=getOne(obj.vectorizedOperations,'instructions',operands(2));

    if isempty(osize2)
        % division by scalar
        %fprintf('sparsity_rdivide: division by scalar\n');
        if isempty(instrX2)
            error('sparsity_rdivide: structurally zero scalar for x2 not allowed in x1./x2')
        end
        subsY=subsX1;
        %% Determine instructions for Y
        operands=[instrX1';instrX2*ones(1,length(instrX1))];
    elseif isempty(osize1)
        % scalar divided by
        %fprintf('sparsity_rdivide: scalar divided by\n');
        if isempty(instrX2)
            error('sparsity_rdivide: structurally zero scalar for x2 not allowed in x1./x2')
        end
        subsY=subsX2;
        %% Determine instructions for Y
        operands=[instrX1*ones(1,length(instrX2));instrX2'];
    elseif myisequal(osize1,osize2)
        % entry-wise division
        %fprintf('sparsity_rdivide: entry-wise division\n');
        if size(subsX2,2)~=prod(osize)
            error('sparsity_rdivide: structurally zero values for x2 not allowed in x1./x2')
        end
        [kX2,kX1]=ismember(subsX2',subsX1','rows');
        subsY=subsX1(:,kX1(kX2));
        %% Determine instructions for Y
        operands=[instrX1(kX1(kX2))';instrX2(kX2)'];
    else
        osize1,osize2,;
        error('sparsity_rdivide: mismatched sizes in x1./x2');
    end
    if isempty(operands)
        instrY=zeros(0,1);
    else
        operands=mat2cell(operands,2,ones(size(operands,2),1));
        instrY=newInstructions(obj,obj.Itypes.I_div,{[]},operands,thisExp);
    end
end

function [subsX,instrX]=sparsity_lu_d(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'lu_d', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% lu_d(LDL):
% 1) computes D by
%    . extracting the main diagonal entries of LU
%    . adjust sign to account for permutations
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    verboseLevel=0;

    operands=getOne(obj.vectorizedOperations,'operands',thisExp);

    osizeLU=getOne(obj.vectorizedOperations,'osize',operands(1));
    subsLU=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrLU=getOne(obj.vectorizedOperations,'instructions',operands(1));
    pq=getOne(obj.vectorizedOperations,'parameters',operands(1));

    if length(osizeLU)~=2 || osizeLU(1)~=osizeLU(2)
        osizeLU
        error('lu_d: in lu_l(LU), LU must be a square (2D) matrix\n');
    else
        n=double(osizeLU(1));
    end

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsify_lu_d (%3d): LU      size=%-10s, nnz=%d\n',...
                thisExp,['[',index2str(osizeLU),']'],length(instrLU));
    end

    % Compute instructions
    k=find(subsLU(1,:)==subsLU(2,:)); % find diagonal elements
    subsX=subsLU(1,k);
    instrX=instrLU(k);

    % keep subsX in the natural order
    [subsX,k]=sortrows(subsX');
    subsX=subsX';
    instrX=instrX(k);

    % undo sign changes in determinant by permutations
    si=det(sparse(1:n,pq{1},ones(n,1),n,n))*det(sparse(1:n,pq{2},ones(n,1),n,n));
    if si<0
        instructions=num2cell(instrX,2);
        instrX=newInstructions(obj,obj.Itypes.I_sum,{-1},instructions,thisExp);
    end

    if verboseLevel>0
        fprintf('  sparsify_lu_d (%3d): D     size=%-10s, nnz=%d,                          # new instr=%4d (%d..%d) (%.2f sec)\n',...
                thisExp,['[',index2str(osizeLU(1)),']'],length(instrX),...
                instructionsTableHeight()-n0,n0+1,instructionsTableHeight(),etime(clock,t0));
    end
end

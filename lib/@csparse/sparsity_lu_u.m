function [subsX,instrX]=sparsity_lu_u(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'lu_d', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% lu_l(LDL):
% 1) computes U by
%    . extracting the lower diagonal entries of LU (including main diagonal)
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
    q=pq{2}; % column permutations
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
    k=find(subsLU(1,:)<=subsLU(2,:)); % find below diagonal
    subsX=[subsLU(1,k);subsLU(2,k)];
    instrX=instrLU(k);

    % apply q permutation to columns
    subsX(2,:)=q(subsX(2,:));

    % keep subsX in the natural order
    [subsX,k]=sortrows(subsX');
    subsX=subsX';
    instrX=instrX(k);

    if verboseLevel>0
        fprintf('  sparsify_lu_u (%3d): D     size=%-10s, nnz=%d,                          # new instr=%4d (%d..%d) (%.2f sec)\n',...
                thisExp,['[',index2str(osizeLU(1)),']'],length(instrX),...
                instructionsTableHeight()-n0,n0+1,instructionsTableHeight(),etime(clock,t0));
    end

end

function [subsX,instrX]=sparsity_ldl_l(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'ldl_l', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% ldl_l(LDL):
% 1) computes L by
%    . extracting the lower diagonal entries of LDL
%    . adding ones in the main diagonal
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
    q=getOne(obj.vectorizedOperations,'parameters',operands(1));
    q=q{1}; % column (and row) permutation
    if length(osizeLDL)~=2 || osizeLDL(1)~=osizeLDL(2)
        osizeLDL
        error('ldl_d: in ldl_l(LDL), LDL must be a square (2D) matrix\n');
    else
        n=double(osizeLDL(1));
    end

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsify_ldl_d (%3d): LDL      size=%-10s, nnz=%d\n',...
                thisExp,['[',index2str(osizeLDL),']'],length(instrLDL));
    end

    % Compute instructions
    k=find(subsLDL(1,:)>subsLDL(2,:)); % find below diagonal
    subsX=[subsLDL(1,k);subsLDL(2,k)];
    instrX=instrLDL(k);

    % add ones to diagonal
    subsX=[subsX,[1:osizeLDL(1);1:osizeLDL(1)]];
    instrX=[instrX;repmat(newInstructions(obj,obj.Itypes.I_load,...
                                          {1},{[]},thisExp),osizeLDL(1),1)];

    % apply q permutation to rows
    subsX(1,:)=q(subsX(1,:));

    % keep subsX in the natural order
    [subsX,k]=sortrows(subsX');
    subsX=subsX';
    instrX=instrX(k);

    if verboseLevel>0
        fprintf('  sparsify_ldl_d (%3d): D     size=%-10s, nnz=%d,                          # new instr=%4d (%d..%d) (%.2f sec)\n',...
                thisExp,['[',index2str(osizeLDL(1)),']'],length(instrX),...
                instructionsTableHeight()-n0,n0+1,instructionsTableHeight(),etime(clock,t0));
    end
end

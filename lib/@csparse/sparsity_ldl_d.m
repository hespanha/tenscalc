function [subsX,instrX]=sparsity_ldl_d(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'ldl_d', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% ldl_d(LDL):
% 1) computes D by
%    . extracting the main diagonal entries of LDL
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
    k=find(subsLDL(1,:)==subsLDL(2,:)); % find diagonal elements
    subsX=subsLDL(1,k);
    instrX=instrLDL(k);

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

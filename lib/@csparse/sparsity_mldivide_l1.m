function [subsX,instrX]=sparsity_mldivide_l1(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'mldivide_l1', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% mldivide_l1(LU,b):
% 1) computes L by
%    . extracting the stictly lower-triangular entries of LU
%    . adding 1 to the diagonal
% 2) applies any required row-permutation to b
%    (from the LU factorization)
% 3) solves L x = b
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    verboseLevel=1;

    operands=getOne(obj.vectorizedOperations,'operands',thisExp);

    osizeLU=getOne(obj.vectorizedOperations,'osize',operands(1));
    subsLU=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrLU=getOne(obj.vectorizedOperations,'instructions',operands(1));
    pq=getOne(obj.vectorizedOperations,'parameters',operands(1));
    p=pq{1}; % row permutation
    if length(osizeLU)~=2 || osizeLU(1)~=osizeLU(2)
        osizeLU
        error('mldivide_l1: in LU\\B, LU must be a square (2D) matrix\n');
    else
        n=double(osizeLU(1));
    end

    osizeB=getOne(obj.vectorizedOperations,'osize',operands(2));
    subsB=getOne(obj.vectorizedOperations,'subscripts',operands(2));
    instrB=getOne(obj.vectorizedOperations,'instructions',operands(2));

    if length(osizeB)==1
        subsB=[subsB;ones(1,size(subsB,2))];
        m=1;
    else
        m=double(osizeB(2));
    end
    if size(subsB,1)~=2 || osizeB(1)~=n
        size(subsB),osizeB
        error('mldivide: in (I+L)\\B, B must be a vector or a (2D) matrix with the same number of rows as L\n');
    end

    if verboseLevel<=1 && m*n>1e6
        verboseLevel=2;
        fprintf('    computing instructions for (large) mldivide:\n');
    end

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsify_mldivide_l1(%3d): LU      size=%-10s, nnz=%d, B  size=%-7s, nnz=%d\n',...
                thisExp,['[',index2str(osizeLU),']'],length(instrLU),...
                ['[',index2str(osizeB),']'],length(instrB));
    end

    % Reconstruct instructions of L
    k=find(subsLU(1,:)>subsLU(2,:));
    instrL=sparse(double(subsLU(1,k)),double(subsLU(2,k)),double(instrLU(k)),n,n);

    instrB=sparse(double(subsB(1,:)),double(subsB(2,:)),double(instrB),n,m);

    % apply p permutation to B
    instrB=instrB(p,:);

    %% Compute instructions
    instrX=instrB;
    instrXTranspose=instrX';
    cols=find(any(instrXTranspose,2))';
    for col=cols
        for row=2:n
            cond1=instrXTranspose(col,1:row-1);
            if any(cond1)
                cond2=instrL(row,1:row-1);
                if any(cond2)
                    k=find(cond1 & cond2);
                    if ~isempty(k)
                        if instrXTranspose(col,row)
                            operands=[instrL(row,k);instrXTranspose(col,k)];
                            operands=full([instrXTranspose(col,row),operands(:)']);
                            instrXTranspose(col,row)=newInstructions(obj,obj.Itypes.I_plus_minus_dot,...
                                                                     {[]},num2cell(operands,2),thisExp);
                        elseif instrXTranspose(col,row)==0
                            operands=[instrL(row,k);instrXTranspose(col,k)];
                            operands=full([operands(:)']);
                            instrXTranspose(col,row)=newInstructions(obj,obj.Itypes.I_minus_dot,...
                                                                     {[]},num2cell(operands,2),thisExp);
                        end
                    end
                end
            end
        end
    end
    instrX=instrXTranspose';
    % keep subsX in the natural order
    [i,j,instrX]=find(instrX);
    if length(osizeB)==1
        [subsX,k]=sort(uint64(i)');
    else
        [subsX,k]=sortrows([uint64(i),uint64(j)],2:-1:1);
        subsX=subsX';
    end
    instrX=instrX(k);

    if verboseLevel>0
        fprintf('  sparsify_mldivide_l1(%3d): (I+L)\\B size=%-10s, nnz=%d, # new instr=%4d (%d..%d) (%.2f sec)\n',...
                thisExp,['[',index2str(osizeB),']'],length(instrX),...
                instructionsTableHeight()-n0,n0+1,instructionsTableHeight(),etime(clock,t0));
    end
end

function test()

    n=3;
    m=2
    LU=rand(n,n);
    B=rand(n,m);

    [i,j,v]=find(LU);
    k=find(i>j);
    L=sparse(i(k),j(k),v(k),n,n);
    L1=L+speye(n);

    X=B;
    for col=1:m
        for row=2:n
            X(row,col)=X(row,col)-L(row,1:row-1)*X(1:row-1,col)
        end
    end
    X
    L1\B

end
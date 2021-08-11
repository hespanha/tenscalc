function [subsX,instrX]=sparsity_mldivide_d(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'mldivide_d', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% mldivide_d(LDL,b):
% 1) computes D by
%    . extracting the main diagonal entries of LDL
% 2) solves D x = b
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
        error('mldivide_d: in D\\B, LDL must be a square (2D) matrix\n');
    else
        n=double(osizeLDL(1));
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
        error('mldivide: in D\\B, B must be a vector or a (2D) matrix with the same number of rows as U\n');
    end

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsify_mldivide_d (%3d): LDL      size=%-10s, nnz=%d, B  size=%-7s, nnz=%d\n',...
                thisExp,['[',index2str(osizeLDL),']'],length(instrLDL),...
                ['[',index2str(osizeB),']'],length(instrB));
    end

    % Reconstruct instructions of D
    k=find(subsLDL(1,:)==subsLDL(2,:));
    instrD=sparse(double(subsLDL(1,k)),double(subsLDL(2,k)),double(instrLDL(k)),n,n);

    instrB=sparse(double(subsB(1,:)),double(subsB(2,:)),double(instrB),n,m);


    %% Compute instructions
    instrX=instrB;
    if any(diag(instrD)==0)
        error('mldivide_d: zero in the diagonal\n');
    end
    cols=find(any(instrX,1));
    for col=cols
        for row=n:-1:1
            if instrX(row,col)
                operands=full([instrX(row,col),instrD(row,row)]);
                instrX(row,col)=newInstructions(obj,obj.Itypes.I_div,...
                                                {[]},num2cell(operands,2),thisExp);
            end
        end
    end

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
        fprintf('  sparsify_mldivide_d (%3d): U\\B     size=%-10s, nnz=%d,                          # new instr=%4d (%d..%d) (%.2f sec)\n',...
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
    k=find(i==j);
    D=sparse(i(k),j(k),v(k),n,n);

    X=B;
    for col=1:m
        for row=n:-1:1
            X(row,col)=X(row,col)/D(row,row);
        end
    end
    X
    D\B

end
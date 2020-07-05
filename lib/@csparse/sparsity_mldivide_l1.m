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
% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.

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
    
    if verboseLevel<=1 && m*n>20000
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
    for col=1:m
        for row=2:n
            k=find(instrL(row,1:row-1) & instrX(1:row-1,col)');
            if ~isempty(k)
                if instrX(row,col)
                    operands=[instrL(row,k);instrX(k,col)'];
                    operands=full([instrX(row,col),operands(:)']);
                    instrX(row,col)=newInstructions(obj,obj.Itypes.I_plus_minus_dot,...
                                                    {[]},num2cell(operands,2),thisExp);
                elseif instrX(row,col)==0 
                    operands=[instrL(row,k);instrX(k,col)'];
                    operands=full([operands(:)']);
                    instrX(row,col)=newInstructions(obj,obj.Itypes.I_minus_dot,...
                                                    {[]},num2cell(operands,2),thisExp);
                end
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
        fprintf('  sparsify_mldivide_l1(%3d): (I+L)\\B size=%-10s, nnz=%d,                          # new instr=%4d (%d..%d) (%.2f sec)\n',...
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
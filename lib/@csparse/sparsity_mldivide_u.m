function [subsX,instrX]=sparsity_mldivide_u(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'mldivide_u', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% mldivide_u(LU,b):
% 1) computes U by 
%    . extracting the (not strictly) upper-triangular entries of LU
% 2) solves U x = b
% 3) applies any required column-permutation to x 
%    (from the LU factorization)
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

    verboseLevel=0;
    
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);
    
    osizeLU=getOne(obj.vectorizedOperations,'osize',operands(1));
    subsLU=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrLU=getOne(obj.vectorizedOperations,'instructions',operands(1));
    pq=getOne(obj.vectorizedOperations,'parameters',operands(1));
    q=pq{2}; % column permutation
    if length(osizeLU)~=2 || osizeLU(1)~=osizeLU(2)
        osizeLU
        error('mldivide_u: in LU\\B, LU must be a square (2D) matrix\n');
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
        error('mldivide: in U\\B, B must be a vector or a (2D) matrix with the same number of rows as U\n');
    end
    
    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsify_mldivide_u (%3d): LU      size=%-10s, nnz=%d, B  size=%-7s, nnz=%d\n',...
                thisExp,['[',index2str(osizeLU),']'],length(instrLU),...
                ['[',index2str(osizeB),']'],length(instrB));
    end

    % Reconstruct instructions of U
    k=find(subsLU(1,:)<=subsLU(2,:));
    instrU=sparse(double(subsLU(1,k)),double(subsLU(2,k)),double(instrLU(k)),n,n);
    
    instrB=sparse(double(subsB(1,:)),double(subsB(2,:)),double(instrB),n,m);
    
    
    %% Compute instructions
    instrX=instrB;
    for col=1:m
        for row=n:-1:1
            if instrU(row,row)==0
                error('mldivide_u: zero in the diagonal\n');
            end
            k=find(instrU(row,row+1:n) & instrX(row+1:n,col)')+row;
            if instrX(row,col) && ~isempty(k)
                operands=[instrU(row,k);instrX(k,col)'];
                operands=full([instrX(row,col),instrU(row,row),operands(:)']);
                instrX(row,col)=newInstructions(obj,obj.Itypes.I_plus_minus_dot_div,...
                                                {[]},num2cell(operands,2),thisExp);
            elseif instrX(row,col)==0 && ~isempty(k)
                operands=[instrU(row,k);instrX(k,col)'];
                operands=full([instrU(row,row),operands(:)']);
                instrX(row,col)=newInstructions(obj,obj.Itypes.I_minus_dot_div,...
                                                {[]},num2cell(operands,2),thisExp);
            elseif instrX(row,col) && isempty(k)
                operands=full([instrX(row,col),instrU(row,row)]);
                instrX(row,col)=newInstructions(obj,obj.Itypes.I_div,...
                                                {[]},num2cell(operands,2),thisExp);
            end
        end
    end

    % apply q permutation to X
    q(q)=1:length(q);
    instrX=instrX(q,:);

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
        fprintf('  sparsify_mldivide_u (%3d): U\\B     size=%-10s, nnz=%d,                          # new instr=%4d (%d..%d) (%.2f sec)\n',...
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
    k=find(i<=j);
    U=sparse(i(k),j(k),v(k),n,n);
   
    X=B;
    for col=1:m
        for row=n:-1:1
            X(row,col)=(X(row,col)-U(row,row+1:n)*X(row+1:n,col))/U(row,row);
        end
    end
    X
    U\B
    
end
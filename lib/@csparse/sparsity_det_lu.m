function [subsY,instrY]=sparsity_det_lu(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'det_lu', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% det_lu(LU):
% 1) computes det(LU) by 
%    . computing the product of the main diagonal entries of LU
%    . multiplying the result by -1 in case the one of the
%      premutations changes the sign of the det.
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
    p=pq{1}; % row permutation
    q=pq{2}; % col permutation;
    if length(osizeLU)~=2 || osizeLU(1)~=osizeLU(2)
        osizeLU
        error('sparsity_det_lu: in det_lu(LU), LU must be a square (2D) matrix\n');
    else
        n=double(osizeLU(1));
    end
    
    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsity_det_lu (%3d): LU      size=%-10s, nnz=%d\n',...
                thisExp,['[',index2str(osizeLU),']'],length(instrLU));
    end

    % find instructions for diagonal
    k=find(subsLU(1,:)==subsLU(2,:)); % find diagonal elements
    instrX=instrLU(k);
    
    subsY=zeros(0,1);
    if length(instrX)==n
        instrY=newInstruction(obj,obj.Itypes.I_sumprod,[length(instrX),1],instrX(:)',thisExp);
        P(p,1:n)=speye(n);
        Q(q,1:n)=speye(n);
        s=det(P*Q);
        if s<0
            instrY=newInstruction(obj,obj.Itypes.I_sum,[-1],instrY,thisExp);
        end
    else
        error('sparsity_det_lu: LU structurally singular in det_lu(LU)\n');
    end
    
end

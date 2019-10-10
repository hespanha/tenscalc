function [subsY,instrY,instructions]=sparsity_min2(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'plus', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% Copyright 2012-2019 Joao Hespanha

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

osize=getOne(obj.vectorizedOperations,'osize',thisExp);
operands=getOne(obj.vectorizedOperations,'operands',thisExp);

%% Compute sparsity pattern
subsY=zeros(0,length(osize),'uint64');
instrXs=zeros(0,length(operands),'uint64'); % one operand per column
for i=1:length(operands)
    subsX=getOne(obj.vectorizedOperations,'subscripts',operands(i))';
    instrX=getOne(obj.vectorizedOperations,'instructions',operands(i));
    nnzX{i}=length(instrX);
    [liX,kY]=ismember(subsX,subsY,'rows');
    instrXs(kY(liX),i)=instrX(liX);
    subsY=[subsY;subsX(~liX,:)];
    instrXs(end+1:size(subsY,1),i)=instrX(~liX); % grow instrXs
end
[~,k]=sortrows(subsY,size(subsY,2):-1:1);
subsY=subsY(k,:)';
instrXs=instrXs(k,:);

%% Determine instructions for Y
instrY=nan(size(subsY,2),1);
nTerms=sum(instrXs~=0,2);
[instrXs,k]=sort(instrXs,2,'ascend'); % sort terms by instructions

%% Compute instructions
for n=1:max(nTerms)
    k=find(nTerms==n);
    if ~isempty(k)
        instructions=num2cell(instrXs(k,end-n+1:end),2);
        if n<length(operands)
            % less than all terms -> min0 (to account for structural zero)
            instrY(k)=newInstructions(obj,obj.Itypes.I_min0,{[]},instructions,thisExp);
        else
            instrY(k)=newInstructions(obj,obj.Itypes.I_min,{[]},instructions,thisExp);
        end
    end
end

if verboseLevel>0
    fprintf('  sparsify: size=%-10s, nnz=%4d (%6.2f%%) <- min2(',...
            ['[',index2str(osize),']'],length(instrY),100*length(instrY)/prod(osize));
    for i=1:length(operands)
        fprintf('size=%-10s, nnz=%4d; ',['[',index2str(osize),']'],nnzX{i});
    end
    fprintf(')\n');
end


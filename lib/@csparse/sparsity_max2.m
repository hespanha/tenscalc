function [subsY,instrY,instructions]=sparsity_max2(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'plus', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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
    k=find(nTerms==1);
    if ~isempty(k)
        instrY(k)=max(instrXs(k,:),[],2); % get the onky nonzero instruction
    end
    for n=2:max(nTerms)
        k=find(nTerms==n);
        if ~isempty(k)
            instructions=num2cell(instrXs(k,end-n+1:end),2);
            if n<length(operands)
                % less than all terms -> min0 (to account for structural zero)
                instrY(k)=newInstructions(obj,obj.Itypes.I_max0,{[]},instructions,thisExp);
            else
                instrY(k)=newInstructions(obj,obj.Itypes.I_max,{[]},instructions,thisExp);
            end
        end
    end

    if verboseLevel>0
        fprintf('  sparsify: size=%-10s, nnz=%4d (%6.2f%%) <- max2(',...
                ['[',index2str(osize),']'],length(instrY),100*length(instrY)/prod(osize));
        for i=1:length(operands)
            fprintf('size=%-10s, nnz=%4d; ',['[',index2str(osize),']'],nnzX{i});
        end
        fprintf(')\n');
    end

end

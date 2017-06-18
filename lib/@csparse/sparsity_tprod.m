function [subsY,instrY]=sparsity_tprod(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'plus', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
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

if verboseLevel>1
    fprintf('Computing instructions for tprod...\n');
end

t0=clock;

operands=getOne(obj.vectorizedOperations,'operands',thisExp);
parameters=getOne(obj.vectorizedOperations,'parameters',thisExp);
parameters=parameters(2:end); % only keep indices 

%% Compute sparsity pattern
dimY=0;
nSums=0;
for i=1:length(operands)
    osizeX{i}=getOne(obj.vectorizedOperations,'osize',operands(i));
    subsX{i}=getOne(obj.vectorizedOperations,'subscripts',operands(i));
    instrX{i}=getOne(obj.vectorizedOperations,'instructions',operands(i));
    nnzX{i}=length(instrX{i});

    if verboseLevel>1
        fprintf('  Operand %d: size=[%s], nnz=%d, indices=[%s]\n',...
                i,index2str(osizeX{i}),length(instrX{i}),index2str(parameters{i}));
    end
    
    if verboseLevel>2
        fprintf('    Subscripts (pre)\n');
        disp(subsX{i})
    end
    % get subs in the right order for CStprodindices_raw
    % (least-to-most significant index)
    % UNCLEAR IF NEEDED NOW THAT WE DO NOT USE CStprodindices_raw
    if size(subsX{i},1)>=1
        [subsX{i},k]=sortrows(subsX{i}',size(subsX{i},1):-1:1);
        subsX{i}=subsX{i}';
        instrX{i}=instrX{i}(k);
        if verboseLevel>2
            fprintf('  Subscripts (post)\n');
            disp(subsX{i})
        end
    end
    
    dimY=max([dimY,parameters{i}]);
    nSums=max([nSums,-parameters{i}]);
end

dimYS=dimY+nSums;
    
%% Set parameters with all indices starting from 1, starting with the Sums
%% Compute sizes of each indice
indSum=1:nSums;
indOut=nSums+(1:dimY);
tprodSizes=zeros(nSums+dimY,1);
for i=1:length(operands)
    k=find(parameters{i}>0);
    parameters{i}(k)=indOut(parameters{i}(k));
    k=find(parameters{i}<0);
    parameters{i}(k)=indSum(-parameters{i}(k));
    tprodSizes(parameters{i},1)=osizeX{i};
end
osizeY=tprodSizes(nSums+1:end)';

%% Compute # non-zero elements if operand was to be expanded
nnzExpansion=zeros(1,length(operands));
for i=1:length(operands)
    % Which indices/sizes do not appear?
    sz=tprodSizes;
    sz(parameters{i})=[];
    % expansion requires all combinations of missing indices
    nnzExpansion(i)=size(subsX{i},2)*prod(sz);
end

[~,operandOrder]=sort(nnzExpansion,'ascend');

%% Expand "smallest" operand
i=operandOrder(1);
subsYS=zeros(dimYS,nnzExpansion(i));
ind=1:dimYS;
ind(parameters{i})=[];
sz=tprodSizes;
sz(parameters{i})=[];
subs=memory2subscript(sz',1:prod(sz));
subsYS(ind,:)=repmat(subs,1,size(subsX{i},2));
subsYS(parameters{i},:)=kron(double(subsX{i}),ones(1,size(subs,2)));
instrYS=zeros(length(operands),size(subsYS,2));
instrYS(i,:)=kron(double(instrX{i}(:)'),ones(1,size(subs,2)));

%% Expand remaining operators
for i=operandOrder(2:end)
    % Look only within existing expansion
    %[subsX{i}',instrX{i}]
    %subsYS(parameters{i},:)'
    [lia,locb]=ismember(subsYS(parameters{i},:)',subsX{i}','rows');
    subsYS=subsYS(:,lia);
    instrYS=instrYS(:,lia);
    if any(lia)
        instrYS(i,:)=instrX{i}(locb(lia));
    end
end

% YS = matrix with the "expanded" product of all the factors
%      columns: 1:nSums     -- columns corresponding to indices to be summed
%               nSums+1:end -- columns corresponding to the result of tprod

%subsYS
%subsYS(nSums+1:end,:)

% find nonzero entries of tprod
[subsY,ia,ic]=unique(subsYS(nSums+1:end,:)','rows');
subsY=subsY';

if verboseLevel<=1 && length(ia)>10000
    verboseLevel=2;
    fprintf('Computing instructions for (large) tprod with %d nonzero entries... ',length(ia));
end

%% Compute instructions
instrY=nan(size(subsY,2),1);
% ATTENTION: very slow, needs to be sped up
for i=1:length(ia)
    if verboseLevel>1 && mod(i,2000)==0
        fprintf('%d ',i);
    end
    k=find(ic==i);
    if size(instrYS,1)==1 && length(k)==1
        % no need to any sum_prod
        instrY(i)=instrYS(1,k);
    else
        %fprintf('sumprod %d prods, %d sums\n',size(instrYS,1),length(k));
        instrY(i)=newInstruction(obj,obj.Itypes.I_sumprod,...
                                 [size(instrYS,1),length(k)],instrYS(:,k),thisExp);
    end
end

%% Sort subscripts
[subsY,k]=sortrows(subsY',size(subsY,1):-1:1);
subsY=subsY';
instrY=instrY(k);

if verboseLevel>1
    fprintf('done %d nnzYS, %d nnzY (%.2f ms)\n',size(subsYS,2),size(subsY,2),etime(clock,t0));
end

if verboseLevel>0
    fprintf('  sparsify_tprod(%3d): size=%-10s, nnz=%4d (%6.2f%%) <- tprod(',...
            thisExp,['[',index2str(osizeY),']'],length(instrY),100*length(instrY)/prod(osizeY));
    for i=1:length(operands)
        fprintf('size=%-10s, nnz=%4d; ',['[',index2str(osizeX{i}),']'],nnzX{i});
    end
    fprintf(')\n');
end



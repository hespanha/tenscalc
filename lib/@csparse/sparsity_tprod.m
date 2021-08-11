function [subsY,instrY]=sparsity_tprod(obj,thisExp)
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

    if verboseLevel>1
        fprintf('Computing instructions for tprod...\n');
    end

    t0=clock;

    operands=getOne(obj.vectorizedOperations,'operands',thisExp);
    parameters=getOne(obj.vectorizedOperations,'parameters',thisExp);
    parameters=parameters(2:end); % only keep indices

    if verboseLevel>0
        fprintf('  sparsify_tprod(%3d): tprod(',thisExp);
    end

    %% Compute sparsity pattern
    dimY=0;
    nSums=0;
    osizeX=cell(length(operands),1);
    subsX=cell(length(operands),1);
    instrX=cell(length(operands),1);
    nnzX=cell(length(operands),1);
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

        % handle repeated indices:
        % 1) replace operand by diagonal
        % 2) erase repeated indiced from parameters{i}
        uindi=unique(parameters{i});
        if length(uindi)~=length(parameters{i})
            k=1;
            kk=find(parameters{i}(k)==parameters{i},1,'last');
            while kk>k
                % k & kk are equal
                % parameters{i},k,kk

                % find entries with matching indices
                keq=find(subsX{i}(k,:)==subsX{i}(kk,:));
                %subsX{i}(:,keq)

                % order entries (by subscript without the one to be removed)
                [~,kx]=sortrows(subsX{i}([end:-1:kk+1,kk-1:-1:1],keq)');
                keq=keq(kx);

                if verboseLevel>2
                    disp('before');
                    fprintf('  osizeX{i} = ');
                    disp(osizeX{i});
                    fprintf('  parameters{i} = ');
                    disp(parameters{i})
                    fprintf('  subsX{i} =\n');
                    disp(subsX{i});
                    fprintf('  instrX{i} =\n');
                    disp(instrX{i}');
                end

                subsX{i}=subsX{i}(:,keq);
                subsX{i}(kk,:)=[];
                instrX{i}=instrX{i}(keq);
                osizeX{i}(kk)=[];
                parameters{i}(kk)=[];

                if verboseLevel>2
                    disp('after');
                    fprintf('  osizeX{i} = ');
                    disp(osizeX{i});
                    fprintf('  parameters{i} = ');
                    disp(parameters{i})
                    fprintf('  subsX{i} =\n');
                    disp(subsX{i});
                    fprintf('  instrX{i} =\n');
                    disp(instrX{i}');
                end

                % prepare for next one
                k=k+1;
                if k<length(parameters{i})
                    kk=find(parameters{i}(k)==parameters{i},1,'last');
                else
                    kk=-inf;
                end
            end % while
            nnzX{i}=length(instrX{i});

            if verboseLevel>1
                fprintf('  Operand %d: size=[%s], nnz=%d, indices=[%s]\n',...
                        i,index2str(osizeX{i}),length(instrX{i}),index2str(parameters{i}));
            end
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

    if verboseLevel>0
        for i=1:length(operands)
            fprintf('size=%-10s, nnz=%4d; ',['[',index2str(osizeX{i}),']'],nnzX{i});
        end
        fprintf('  ) size=%-10s\n',['[',index2str(osizeY),']']);
    end

    %% Compute # non-zero elements if operand were to be expanded
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

    %% find nonzero entries of tprod
    [subsY,ia,ic]=unique(subsYS(nSums+1:end,:)','rows');
    subsY=subsY';

    if verboseLevel<=1 && length(ia)>20000
        verboseLevel=2;
        fprintf('    computing instructions for (large) tprod with size [%s] and %d nonzero entries... ',index2str(osizeY),length(ia));
    end

    % sort by ic to simplify going over variables with same ic
    if ~isequal(subsYS(nSums+1:end,:),subsY(:,ic)) || ~isequal(subsY(:,:),subsYS(nSums+1:end,ia))
        size(subsYS(nSums+1:end,:)),size(subsY(nSums+1:end,ic)),;
        size(subsY(nSums+1:end,:)),size(subsYS(nSums+1:end,ia)),;
        error('before');
    end
    % before: subsYS=subsY(:,ic), subsY=subsYS(:,ia)
    [ic,k]=sort(ic);
    % after: subsYS(:,k)=subsY(:,ic), subsY=subsYS(:,ia)
    if ~isequal(subsYS(nSums+1:end,k),subsY(:,ic)) || ~isequal(subsY(:,:),subsYS(nSums+1:end,ia))
        error('after');
    end
    subsYS=subsYS(:,k);
    instrYS=instrYS(:,k);

    %% Compute instructions
    instrY=nan(size(subsY,2),1);
    if 0
        % ATTENTION: find(ic==i) is very slow, needs to be sped up
        for i=1:length(ia)
            if verboseLevel>1 && mod(i,5000)==0
                fprintf('%d ',i);
            end
            k=find(ic==i);
            if size(instrYS,1)==1 && length(k)==1
                % no need for any sum_prod
                instrY(i)=instrYS(1,k);
            else
                %fprintf('sumprod %d prods, %d sums\n',size(instrYS,1),length(k));
                instrY(i)=newInstruction(obj,obj.Itypes.I_sumprod,...
                                         [size(instrYS,1),length(k)],instrYS(:,k),thisExp);
            end
        end
    else
        % ATTENTION, requires ic to be sorted (must be done above)
        k1=1;
        for i=1:length(ia)
        if verboseLevel>1 && mod(i,25000)==0
            fprintf('%d ',i);
        end
        k2=k1+1;
        while (k2<=length(ic) && ic(k2)==i); k2=k2+1;end
        % k=find(ic==i);
        % if ~myisequal(k',k1:k2-1)
        %     i, ic(k1:k2-1),
        %     k,k1:k2-1,;
        %     error('mismatch');
        % end
        if size(instrYS,1)==1 && k2==k1+1
            % no need for any sum_prod
            instrY(i)=instrYS(1,k1);
        else
            %fprintf('sumprod %d prods, %d sums\n',size(instrYS,1),length(k1:k2-1));
            instrY(i)=newInstruction(obj,obj.Itypes.I_sumprod,...
                                     [size(instrYS,1),k2-k1],instrYS(:,k1:k2-1),thisExp);
        end
        k1=k2;
        end
    end

    %% Sort subscripts
    [subsY,k]=sortrows(subsY',size(subsY,1):-1:1);
    subsY=subsY';
    instrY=instrY(k);

    if verboseLevel>1
        fprintf('tprod nnzYS=%d, nnzY=%d',size(subsYS,2),size(subsY,2));
    end

    if verboseLevel>0
        fprintf('    -> size=%-10s, nnz=%4d (%6.2f%%) <- tprod(',...
                ['[',index2str(osizeY),']'],length(instrY),100*length(instrY)/prod(osizeY));
        fprintf(')\n');
    end

    if verboseLevel>1
        fprintf('tprod done in %.2f sec  ',etime(clock,t0));
    end
end

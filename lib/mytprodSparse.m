function M=mytprodSparse(tprod_size,varargin)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    verboseLevel=0;

    dimY=0;
    nSums=0;
    for l=1:2:length(varargin)
        i=(l+1)/2;
        osize{i}=size(varargin{l});
        parameters{i}=varargin{l+1};

        if length(osize{i})==2
            % matrix
            [ii,jj,instrX{i}]=find(varargin{l});
            subsX{i}=[ii(:)';jj(:)'];
        else
            % tensor
            ii=find(varargin{l});
            sub=cell(length(osize{i}),1);
            [sub{:}]=ind2sub(osize{i},ii(:));
            subsX{i}=cat(2,sub{:})';
            instrX{i}=nonzeros(varargin{l});
        end
        instrX{i}=instrX{i}(:)';

        if verboseLevel>2
            fprintf('Operand %d (size = [%s], indices = [%s])\n',...
                    i,index2str(osize{i}),index2str(parameters{i}));
            fprintf('  Subscripts (pre) (size=[%s])\n',index2str(size(subsX{i})));
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
                fprintf('  Subscripts (post) (size=[%s])\n',index2str(size(subsX{i})));
                disp(subsX{i})
            end
        end

        if length(osize{i})>length(parameters{i})
            if all(osize{i}(length(parameters{i})+1:end)==1)
                % extra dimensions are singletons
                osize{i}(length(parameters{i})+1:end)=[];
                subsX{i}(length(parameters{i})+1:end,:)=[];
            else
                error('mytprodSparse: size of %d argument ([%s]) does not match size of parameters ([%s])\n',...
                      index2str(osize{i}),index2str(parameters));
            end
        end

        if 0
            fprintf('parameter %d: osize=[%s], parameters=[%s], size(subsX)=[%s], size(instrX)=[%s]\n',...
                    i,index2str(osize{i}),index2str(parameters{i}),...
                    index2str(size(subsX{i})),index2str(size(instrX{i})));
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
    for i=1:length(subsX)
        k=find(parameters{i}>0);
        parameters{i}(k)=indOut(parameters{i}(k));
        k=find(parameters{i}<0);
        parameters{i}(k)=indSum(-parameters{i}(k));
        tprodSizes(parameters{i},1)=osize{i};
    end

    %% Compute # non-zero elements if operand was to be expanded
    nnzExpansion=zeros(1,length(subsX));
    for i=1:length(subsX)
        % Which indices/sizes do not appear?
        sz=tprodSizes;
        sz(parameters{i})=[];
        % expansion requires all combinations of missing indices
        nnzExpansion(i)=size(subsX{i},2)*prod(sz);
    end

    [~,operandOrder]=sort(nnzExpansion,'ascend');

    %% Expand "smallest" operand
    i=operandOrder(1);
    if verboseLevel>3
        fprintf('Expanding smallest operand %d\n',i);
    end
    subsYS=zeros(dimYS,nnzExpansion(i));
    ind=1:dimYS;
    ind(parameters{i})=[];
    sz=tprodSizes;
    sz(parameters{i})=[];
    subs=memory2subscript(sz,1:prod(sz));
    subsYS(ind,:)=repmat(subs,1,size(subsX{i},2));
    subsYS(parameters{i},:)=kron(subsX{i},ones(1,size(subs,2)));
    instrYS=kron(instrX{i}(:)',ones(1,size(subs,2)));

    if verboseLevel>3
        subsYS
        instrYS
        subsYS(nSums+1:end,:)
    end

    %% Expand remaining operators
    for i=operandOrder(2:end)
        if verboseLevel>3
            fprintf('Expanding operand %d\n',i);
        end
        % Look only within existing expansion
        if verboseLevel>3
            [subsX{i}',instrX{i}']
            subsYS(parameters{i},:)'
        end
        [lia,locb]=ismember(subsYS(parameters{i},:)',subsX{i}','rows');
        subsYS=subsYS(:,lia);
        if isempty(subsYS)
            break;
        end
        instrYS=instrYS(:,lia).*instrX{i}(1,locb(lia));
    end

    if verboseLevel>3
        subsYS
        instrYS
        subsYS(nSums+1:end,:)
    end

    while length(tprod_size)<2
        tprod_size(end+1)=1;
    end
    M=zeros(tprod_size);

    if isempty(subsYS)
        %fprintf('EMPTY subsYS\n');
        return
    end

    [subsY,ia,ic]=unique(subsYS(nSums+1:end,:)','rows');
    subsY=subsY';

    instrY=full(sparse(ic,ones(size(ic)),instrYS))';

    ndx=[1,cumprod(tprod_size)];
    ndx(end)=[];
    ndx=ndx*(subsY-1)+1;

    M(ndx)=instrY;

end

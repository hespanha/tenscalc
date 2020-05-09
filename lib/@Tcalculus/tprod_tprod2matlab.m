function obj=tprod_tprod2matlab(obj)
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

    global TCsymbolicExpressions;

    % directly access osize to improve performance
    %old_size=size(obj;
    old_size=obj.osize;

    %% collects operands and convert nested tprods

   
    % directly access operands to improve performance
    %ops=operands(obj);
    ops=TCsymbolicExpressions(obj.TCindex).operands;
    
    objs=cell(length(ops),1);
    changed=false;
    for i=1:length(ops)
        % provide size to speed up object creation
        % objs{i}=Tcalculus(ops(i));
        objs{i}=Tcalculus(ops(i),TCsymbolicExpressions(ops(i)).osize);
        objs{i}=tprod_tprod2matlab(objs{i});
        if ops(i)~=objs{i}.TCindex
            changed=true;
            ops(i)=objs{i}.TCindex;
        end
    end
    
    % directly access type to improve performance
    %if ~strcmp(type(obj),'tprod')
    if ~strcmp(TCsymbolicExpressions(obj.TCindex).type,'tprod')
        if changed
            obj=Tcalculus(type(obj),size(obj),parameters(obj),ops,op_parameters(obj),file_line(obj));
        end
        % directly access osize to improve performance
        %if ~myisequal(size(obj),old_size)
        %if ~myisequal(obj.osize,old_size)
        if length(obj.osize)~=length(old_size) || any(obj.osize~=old_size)
            error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
        end
        return
    end
            
    inds=op_parameters(obj);
    %% sort indices
    [inds,objs]=sortIndices(inds,objs);

    % directly access osize to improve performance
    %tprod_size=size(obj);
    tprod_size=obj.osize;
    sums_size=parameters(obj);
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collapse (positive) dimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % find all possible sequences of (consecutive) indices to collapse
    mink=min(cellfun(@(x)min([x,0]),inds,'UniformOutput',1));
    maxk=max(cellfun(@(x)max([x,0]),inds,'UniformOutput',1));
    potentialCollapse=zeros(0,2);
    if maxk>1
        [i,j]=ind2sub([maxk,maxk],find(triu(true(maxk),1)));
        potentialCollapse=[potentialCollapse;i,j];
        [i,j]=ind2sub([maxk,maxk],find(tril(true(maxk),-1)));
        potentialCollapse=[potentialCollapse;i,j];
    end   
    if mink<-1
        [i,j]=ind2sub([-mink,-mink],find(triu(true(-mink),1)));
        potentialCollapse=[potentialCollapse;-i,-j];
        [i,j]=ind2sub([-mink,-mink],find(tril(true(-mink),-1)));
        potentialCollapse=[potentialCollapse;-i,-j];
    end
    % sort by length
    [~,k]=sort(abs(potentialCollapse(:,2)-potentialCollapse(:,1)),'descend');
    potentialCollapse=potentialCollapse(k,:);
    collapse=[];

    % Check if any can collapse
    for k=1:size(potentialCollapse,1)
        potentialCollapse(k,:);
        if potentialCollapse(k,1)>0
            % always increasing order for positive indices
            if potentialCollapse(k,1)<potentialCollapse(k,2)
                collapse=potentialCollapse(k,1):potentialCollapse(k,2);
            else
                collapse=potentialCollapse(k,2):potentialCollapse(k,1);
            end
        else
            % any order for negative indices (summation)
            if potentialCollapse(k,1)<potentialCollapse(k,2)
                collapse=potentialCollapse(k,1):potentialCollapse(k,2);
            else
                collapse=potentialCollapse(k,1):-1:potentialCollapse(k,2);
            end            
        end
        for i=1:length(inds)
            ndx=inds{i};
            where=min(strfind(ndx,collapse));
            if ~isempty(where)
                % appearing once is okay
                ndx(where:where+length(collapse)-1)=[];
            end
            if ~isempty(intersect(ndx,collapse))
                % but some elements appear elsewhere -- bad
                collapse=[];
                break
            end
        end
        if ~isempty(collapse)
            if verboseLevel>3
                fprintf('tprod_tprod2matlab: can collapse index [%s]: ',index2str(collapse));
                for i=1:length(inds)
                    fprintf('size=[%s],inds=[%s],',index2str(size(objs{i})),index2str(inds{i}));
                end
                fprintf('\n');
            end
            break
        end
    end
    
    if ~isempty(collapse)
        for i=1:length(inds)
            if verboseLevel>4
                fprintf('before collapse (%d):\n',i)
                disp(objs{i},0,1)
                disp(inds{i})
            end
            ik1=strfind(inds{i},collapse);
            if ~isempty(ik1)
                % reshape if sequence appears
                ik2=ik1+length(collapse)-1;
                osizei=size(objs{i});
                objs{i}=reshape(objs{i},[osizei(1:ik1-1),prod(osizei(ik1:ik2)),osizei(ik2+1:end)]);
                if collapse(1)>0
                    inds{i}=[inds{i}(1:ik1-1),min(collapse),inds{i}(ik2+1:end)];
                else
                    inds{i}=[inds{i}(1:ik1-1),max(collapse),inds{i}(ik2+1:end)];
                end
            end
            % rename affected indices
            if collapse(1)>0
                % collapsing positive indices
                k=inds{i}>max(collapse);
                inds{i}(k)=inds{i}(k)-(max(collapse)-min(collapse));
            else
                % collapsing negative indices (summations)
                k=inds{i}<min(collapse);
                inds{i}(k)=inds{i}(k)+(max(collapse)-min(collapse));
            end
            if verboseLevel>4
                fprintf('after collapse (%d):\n',i)
                disp(objs{i},0,1)
                disp(inds{i})
            end
        end
        if verboseLevel>3
            fprintf('  after collapse: tprod(');
            for i=1:length(inds)
                fprintf('size=[%s],inds=[%s],',index2str(size(objs{i})),index2str(inds{i}));
            end
            fprintf(')\n');
        end
        if collapse(1)>0
            % collapsing positive indices
            if verboseLevel>4
                fprintf('  before collapse: tprod_size = [%s], sums_size = [%s]\n',...
                        index2str(tprod_size),index2str(sums_size));
            end
            old_size=size(obj);
            tprod_size=[tprod_size(1:min(collapse)-1),...
                        prod(tprod_size(collapse)),...
                        tprod_size(max(collapse)+1:end)];
            if verboseLevel>4
                fprintf('  after collapse: tprod_size = [%s], sums_size = [%s]\n',...
                        index2str(tprod_size),index2str(sums_size));
            end
            [ops,inds]=Tcalculus.tprod_sort_operands(objs,inds);
            obj=Tcalculus('tprod',tprod_size,sums_size,ops,inds,1);
            obj=tprod_tprod2matlab(obj);
            obj=reshape(obj,old_size);
            % fprintf('end\n');
            % obj
        else
            % collapsing negative indices (summations)
            if verboseLevel>4
                fprintf('  before collapse: tprod_size = [%s], sums_size = [%s]\n',...
                        index2str(tprod_size),index2str(sums_size));
            end
            sums_size=[sums_size(1:min(-collapse)-1),...
                       prod(sums_size(-collapse)),...
                       sums_size(max(-collapse)+1:end)];
            if verboseLevel>4
                fprintf('  after collapse: tprod_size = [%s], sums_size = [%s]\n',...
                        index2str(tprod_size),index2str(sums_size));
            end
            [ops,inds]=Tcalculus.tprod_sort_operands(objs,inds);
            obj=Tcalculus('tprod',tprod_size,sums_size,ops,inds,1);
            obj=tprod_tprod2matlab(obj);
        end
        if ~myisequal(size(obj),old_size)
            error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
        end
        return
    elseif verboseLevel>4
        fprintf('tprod_tprod2matlab: cannot collapse: ');
        for i=1:length(inds)
            fprintf('size=[%s],inds=[%s],',index2str(size(objs{i})),index2str(inds{i}));
        end
        fprintf('\n');
    end

    % only gets here is no change
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % take out repeated indices by using .*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=length(inds)-1:-1:1
        if isequal(inds{i},inds{i+1})
            objs{i}=times(objs{i},objs{i+1},0);
            objs(i+1)=[];
            inds(i+1)=[];
            if verboseLevel>3
                fprintf('tprod_tprod2matlab: replacing repeated indices by .*\n');
                fprintf('  after: tprod(');
                for j=1:length(inds)
                    fprintf('size=[%s],inds=[%s],',index2str(size(objs{j})),index2str(inds{j}));
                end
                fprintf(')\n');
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find trivial tprods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if length(objs)==1 && isequal(inds{1},1:length(inds{1}))
        obj=objs{1};
        if ~myisequal(size(obj),old_size)
            error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
        end
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find multiplication by scalars - LEVEL 1 BLAS (_SCAL)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:length(objs)
        if isempty(inds{i})
            c=objs{i};
            objs(i)=[];
            inds(i)=[];
            if isempty(objs)
                obj=c;
                return
            else
                % should find most economical way to do multiplication by scalar
                if verboseLevel>3
                    fprintf('tprod_tprod2matlab: doing multiplication by scalars at end\nbefore:\n');
                    disp(obj,0,1);
                end
                [ops,inds]=Tcalculus.tprod_sort_operands(objs,inds);
                obj=Tcalculus('tprod',tprod_size,sums_size,ops,inds,1);
                obj=tprod_tprod2matlab(obj);
                obj=mtimes(c,obj,0);
                if verboseLevel>3
                    fprintf('tprod_tprod2matlab: doing multiplication by scalars at end\nafter:\n');
                    disp(obj,0,1);
                end
            end
            if ~myisequal(size(obj),old_size)
                error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
            end
            return
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUMMATION RULES - MEMORY SAVING!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find summation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nsums=min([inds{:}]);
    for sumind=nsums:-1
        ks=[];
        for i=1:length(objs)
            if ismember(sumind,inds{i}) 
                if sum(inds{i}==sumind)==1 && isempty(ks)
                    ks=i;
                    continue
                end
                ks=[];
                break
            end
        end
        if length(ks)==1
            % replace by sum
            %fprintf('tprod_tprod2matlab: replacing index %d by summation in object %d\n',sumind,ks);
            k=find(inds{ks}==sumind);
            objs{ks}=sum(objs{ks},k,0);
            inds{ks}(k)=[];
            % remove summation index
            for j=1:length(objs)
                inds{j}(inds{j}<sumind)=inds{j}(inds{j}<sumind)+1;
            end
            vars=[objs,inds]';
            vars=vars(:);
            obj=tprod_tprod2matlab(tprod(vars{:}));
            if ~myisequal(size(obj),old_size)
                error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
            end
            return
        end
    end        
    
    %disp('tprod_tprod2matlab 1:')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find dot products - LEVEL 1 BLAS (_DOT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nsums=min([inds{:}]);
    for sumind=nsums:-1
        ks=[];
        for i=1:length(objs)
            if ismember(sumind,inds{i}) 
                if sum(inds{i}==sumind)==1 && length(inds{i})==1 && length(ks)<2
                    ks=[ks,i];
                    continue
                end
                ks=[];
                break
            end
        end
        if length(ks)==2
            % replace by dot product
            objs{ks(1)}=mtimes(objs{ks(1)},objs{ks(2)},0);
            inds{ks(1)}=[];
            objs(ks(2))=[];
            inds(ks(2))=[];
            % remove summation index
            for j=1:length(objs)
                inds{j}(inds{j}<sumind)=inds{j}(inds{j}<sumind)+1;
            end
            vars=[objs,inds]';
            vars=vars(:);
            obj=tprod_tprod2matlab(tprod(vars{:}));
            if ~myisequal(size(obj),old_size)
                error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
            end
            return
        end
    end
    
    %disp('tprod_tprod2matlab 2:')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find matrix-vector products - LEVEL 2 BLAS (_TRMV)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i2=[];
    for i1=1:length(objs) 
        if length(inds{i1})==1 && inds{i1}<0
            sumind=inds{i1};
            i2=[];
            for j=1:length(objs)
                if j~=i1 && ismember(sumind,inds{j})
                    if isempty(i2) && length(inds{j})==2
                        i2=j;
                    else
                        i2=[];
                        break;
                    end
                end
            end
            if ~isempty(i2)
                break;
            end
        end
    end
    if ~isempty(i2)
        if verboseLevel>2
            fprintf('matrix-vector product: vector=%d, matrix=%d, sumind=%d\n',i1,i2,sumind);
            fprintf('  before: tprod(');
            for i=1:length(inds)
                fprintf('size=[%s],inds=[%s],',index2str(size(objs{i})),index2str(inds{i}));
            end
            fprintf(')\n');
        end
        % replace vector by product
        if inds{i2}(2)==sumind
            objs{i1}=mtimes(objs{i2},objs{i1},0);
            inds{i1}=inds{i2}(1);
        else
            objs{i1}=mtimes(objs{i2}',objs{i1},0);
            inds{i1}=inds{i2}(2);
        end
        % remove matrix
        objs(i2)=[];
        inds(i2)=[];
        % remove summation index
        for i=1:length(objs)
            inds{i}(inds{i}<sumind)=inds{i}(inds{i}<sumind)+1;
        end
        if verboseLevel>3
            fprintf('  after:  tprod(');
            for i=1:length(inds)
                fprintf('size=[%s],inds=[%s],',index2str(size(objs{i})),index2str(inds{i}));
            end
            fprintf(')\n');
        end
        vars=[objs,inds]';
        vars=vars(:);
        obj=tprod_tprod2matlab(tprod(vars{:}));
        if ~myisequal(size(obj),old_size)
            error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
        end
        return
    end
    
    %disp('tprod_tprod2matlab 3:')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find matrix-matrix products - LEVEL 3 BLAS (_GEMM)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nsums=min([inds{:}]);
    for sumind=nsums:-1
        ks=[];
        for i=1:length(objs)
            if ismember(sumind,inds{i}) 
                if sum(inds{i}==sumind)==1 && length(inds{i})==2 && length(ks)<2
                    ks=[ks,i];
                    continue
                end
                ks=[];
                break
            end
        end
        if length(ks)==2
            % replace by matrix product
            if inds{ks(1)}(2)==sumind && inds{ks(2)}(1)==sumind
                objs{ks(1)}=mtimes(objs{ks(1)},objs{ks(2)},0);
                inds{ks(1)}=[inds{ks(1)}(1),inds{ks(2)}(2)];
            elseif inds{ks(1)}(1)==sumind && inds{ks(2)}(2)==sumind
                objs{ks(1)}=mtimes(objs{ks(2)},objs{ks(1)},0);
                inds{ks(1)}=[inds{ks(2)}(1),inds{ks(1)}(2)];
            elseif inds{ks(1)}(2)==sumind && inds{ks(2)}(2)==sumind
                objs{ks(1)}=mtimes(objs{ks(1)},objs{ks(2)}',0);
                inds{ks(1)}=[inds{ks(1)}(1),inds{ks(2)}(1)];
            else
                objs{ks(1)}=mtimes(objs{ks(1)}',objs{ks(2)},0);
                inds{ks(1)}=[inds{ks(1)}(2),inds{ks(2)}(2)];
            end
            objs(ks(2))=[];
            inds(ks(2))=[];
            % remove summation index
            for j=1:length(objs)
                inds{j}(inds{j}<sumind)=inds{j}(inds{j}<sumind)+1;
            end
            vars=[objs,inds]';
            vars=vars(:);
            obj=tprod_tprod2matlab(tprod(vars{:}));
            if ~myisequal(size(obj),old_size)
                error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
            end
            return
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXPANSION RULES - MEMORY CONSUMING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %disp('tprod_tprod2matlab 4:')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find vector-vector (external) product - a(column) * b(row)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ks=[];
    for i1=1:length(objs)
        for i2=i1+1:length(objs)
            if length(inds{i1})==1 && length(inds{i2})==1 && inds{i1}~=inds{i2}
                ks=[i1,i2];
                break;
            end
        end
        if ~isempty(ks)
            break
        end
    end
    if length(ks)==2
        % replace by product
        objs{ks(1)}=mtimes(reshape(objs{ks(1)},size(objs{ks(1)}),1),...
                           reshape(objs{ks(2)},1,size(objs{ks(2)})),0);
        inds{ks(1)}=[inds{ks(1)},inds{ks(2)}];
        objs(ks(2))=[];
        inds(ks(2))=[];
        vars=[objs,inds]';
        vars=vars(:);
        obj=tprod_tprod2matlab(tprod(vars{:}));
        if ~myisequal(size(obj),old_size)
            error('Size changed [%s] -> [%s]',index2str(old_size),index2str(size(obj)))
        end
        return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CUSTOM LIBRARY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    library={;
    % length(size(obj)), { obj.inds{1}, obje/inds{2}, ...}, function hadle;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 0D-vector output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % y_{}=A_{k2k3k1}(2,2,9) * B_{k3k1}(2,9) * C_{k2k1}(2,9)
        0, {[-2,-3,-1],[-3,-1],[-2,-1]},      '' ;
             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D-vector output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % y_{i}=A_{k1k2i} * B_{k1k2}
        1, {[-1,-2,1],[-1,-2]},      @(A,B) reshape(A,A.size(1)*A.size(2),A.size(3))'*reshape(B,prod(B.size));
        % y_{i}=A_{k2k1i} * B_{k2k1}
        1, {[-2,-1,1],[-2,-1]},      @(A,B) reshape(A,A.size(1)*A.size(2),A.size(3))'*reshape(B,prod(B.size));
        % y_{i}=A_{k1k2i} * B_{k2k1} * C_{k1k2}
        1, {[-1,-2,1],[-2,-1],[-1,-2]},      @(A,B,C) tprod(A,[-1,-2,1],times(B',C,0),[-1,-2]);
        
        % y_{i}=A_{k1k2i}(2,7,54) * B_{k1k2}(2,7) * c_{k2}(7)
        1, {[-1,-2,1],[-1,-2],[-2]},      '' ;

        % y_{i}=A_{k2k3k1}(2,2,9) * B_{k2k1i}(2,9,54) * C_{k3k1}(2,9)
        1, {[-2,-3,-1],[-2,-1,1],[-3,-1]},      '' ;
        % y_{i}=A_{k3k1i}(2,9,54) * B_{k2k3k1}(2,2,9) * C_{k2k1}(2,9)
        1, {[-3,-1,1],[-2,-3,-1],[-2,-1]},      '' ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D-matrix output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % y_{ij}=A_{ji}
        2, {[2,1]},      @(A) A';
        % y_{ij}=A_{ji} * b_{j}
        2, {[2,1],[2]},           @(A,b) timesrep(b,A') ;
        % y_{ij}=A_{ij} * b_{i}
        2, {[1,2],[1]},           @(A,b) timesrep(A,b);
        % y_{ij}=A_{k1i} * B_{k1j} * c_{k1}
        2, {[-1,1],[-1,2],[-1]},      @(A,B,c) mtimes(A',timesrep(B,c),0);
        % y_{ij}=A_{k1j} * B_{ik1} * c_{k1}
        2, {[-1,2],[1,-1],[-1]},      @(A,B,c) mtimes(B,timesrep(A,c));
        
        % y_{ij}=A_{k1k2i} * B_{k1k2j}
        2, {[-1,-2,1],[-1,-2,2]},      @(A,B) reshape(A,A.size(1)*A.size(2),A.size(3))'*reshape(B,B.size(1)*B.size(2),B.size(3));
        % y_{ij}=A_{k2k1i} * B_{k2k1j}
        2, {[-2,-1,1],[-2,-1,2]},      @(A,B) reshape(A,A.size(1)*A.size(2),A.size(3))'*reshape(B,B.size(1)*B.size(2),B.size(3));
        % y_{ij}=A_{k1k2ij} * B_{k1k2}
        2, {[-1,-2,1,2],[-1,-2]},      @(A,B) reshape(reshape(B,prod(B.size),1)'*reshape(A,prod(A.size(1:2)),prod(A.size(3:4))),A.size(3:4));
        % y_{ij}=A_{k3k1j} * B_{k2k1i} * C_{k3k2}
        2, {[-3,-1,2],[-2,-1,1],[-3,-2]},      @(A,B,C) tprod(reshape(C'*reshape(A,A.size(1),A.size(2)*A.size(3)),C.size(2),A.size(2),A.size(3)),[-2,-1,2],B,[-2,-1,1]);
        % y_{ij}=A_{k2k3k1} * B_{k2k3j} * C_{k1i}
        2, {[-2,-3,-1],[-2,-3,2],[-1,1]},      @(A,B,C) tprod(reshape(reshape(A,prod(A.size(1:2)),A.size(3))*C,A.size(1),A.size(2),C.size(2)),[-1,-2,1],B,[-1,-2,2]);
        
        % y_{ij}=A_{ik1j} * B_{ik1}
        2, {[1,-1,2],[1,-1]},      '';
        % y_{ij}=A_{k1k2j}(2,1,54) * B_{ik1k2}(1,2,1)
        2, {[-1,-2,2],[1,-1,-2]},      '' ;
        % y_{ij}=A_{k3k1j}(2,9,54) * B_{k2k3k1}(2,2,9) * C_{k2k1i}(2,9,54)
        2, {[-3,-1,2],[-2,-3,-1],[-2,-1,1]},      '' ;

        % y_{ij}=A_{k1k2i} * B_{k1k2j} * C_{k1k2}
        2, {[-1,-2,1],[-1,-2,2],[-1,-2]},      @(A,B,C) tprod(A,[-1,-2,1],B.*repmat(reshape(C,[C.size,1]),[1,1,B.size(3)]),[-1,-2,2]);

        %% THE NEXT TWO RULES SEEM TO LEAD TO ERRORS
        % y_{ij}=A_{k1ij}(2,2,9) * B_{k1j}(2,9)
        2, {[-1,1,2],[-1,2]},      '';
        %     @(A,B) reshape(...
        %     tprod_matlab(reshape(A,[A.size(1),A.size(2)*A.size(3)]),[-1,1],...
        %           reshape(repmat(reshape(B,[B.size,1]),[1,A.size(2),1]),[A.size(1),A.size(2)*A.size(3)]),[-1,1]),...
        %     A.size(2),A.size(3));
        % y_{ij}=A_{k1k2ij}(2,1,1000,2) * B_{k1k2i}(2,1,1000)
        2, {[-1,-2,1,2],[-1,-2,1]},   '';   
        %     @(A,B) reshape(...
        %     tprod_matlab(reshape(A,A.size(1)*A.size(2),A.size(3)*A.size(4)),[-1,1],...
        %           reshape(repmat(reshape(B,[B.size,1]),[1,1,1,A.size(4)]),A.size(1)*A.size(2),A.size(3)*A.size(4)),[-1,1]),...
        %     A.size(3),A.size(4));
        % y_{ij}=A_{k3k1i}(2,9,54) * B_{k2k3k1}(2,2,9) * C_{k2k1j}(2,9,54)
        2, {[-3,-1,1],[-2,-3,-1],[-2,-1,2]},      '' ;
        % y_{ij}=A_{k1k2i}(2,8,54) * B_{k1k2j}(2,8,54) * c_{k2}(8)
        2, {[-1,-2,1],[-1,-2,2],[-2]},      '' ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D-matrix output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % y_{ijk}=A_{k1jk} * B_{ik1} * c_{i}
        3, {[-1,2,3],[1,-1],[1]},      @(A,B,c) tprod(reshape(B*reshape(A,A.size(1),A.size(2)*A.size(3)),B.size(1),A.size(2),A.size(3)),[1,2,3],c,1);
        % y_{ijk}=A_{ik1j} * B_{ik1k}
        3, {[1,-1,2],[1,-1,3]},      '';
        % y_{ijk}=A_{ik1k} * B_{k1j}
        3, {[1,-1,3],[-1,2]},      '';
        % y_{ijk}=A_{ij}(1,1) * B_{ik}(1,1) * c_{i}(1)
        3, {[1,2],[1,3],[1]},      @(A,B,c)  tprod(timesrep(A,c),[1,2],B,[1,3]);
        % y_{ijk}=A_{ij}(1,1) * B_{ik}(1,1)
        3, {[1,2],[1,3]},      '';

        % y_{ijk}=A_{ik}(2,2) * b_{j}(1000)
        3, {[1,3],[2]},      '' ;
        % y_{ijk}=A_{ik1k2k3k4}(1,2,1,2,1) * B_{k3k4k}(2,1,54) * C_{k1k2j}(2,1,54)
        3, {[1,-1,-2,-3,-4],[-3,-4,3],[-1,-2,2]},      '' ;
        % y_{ijk}=A_{k1j}(1000,2) * B_{k1k}(1000,2) * C_{ik1}(1,1000) * d_{k1}(1000)
        3, {[-1,2],[-1,3],[1,-1],[-1]},      '' ;
        % y_{ijk}=A_{ijk}(2,8,54) * b_{j}(8)
        3, {[1,2,3],[2]},      '' ;
        % y_{ijk}=A_{k1ij}(2,2,9) * B_{k1jk}(2,9,54)
        3, {[-1,1,2],[-1,2,3]},      '' ;
        % y_{ijk}=A_{k1k2ij}(2,1,1000,2) * B_{k1k2ik}(2,1,1000,2)
        3, {[-1,-2,1,2],[-1,-2,1,3]},      '' ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 4D-matrix output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % y_{ijkl}=A_{k1k2ij} * B_{k1k2kl}
        %4, {[-1,-2,1,2],[-1,-2,3,4]},      @(A,B) reshape(tprod(reshape(A,prod(A.size(1:2)),prod(A.size(3:4))),[-1,1],reshape(B,prod(B.size(1:2)),prod(B.size(3:4))),[-1,2]),A.size(3),A.size(4),B.size(3),B.size(4));
        
        % y_{ijkl}=A_{ik}(30,30) * B_{jl}(2,2)
        4, {[1,3],[2,4]},      '';
        % y_{ijkl}=A_{k1jl}(2,11,25) * B_{ik1k}(1,2,25)
        4, {[-1,2,4],[1,-1,3]},      '';

            };

    % disp('tprod_tprod2matlab - start')
    % obj.type
    % obj
    % objs{:}
    % inds{:}

    found=false;
    for i=1:size(library,1)
        if length(tprod_size)==library{i,1} && length(inds)==length(library{i,2})
            found=true;
            trans=zeros(length(inds),1);
            for k=1:length(inds)
                if isequal(inds{k},library{i,2}{k})
                    continue;
                elseif length(inds{k})==2 && isequal(inds{k},library{i,2}{k}(end:-1:1))
                    trans(k)=1;
                else
                    found=false;
                    break
                end
            end
            if found 
                break
            end
        end
    end

    if found
        if ~isempty(library{i,3})
            %i,library{i,3}
            str=sprintf('library{i,3}(');
            for k=1:length(inds)
                if trans(k)
                    str=sprintf('%sobjs{%d}''',str,k);
                else
                    str=sprintf('%sobjs{%d}',str,k);
                end
                if k<length(inds)
                    str=[str,','];
                end
            end
            obj=eval([str,')']);
            obj=tprod_tprod2matlab(obj);
            if ~myisequal(size(obj),old_size)
                error('Size changed [%s] -> [%s] (library entry %d)',...
                      index2str(old_size),index2str(size(obj)),i)
            end
        else
            obj=Tcalculus('tprod_matlab',size(obj),parameters(obj),ops,op_parameters(obj),file_line(obj));
            if ~myisequal(size(obj),old_size)
                error('Size changed [%s] -> [%s] (library entry %d)',...
                      index2str(old_size),index2str(size(obj)),i)
            end
        end
    else
        if verboseLevel>0
            fprintf('ATTENTION: The tensor product:\n        %% y_{%s}=','h'+(1:length(tprod_size)))
            for k=1:length(inds)
                if length(inds{k})>1
                    fprintf('%c_{','A'+k-1);
                else
                    fprintf('%c_{','a'+k-1);
                end
                for j=1:length(inds{k})
                    if inds{k}(j)>0
                        fprintf('%c','h'+inds{k}(j));
                    else
                        fprintf('k%d',-inds{k}(j));
                    end
                end
                fprintf('}(%s)',index2str(objs{k}.size));
                if k<length(inds)
                    fprintf(' * ');
                end
            end
            fprintf('\n        %d, {',length(tprod_size));
            for k=1:length(inds)
                fprintf('[%s]',index2str(inds{k}));
                if k<length(inds)
                    fprintf(',');
                end
            end
            fprintf('},      @(');
            for k=1:length(inds)
                if length(inds{k})>1
                    fprintf('%c','A'+k-1);
                else
                    fprintf('%c','a'+k-1);
                end
                if k<length(inds)
                    fprintf(',');
                end
            end
            fprintf(') ;\n');
            fprintf('  needs to be included in the library\n')
            fprintf('  Brute-force rule used for now.\n\n');
            %disp('paused');pause
        end
        obj=Tcalculus('tprod_matlab',size(obj),parameters(obj),ops,op_parameters(obj),file_line(obj));
    end

end  

function Z=timesrep(X,y)
    if length(y.size)==1
        Z=mtimes(diag(y,0),X,0);
    elseif length(size(X))==1
        Z=mtimes(y,diag(X,0),0);
    else
        y.size
        error('unexpected call to timesrep');
    end
end

function [inds,objs]=sortIndices(inds,objs)
    order=zeros(length(inds),20);
    for i=1:length(inds)
        order(i,1)=-length(inds{i});
        order(i,2:1-order(i,1))=inds{i};
    end
    [~,k]=sortrows(order);
    inds=inds(k);
    objs=objs(k);
end

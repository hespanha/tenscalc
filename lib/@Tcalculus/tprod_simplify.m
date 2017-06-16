function obj=tprod_simplify(obj)
% obj=tprod_simplify(obj)
%
% Simplify tprod object by 
% 1) removing operands of type transpose, ctranspose
% 2) removing operands of type eye, zero
% 3) removing operands of ones
% 4) removing trival tprod
% 5) unflattening tprod: parent tprod with no summations
%                        child tprod with summations for every
%                        operand
% 6) sort operands to simplify tprod matching
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

    if ~strcmp(type(obj),'tprod')
        return
    end

    %% Compile operatos in indices into cell arrays
    ops=operands(obj);
    inds=op_parameters(obj);
    objs=cell(length(ops),1);
    tprod_size=size(obj);
    sums_size=parameters(obj);
    for i=1:length(ops)
        objs{i}=Tcalculus(ops(i));
    end

    msg='';
    
    %% Process transpose() - exchange indices
    for i=length(objs):-1:1
        if ismember(type(objs{i}),{'ctranspose','transpose'})
            objs{i}=Tcalculus(operands(objs{i}));
            inds{i}=inds{i}(end:-1:1);
        end
    end
    
    %% Process zeros(): -> zeros()
    %          eye():   -> remove summations/indices
    i=1;
    while i<=length(objs)
        switch (type(objs{i}))
          case 'zeros'
            obj=Tzeros(tprod_size);
            return
          case 'eye'
            % i
            % objs{i}
            osizei=size(objs{i});
            for j1=1:length(inds{i})
                j2=mod(j1+length(osizei)/2-1,length(osizei))+1;
                ind1=inds{i}(j1);
                ind2=inds{i}(j2);
                if ind1>=ind2
                    continue
                end
                if ind1<0 
                    % one summation (most negative) must disappear,
                    % but first...  add a Tones with the size of the
                    % index removed to make sure index will not
                    % disappear, if possible Tones will be removed
                    % below
                    objs{end+1,1}=Tones(sums_size(-ind1));
                    inds{end+1,1}=ind2;
                    inds=replaceIndex(inds,ind1,ind2);
                    for ind=ind1-1:-1:-length(sums_size)
                        inds=replaceIndex(inds,ind,ind+1);
                    end
                    % objs{i}
                    % inds{i}
                    % sums_size
                    % ind1
                    %fprintf('Will remove %d, %d and sum %d\n',j1,j2,sums_size(-ind1));
                    sums_size(-ind1)=[];
                    % sums_size
                    % remove corresponding dimensions from eye
                    osizei([j1,j2])=[];
                    objs{i}=Teye(osizei);
                    inds{i}([j1,j2])=[];
                    % try again for same eye
                    i=i-1;
                    break
                else
                    % replace everywhere else highest index (ind2) by lowest (ind1)
                    inds=replaceIndex(inds,ind2,ind1);
                    % but keep ind2 in eye()
                    inds{i}(j2)=ind2;
                    if length(inds)>2
                        msg=horzcat(msg,sprintf('tprod: unable to remove product by eye([%s]) with 2 positive indices [%d,%d] - should be replaced by parent tprod to reduce computation\n',index2str(osizei),ind1,ind2));
                    end
                end
            end
        end
        i=i+1;
    end
    % remove "empty" eye()
    for i=length(objs):-1:1
        if strcmp(type(objs{i}),'eye') && isempty(size(objs{i}))
            objs(i)=[];
            inds(i)=[];
        end
    end

    %% Process ones() - remove multiplications by ones
    c=1;                               % may need to multiply by this constant
    repma=nan(1,length(tprod_size));  % may need repmat
    toremove=[];
    for i=length(objs):-1:1
        if strcmp(type(objs{i}),'ones')
            onessize=size(objs{i});
            % loop over ones dimensions
            for j=length(inds{i}):-1:1
                % check if index appears elsewhere
                appear=false;
                for l=1:length(inds)
                    if l~=i && any(inds{l}==inds{i}(j))
                        % yes, can simply remove index from ones
                        appear=true;
                        break
                    end
                end
                if ~appear
                    if inds{i}(j)<0
                        % does not appear and is a summation -> multiplication by scalar constant
                        c=c*onessize(j);
                    else
                        % does not appear and is positive index -> repmat
                        repma(inds{i}(j))=onessize(j);
                    end
                end
            end
            toremove(end+1)=i;
        end
    end
    % remove ones
    objs(toremove)=[];
    inds(toremove)=[];
    % process multiplication by scalar
    if c~=1
        % add multiplication by constant
        objs{end+1}=Tconstant(c,[]);
        inds{end+1}=[];
    end
    % process repmats
    ks=find(~isnan(repma)); % these dimensions are repeated
    repma(isnan(repma))=1; % these dimensions are maintained
    if ~isempty(ks)
        reshap=tprod_size;
        for k=length(ks):-1:1
            % remove ks(k) index
            tprod_size(ks(k))=[];
            if ks(k)<length(tprod_size);
                % old index ks(k)+1 -> ks(k)
                inds=replaceIndex(inds,ks(k)+1,ks(k));
            end
            % remember for reshape;
            reshap(ks(k))=1;
        end
        [ops,inds]=Tcalculus.tprod_sort_operands(objs,inds);
        obj=Tcalculus('tprod',tprod_size,sums_size,ops,inds,1);
        % simplify what is left
        obj=tprod_simplify(obj);
        % add singleton dimensions to permit repmat
        obj=reshape(obj,reshap);
        % do repmat
        obj=repmat(obj,repma);
        return
    end
    
    %% find trivial tprods
    if length(objs)==1 && myisequal(inds{1},1:length(inds{1}))
        obj=objs{1};
        msg='';
        return
    end

    %% Break away operands just with positive indices (no sums) -- not fully tested
    noSums=cellfun(@(x)all(x>0),inds);
    if ~all(noSums) && any(noSums)
        kNoSums=find(noSums);
        if 0
            for i=1:length(kNoSums)
                msg=horzcat(msg,sprintf('indices %d [%s] have no summations, are being unflattened\n',...
                                        i,index2str(inds{kNoSums(i)})));
            end
            [ops,inds]=Tcalculus.tprod_sort_operands(objs,inds);
            Tcalculus('tprod',tprod_size,sums_size,ops,inds,1)
        end
        % outer tprod
        outer_objs=objs(noSums);
        outer_inds=inds(noSums);
        outer_tprod_size=tprod_size;
        outer_sums_size=[];
        
        inner_objs=objs(~noSums);
        inner_inds=inds(~noSums);
        inner_tprod_size=tprod_size;
        inner_sums_size=sums_size;
        
        % check if some indices not in inner tprod
        exist_inner=false(size(tprod_size));
        for i=1:length(inner_inds)
            ks=inner_inds{i};
            ks(ks<0)=[];
            exist_inner(ks)=true;
        end
        % remove indices that do not exists in inner tprod
        inner_tprod_size(~exist_inner)=[];
        old=1:length(exist_inner);
        old(~exist_inner)=[];
        new=1:length(inner_tprod_size);
        for i=1:length(old)
            inner_inds=replaceIndex(inner_inds,old(i),new(i));
        end
        
        % create inner tprod
        [ops,inner_inds]=Tcalculus.tprod_sort_operands(inner_objs,inner_inds);
        obj_inner=Tcalculus('tprod',inner_tprod_size,inner_sums_size,ops,inner_inds,1);
        
        outer_objs{end+1,1}=obj_inner;
        outer_inds{end+1,1}=old;
        objs=outer_objs;
        inds=outer_inds;
        sums_size=outer_sums_size;
    end
    
    % convert objects to TCindex and sort them
    [ops,inds]=Tcalculus.tprod_sort_operands(objs,inds);

    %global TCsymbolicExpressions;
    %bef=length(TCsymbolicExpressions);
    obj=Tcalculus('tprod',tprod_size,sums_size,ops,inds,1);
    % if bef~=length(TCsymbolicExpressions)
    %     disp('tprod_simplify: tprod');
    %     disp(obj,0,1)
    %     disp('tprod_simplify: matlab');
    %     disp(obj,1,1)
    % end

    if ~isempty(msg)
        fprintf(msg)
        disp(obj)
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inds=replaceIndex(inds,old,new)
%fprintf('replacing index %d->%d\n',old,new);
    for j=1:length(inds)
        inds{j}(inds{j}==old)=new;
    end
end


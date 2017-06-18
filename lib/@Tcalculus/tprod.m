function obj=tprod(varargin)
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
    
    if length(varargin)<2
        error('tprod: t least 2 input arguments needed (%d found)\n',length(varargin));
    end
    
    if isequal(varargin{end},'associate')
        associate=1;
        varargin(end)=[];
    elseif isequal(varargin{end},'noassociate')
        associate=0;
        varargin(end)=[];
    else
        associate=1;  % often helps Hessian being symmetric, leading
                      % to smaller code/scratchbook, but is still
                      % slower;
                      % UPDATE: for sparsity_full seems to lead to
                      % same size & faster in latest version with
                      % unflating of tprod

        associate=0;  % often faster because it avoids very time
                      % consuming tprod
    end
    
    if isequal(varargin{end},'distribute')
        distribute=1;
        varargin(end)=[];
    elseif isequal(varargin{end},'nodistribute')
        distribute=0;
        varargin(end)=[];
    else
        distribute=0;
    end
    
    %% get operands
    [tprod_size,sums_size,objs,inds]=tprod_argin2operands(varargin{:});

    if any(sums_size==0)
        obj=Tzeros(tprod_size);
        return
    end
    
    if isempty(objs)
        varargin{:}
        error('trying to create tprod() with no arguments\n');
    end
    
    %% Apply associative rule for tprod(tprod)
    if associate
        % flatten nested tprod's 
        look4tprods=true;
        while look4tprods 
            look4tprods=false;
            for i=1:length(objs)
                if strcmp(type(objs{i}),'tprod')
                    [objs,inds,sums_size]=flattenTprod(objs,inds,sums_size,i);
                    look4tprods=true;
                    break;
                end
            end
        end
        if isempty(objs)
            varargin{:}
            error('associative rule resulted in tprod() with no arguments\n');
        end
    end
    
    %% Apply distributive rule for tprod(plus) 
    if distribute
        error('tprod distributive rule not implemented\n');
        look4pluses=true;
        while look4pluses 
            look4pluses=false;
            for i=1:length(objs)
                if strcmp(objs{i}.type,'plus')
                    [obj,changed]=distributeTprod(obj,i);
                    if changed
                        return
                    end
                end
            end
        end
        if isempty(objs)
            varargin{:}
            error('distributive rule resulted in tprod() with no arguments\n');
        end
    end
    
    %% Create object
    [ops,inds]=Tcalculus.tprod_sort_operands(objs,inds);
    obj=Tcalculus('tprod',tprod_size,sums_size,ops,inds,1);

    %% Simplify tprod
    obj=tprod_simplify(obj);
    updateFile2table(obj,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Expand: tprod(tprod) -> tprod()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [objs,inds,sums_size]=flattenTprod(objs,inds,sums_size,k)

    global TCsymbolicExpressions;

    nsums=length(sums_size);
    % fix indices of object to expanded
    ops_k=operands(objs{k});
    inds_k=op_parameters(objs{k});
    for i=1:length(ops_k)
        for j=1:length(inds_k{i})
            if inds_k{i}(j)>0
                inds_k{i}(j)=inds{k}(inds_k{i}(j));
            else
                % old summation variable now in inner sum
                inds_k{i}(j)=inds_k{i}(j)-nsums;
                sums_size(-inds_k{i}(j))=TCsymbolicExpressions(ops_k(i)).osize(j);
            end
        end
    end
    % convert operands to Tcalculus objects
    objs_k=cell(length(ops_k),1);
    for i=1:length(ops_k)
        objs_k{i}=Tcalculus(ops_k(i));
    end
    % add new objects
    objs=[objs;objs_k];
    inds=[inds;inds_k];
    % remove old objects
    objs(k)=[];
    inds(k)=[];
end  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Expand: tprod(sum) -> sum(tprod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newobj,changed]=distributeTprod(obj,k)
    
    if ~strcmp(obj.type,'tprod')
        error('distributeTprod can only expand a tprod(), not %s\n',obj.type);
    end
    
    if ~strcmp(obj.objs{k}.type,'plus')
        error('distributeTprod can only expand a tprod(), not %s\n',obj.objs{k}.type);
    end
    
    if length(size(obj.objs{k}))>=0 % 3?
        fprintf('\ndistributeTprod: tprod size=[%s], plus size=[%s]\n',index2str(size(obj)),index2str(size(obj.objs{k})));
        
        changed=true;

        objs=cell(length(obj.objs{k}.objs),1);
        for i=1:length(obj.objs{k}.objs)
            oobjs=obj.objs;
            oobjs{k}=obj.objs{k}.objs{i};
            vars=[oobjs,obj.inds]';
            vars=vars(:);
            objs{i}=tprod(vars{:});
        end
        newobj=Tcalculus('plus',size(obj));
        newobj.objs=objs;
        newobj.inds=obj.objs{k}.inds;
        newobj.isStatic=obj.isStatic;
        %fprintf('distributeTprod: before\n%s',spy(obj));
        %fprintf('distributeTprod: after\n%s',spy(newobj));
    else
        changed=false;
        newobj=obj;
    end
end


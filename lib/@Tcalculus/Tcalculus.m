classdef Tcalculus
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
        
    properties (GetAccess = {},SetAccess = 'immutable')
        tprod_implementation='mytprod'; % mytprod, Ctprod, mdsptprod
        update_file_line=false;
    end

    properties (GetAccess = {?csparse},SetAccess = 'immutable')
    %properties (GetAccess = {?csparse},SetAccess = 'private')
        TCindex    % Index of symbolic expression in the table of
                   % symbolic expressions 
        osize      % size for quick access without reading the table
    end

    properties (GetAccess = {?csparse},SetAccess = {})
    end
    
    properties (GetAccess = {},SetAccess = 'private');
    end

    methods (Static)
        function clear()
        % ATTENTION: very dangerous to call this function if TC
        % symbolic variables remain in the work space because the
        % references will point to the wrong variable
            global TCsymbolicExpressions TCsymbolicExpressionsHash
            TCsymbolicExpressions=[];
            TCsymbolicExpressionsHash=[];
        end
            
        function [ops,inds]=tprod_sort_operands(objs,inds)
        % sort: 1st eye, 2nd ones, 3rd constanst, ... , last variables & by TCindex
        % order seems to make a difference in terms of code size (?)
            
            ops=cellfun(@(x)x.TCindex,objs);
            p=[zeros(length(objs),1),ops];
            for i=1:length(objs)
                switch objs{i}.type
                  case {'eye'}
                    p(i,1)=1;
                  case {'ones'}
                    p(i,1)=2;
                  case {'constant'}
                    p(i,1)=3;
                  case {'variable'}
                    p(i,1)=10;
                  otherwise
                    p(i,1)=5;
                end
            end
            [~,k]=sortrows(p);
            inds=inds(k);
            ops=ops(k);
         end
    end
    
    methods
        function obj=Tcalculus(type,osize,parameters,operands,op_parameters,file_line)
        %% Add entry to table TCsymbolicExpressions of symbolic expressions with fields:
        % type       string
        % osize      row vector
        % parameters [] is no parameters
        %            name     for type='variable'
        %            value    for type='constant'
        %            dim      for type='cat'
        %            cell array of fun for type='compose'
        %            sum_size for type='tprod'
        %            sizes    for type='repmat'
        %            method   for type='Xinterpolate'
        %            struct('type',...,'subs',...) for type='subsref'
        %            struct('typical_subscripts',...,'typical_values,...)
        %                for type='chol','lu','ldl' 
        % operands - column array of operands
        % op_parameters - cell array (indices for tprod
        %                             signs for plus)
        % file_line - string with file/line number where created 
        %             when caller is a double the file and line
        %             number is retrieved from dbstack():
        %             0  - file_line from function that called  addTC2table
        %             1  - file_line from caller of function that
        %                  called  addTC2table
        %             etc
        % hash - for faster search
            if nargin==0
                error('Tcalculus cannot be created without arguments (nargin=%d, nargout=%d)',...
                      nargin,nargout);
            end
            
            if nargin==1
                % create object for entry already in TCsymbolicExpressions
                global TCsymbolicExpressions
                if type>length(TCsymbolicExpressions)
                    error('Reference to symbolic variable erased by Tcalculus.clear()\n')
                end
                obj.TCindex=type;
                obj.osize=TCsymbolicExpressions(type).osize;
                return
            end
            
            if obj.update_file_line && ~ischar(file_line)
                st=dbstack();
                if length(st)<file_line+2
                    file_line='command line';
                else
                    file_line=sprintf('file ''%s'', function ''%s'', line number %d',...
                                      st(file_line+2).file,st(file_line+2).name,st(file_line+2).line);
                end
            end

            % add new entry to TCsymbolicExpressions

            global TCsymbolicExpressions TCsymbolicExpressionsHash

            if isempty(osize) 
                % make size uniform to make sure isequal works when
                % comparing size (technically not needed since we
                % use myisequal)
                osize=[];
            end

            if ~isnumeric(operands)
                operands
                error('internal error, operands should be numeric');
            end
                
            if ~isempty(TCsymbolicExpressions) && ~strcmp(class(TCsymbolicExpressions),'struct')
                error('unexpected global variable TCsymbolicExpressions\n')
            end
            ohash=sum(type)+sum(osize)+sum(operands);
            if ~isempty(TCsymbolicExpressions)
                if isequal(type,'variable')
                    % variable with same name already exists?
                    types={TCsymbolicExpressions(:).type};
                    for i=find(strcmp(types,'variable'))
                        if myisequal(TCsymbolicExpressions(i).parameters,parameters)
                            warning('variable ''%s'' [%s] already exists, new variables ''%s'' [%s] will be created\n',...
                                    parameters,index2str(TCsymbolicExpressions(i).osize),parameters,index2str(osize));
                        end
                    end
                end
                for i=find(TCsymbolicExpressionsHash==ohash)
                    if strcmp(TCsymbolicExpressions(i).type,type) && ...
                            myisequal(TCsymbolicExpressions(i).parameters,parameters) && ...
                            myisequal(TCsymbolicExpressions(i).op_parameters,op_parameters) && ...
                            myisequal(TCsymbolicExpressions(i).osize,osize) && ...
                            myisequal(TCsymbolicExpressions(i).operands,operands)
                        %1) file_line need not match (keeps the file_line
                        %   that originally created the expression)
                        %2) derivatives_cache need not match (keep existing)
                        obj.TCindex=i;
                        obj.osize=osize;
                        return
                    end
                end
            end
            TCsymbolicExpressions(end+1).type=type;
            TCsymbolicExpressions(end).osize=osize;
            TCsymbolicExpressions(end).parameters=parameters;
            TCsymbolicExpressions(end).operands=operands;
            TCsymbolicExpressions(end).op_parameters=op_parameters;
            TCsymbolicExpressions(end).file_line=file_line;
            TCsymbolicExpressions(end).derivatives_cache=zeros(0,2);
            TCsymbolicExpressionsHash(end+1)=ohash;
            TCindex=length(TCsymbolicExpressions);
            obj.TCindex=TCindex;
            obj.osize=osize;
        end
        
        function updateFile2table(obj,file_line)
            if obj.update_file_line && ~ischar(file_line)
                st=dbstack();
                if length(st)<file_line+2
                    file_line='command line';
                else
                    file_line=sprintf('file ''%s'', function ''%s'', line number %d',...
                                      st(file_line+2).file,st(file_line+2).name,st(file_line+2).line);
                end
            end

            global TCsymbolicExpressions;
            TCsymbolicExpressions(obj.TCindex).file_line=file_line;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                          Display                         %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function disp(obj,tprod2mat,maxDepth)
            
            if nargin<2
                tprod2mat=false;
            end

            if nargin<3
                maxDepth=inf;
            end

            if tprod2mat
                obj=tprod_tprod2matlab(obj);
            end

            s=str(obj,tprod2mat,maxDepth);
            fprintf('%s\n%s',s);
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                       Get properties                     %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function varargout=size(obj,dim)
            if 0
                global TCsymbolicExpressions
                osize=TCsymbolicExpressions(obj.TCindex).osize;
            else
                osize=obj.osize;
            end
            if nargin>1
                osize=osize(dim);
            end
            if nargout>1
                if nargout>length(osize)
                    error('too many output arguments for size = [%s]',index2str(osize));
                else
                    varargout=num2cell(osize);
                end
            else
                varargout{1}=osize;
            end
        end
        
        function type=type(obj)
            global TCsymbolicExpressions;
            type=TCsymbolicExpressions(obj.TCindex).type;
        end
        
        function parameters=parameters(obj)
            global TCsymbolicExpressions;
            parameters=TCsymbolicExpressions(obj.TCindex).parameters;
        end
        
        function parameters=name(obj)
            global TCsymbolicExpressions;
            if ~strcmp(TCsymbolicExpressions(obj.TCindex).type,'variable')
                disp(obj)
                error('only varibales have names');
            end
            parameters=TCsymbolicExpressions(obj.TCindex).parameters;
        end
        
        function operands=operands(obj)
            global TCsymbolicExpressions;
            operands=TCsymbolicExpressions(obj.TCindex).operands;
        end
        
        function op_parameters=op_parameters(obj)
            global TCsymbolicExpressions;
            op_parameters=TCsymbolicExpressions(obj.TCindex).op_parameters;
        end
        
        function file_line=file_line(obj)
            global TCsymbolicExpressions;
            file_line=TCsymbolicExpressions(obj.TCindex).file_line;
        end
        
        function derivatives_cache=derivatives_cache(obj)
            global TCsymbolicExpressions;
            derivatives_cache=TCsymbolicExpressions(obj.TCindex).derivatives_cache;
        end
        
        function add2derivatives_cache(obj,var,grad)
            global TCsymbolicExpressions;
           TCsymbolicExpressions(obj.TCindex).derivatives_cache(end+1,:)=[var.TCindex,grad.TCindex];
        end
        
        function bool=isempty(obj)
            bool=any(size(obj)==0);
        end
        
        function len=length(obj)
            len=prod(size(obj));
        end

        function msize=msize(obj,dim)
        % matlab compatible size (i.e., always length >=2)
            msize=size(obj);
            while length(msize)<2
                msize(end+1)=1;
            end
            if nargin==2
                msize=msize(dim);
            end
        end
        
        function ind=end(obj,k,n)
            osize=size(obj);
            if n~=length(osize)
                error('mismatch between object length (%d) and indexing length (%d)',length(osize),n);
            end
            ind=osize(k);
        end
            
        function [children,depth]=children(obj,maxDepth)
        % [children,depth]=children(obj,maxDepth)     
        %
        % Determine children of an obj.  
            global TCsymbolicExpressions;

            depth=inf*ones(1,length(TCsymbolicExpressions));
            depth(obj.TCindex)=0;
            order=inf*ones(1,length(TCsymbolicExpressions));
            order(obj.TCindex)=0;
            added=true;
            k=1;
            while added 
                added=false;
                for i=find(isfinite(depth))
                    children=TCsymbolicExpressions(i).operands;
                    for j=1:length(children)
                        if ~isfinite(depth(children(j))) && depth(i)<maxDepth
                            added=true;
                            depth(children(j))=depth(i)+1;
                            order(children(j))=k;
                            k=k+1;
                            break; % break for depth 1st
                        end
                    end
                        if added
                            break; % break for depth 1st
                        end
                end
            end
            children=find(isfinite(depth));
            depth=depth(children);
            order=order(children);
                
            [order,k]=sort(order);
            depth=depth(k);
            children=children(k);
            %[order,children]=sort(order);
            %depth=depth(children);
        end
                    
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   Referencing, rashaping                 %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=subsref(obj1,S)
            % very ineficient since this method is called for every obj.osize, obj.objs, etc.
            if length(S)==1 && isequal(S.type,'()') 
                if ismember(obj1.type,{'plus','tprod','sum','compose','repmat'})
                    %warning('subsref of ''%s'' should be optimized\n',obj1.type)
                end
                osize=size(obj1);
                if length(S.subs)~=length(osize)
                    S.subs
                    osize
                    error(['subsref: mismatch between object length (%d) and ' ...
                           'indexing length (%d)'],length(osize),length(S.subs));
                end
                trivial=true;
                for i=1:length(S.subs)
                    if ischar(S.subs{i}) && S.subs{i}==':'
                        S.subs{i}=1:osize(i);
                    else
                        if any(S.subs{i}<1)
                            S.subs{i}
                            error('subsref: index smaller than 1\n',i)
                        end
                        if any(S.subs{i}>osize(i))
                            S.subs{i}
                            error('subsref: index exceeds tensor dimension (%d)\n',osize(i))
                        end
                        if ~myisequal(S.subs{i}(:)',1:osize(i))
                            trivial=false;
                        end
                    end
                    osize(i)=length(S.subs{i});
                end
                if isequal(obj1.type,'zeros')
                    obj=Tzeros(osize);
                    updateFile2table(obj,1);
                elseif trivial
                    obj=obj1;
                    updateFile2table(obj,1);
                else
                    obj=Tcalculus('subsref',osize,S,obj1.TCindex,{},1);
                end
            else
                obj=builtin('subsref',obj1,S);
            end
        end
        
        function obj=reshape(obj1,varargin)
            if nargin>2
                osize=[varargin{:}];
            elseif nargin==2
                osize=varargin{1};
            else
                error('reshape require more than one argument\n');
            end
            osize1=size(obj1);
            
            if prod(osize)~=prod(osize1)
                error(['reshape cannot change the number of elements ' ...
                       '(from numel([%s])=%d to numel([%s])=%d)\n'],...
                      index2str(osize1),prod(osize1),index2str(osize),prod(osize));
            end

            if isequal(type(obj1),'reshape')
                % for nested reshapes, only last one counts
                obj1=Tcalculus(operands(obj1));
                osize1=size(obj1);
            end
            if isequal(type(obj1),'zeros')
                obj=Tzeros(osize);
                updateFile2table(obj,1);
            elseif isequal(type(obj1),'ones')
                obj=Tones(osize);
                updateFile2table(obj,1);
            elseif myisequal(osize,osize1)
                % trivial reshape
                obj=obj1;
                updateFile2table(obj,1);
            else
                obj=Tcalculus('reshape',osize,[],obj1.TCindex,{},1);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%               Unitary functions/operators                %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=uplus(obj1)
            obj=plus(obj1,nan,1,nan);
            updateFile2table(obj,1);
        end
        
        function obj=uminus(obj1)
            obj=plus(obj1,nan,-1,nan);
            updateFile2table(obj,1);
        end
        
        function obj=norm1(obj1)
            if strcmp(type(obj1),'zeros')
                obj=Tzeros([]);
                updateFile2table(obj,1);
            else
                obj=Tcalculus('norm1',[],[],obj1.TCindex,{},1);
            end
        end
        
        function obj=norm2(obj1)
            if strcmp(type(obj1),'zeros')
                obj=Tzeros([]);
                updateFile2table(obj,1);
            else
                obj=Tcalculus('norm2',[],[],obj1.TCindex,{},1);
            end
        end
        
        function obj=norminf(obj1)
            if strcmp(type(obj1),'zeros')
                obj=Tzeros([]);
                updateFile2table(obj,1);
            else
                obj=Tcalculus('norminf',[],[],obj1.TCindex,{},1);
            end
        end
        
        function obj=full(obj1)
            %obj1=toCalculus(obj1);
            if isequal(type(obj1),'full')
                obj=obj1;
            else
                obj=Tcalculus('full',size(obj1),[],obj1.TCindex,{},1);
            end
        end
        
        function obj=sum(obj1,ind1,sum2tprod)
        % y=sum(x,dim)
        % 
        % sum a tensor x along dimension dim, resulting is a tensor
        % with the size of x, but with dimension dim, removed.
            if nargin<3
                sum2tprod=true;
            end
            if nargin<2
                error('Tcalculus.sum: must include dimension to sum');
            end
            if sum2tprod
                ind=zeros(1,length(size(obj1)));
                ind(ind1)=-1:-1:-length(ind1);
                ind(ind==0)=1:length(size(obj1))-length(ind1);
                obj=tprod(obj1,ind);
                updateFile2table(obj,1);
            else
                ssize=size(obj1);
                ssize(ind1)=[];
                obj=Tcalculus('sum',ssize,ind1,obj1.TCindex,{},1);
            end
        end
        
        function obj=repmat(obj1,sz,repmat2tprod)
            if nargin<3
                repmat2tprod=false;
            end
            if nargin<2
                error('Tcalculus.repmat: must include size multiplier');
            end
            osize1=size(obj1);
            if repmat2tprod
                error('repmat2prod not implemented');
            else
                % resolve 1-d objects: assume column vector
                if length(osize1)==1 && length(sz)==2
                    osize1=[osize1,1];
                end
                if length(osize1)~=length(sz)
                    obj1,sz
                    error(['repmat''s 2nd argument must match size of object ' ...
                           'length (%d~=%d)'],length(sz),length(osize1))
                end
                if strcmp(type(obj1),'zeros') 
                    obj=Tzeros(osize1.*sz);
                    updateFile2table(obj,1);
                    return
                end
                if strcmp(type(obj1),'ones') 
                    obj=Tones(osize1.*sz);
                    updateFile2table(obj,1);
                    return
                end
                if 0
                    obj=Tcalculus('repmat',osize1.*sz,sz,obj1.TCindex,{},1);
                else
                    % convert repmat to subref
                    subs=cell(length(sz),1);
                    for i=1:length(sz)
                        subs{i}=repmat(1:osize1(i),1,sz(i));
                    end
                    S=struct('type','()','subs',{subs});
                    obj=subsref(obj1,S);
                end
            end
        end

        function obj=reduce(obj1,dimension,type)
        % used by min,max,all,any
            osize=size(obj1);
            osize(dimension)=[];
            obj=Tcalculus(type,osize,dimension,obj1.TCindex,{},2);
        end
        
        function obj=min(obj1,dimension)
            obj=reduce(obj1,dimension,'min');
        end
        function obj=max(obj1,dimension)
            obj=reduce(obj1,dimension,'max');
        end
        function obj=all(obj1,dimension)
            obj=reduce(obj1,dimension,'all');
        end
        function obj=any(obj1,dimension)
            obj=reduce(obj1,dimension,'any');
        end
            
        function obj=diag(obj1,diag2tprod)
            if nargin<2
                diag2tprod=false;
            end
            osize1=size(obj1);
            if diag2tprod
                if length(osize1)==1
                    % vector->matrix
                    obj=tprod(obj1,1,Teye([osize1,osize1]),[1,2]);
                    updateFile2table(obj,1);
                elseif length(osize1)==2 && osize1(1)==osize1(2)
                    % matrix->vector
                    obj=tprod(obj1,[1,1]);
                    updateFile2table(obj,1);
                    warning('diag2prod(matrix->vector) not tested');
                else
                    obj1
                    error('ambigous Tcalculus/diag()');
                end
            else
                if length(osize1)==1
                    % vector->matrix
                    obj=Tcalculus('diag',[osize1,osize1],[],obj1.TCindex,{},1);
                elseif length(osize1)==2 && osize1(1)==osize1(2)
                    % matrix->vector
                    obj=Tcalculus('diag',osize1(1),[],obj1.TCindex,{},1);
                else
                    obj1
                    error('ambigous Tcalculus/diag()');
                end
            end
        end
        
        function obj=transpose(obj1)
            obj=ctranspose(obj1);
            updateFile2table(obj,1);
        end
        
        function obj=ctranspose(obj1)
            %obj1=toCalculus(obj1);
            % if length(size(obj1))==1 % convert vectors to matrices
            %     size(obj1)=[size(obj1),1];
            %     obj1
            %     error
            % end
            osize1=size(obj1);
            type1=type(obj1);
            if isempty(osize1)
                % transpose of scalar = same scalar
                obj=obj1;
                return
            end
            if length(osize1)~=2
                obj1
                error('Transpose on ND array is not defined\n');
            end
            if strcmp(type1,'ctranspose')
                obj=Tcalculus(operands(obj1));
            elseif ismember(type1,{'diag','eye'})
                obj=obj1;
            elseif strcmp(type1,'zeros') 
                obj=Tzeros(osize1(end:-1:1));
                updateFile2table(obj,1);
            elseif strcmp(type1,'ones')
                obj=Tones(osize1(end:-1:1));
                updateFile2table(obj,1);
            else
                obj=Tcalculus('ctranspose',osize1(end:-1:1),[],obj1.TCindex,{},1);
            end
        end
        
        function obj=factor(obj1,type,typical_subscripts,typical_values)
            osize1=size(obj1);
            if length(osize1)~=2
                error('%s is only defined for 2D matrices',type);
            end
            if osize1(1)~=osize1(2)
                error('%s is only defined for square matrices',type);
            end
            if nargin<4
                typical_subscripts='';
                typical_values='';
            end
            pars=struct('typical_subscripts',typical_subscripts,...
                        'typical_values',typical_values);
            obj=Tcalculus(type,osize1,pars,obj1.TCindex,{},1);
        end
        
        
        function obj=chol(obj1,varargin)
            obj=factor(obj1,'chol',varargin{:});
        end
        function obj=ldl(obj1,varargin)
        % LDL=ldl(A) - returns the LDL factorization of a symmetric matrix
            obj=factor(obj1,'ldl',varargin{:});
        end
        function obj=ldl_l(obj1)
        % L=ldl_l(LDL(A)) - returns the L matrix of the LDL factorization of a symmetric matrix
            if ~strcmp(type(obj1),'ldl')
                error('%d can only be called for ldl factorizations (not ''%s'')',type(obj1))
            end
            osize1=size(obj1);
            obj=Tcalculus('ldl_l',osize1,[],obj1.TCindex,{},1);
        end
        function obj=ldl_d(obj1)
        % L=ldl_d(ldl(A)) - returns the D matrix of the LDL factorization of a symmetric matrix
            if ~strcmp(type(obj1),'ldl')
                error('%d can only be called for ldl factorizations (not ''%s'')',type(obj1))
            end
            osize1=size(obj1,1);
            obj=Tcalculus('ldl_d',osize1,[],obj1.TCindex,{},1);
        end
        function obj=lu(obj1,varargin)
        % LU=lu(A) - returns the LU factorization of a general matrix
            obj=factor(obj1,'lu',varargin{:});
        end
        function obj=lu_sym(obj1,varargin)
        % LU=lu_sym(A) - returns the LU factorization of a symmetric matrix
            obj=factor(obj1,'lu_sym',varargin{:});
        end
        
        function obj=pptrs(obj1,obj2)
            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);
            if length(osize1)~=2 || osize1(1)~=osize1(2)
                obj1
                error('pptrs: 1st argument must be a 2D square matrix');
            end
            if length(osize2)~=1
                obj2
                error('pptrs: 2nd argument must be a 1D vector');
            end
            if osize1(1)~=osize2(1)
                obj1
                obj2
                error('pptrs: incompatible sizes %d~=%d',osize1(1),osize2(1));
            end
            obj=Tcalculus('pptrs',osize1,[],[obj1.TCindex;obj2.TCindex],{},1);
        end
        
        function obj=det(obj1)
            %obj1=toCalculus(obj1);
            osize1=size(obj1);
            if length(osize1)~=2
                error('det is only defined for 2D matrices');
            end
            if osize1(1)~=osize1(2)
                error('det is only defined for square matrices');
            end
            if osize1(1)==1
                obj=reshape(obj1,[]);
            elseif osize1(1)==2
                a11=subsref(obj1,struct('type','()','subs',{{1,1}}));
                a12=subsref(obj1,struct('type','()','subs',{{1,2}}));
                a21=subsref(obj1,struct('type','()','subs',{{2,1}}));
                a22=subsref(obj1,struct('type','()','subs',{{2,2}}));
                obj=reshape(a11*a22-a12*a21,[]);
                updateFile2table(obj,1);
            elseif osize1(1)==3
                a11=subsref(obj1,struct('type','()','subs',{{1,1}}));
                a12=subsref(obj1,struct('type','()','subs',{{1,2}}));
                a13=subsref(obj1,struct('type','()','subs',{{1,3}}));
                a21=subsref(obj1,struct('type','()','subs',{{2,1}}));
                a22=subsref(obj1,struct('type','()','subs',{{2,2}}));
                a23=subsref(obj1,struct('type','()','subs',{{2,3}}));
                a31=subsref(obj1,struct('type','()','subs',{{3,1}}));
                a32=subsref(obj1,struct('type','()','subs',{{3,2}}));
                a33=subsref(obj1,struct('type','()','subs',{{3,3}}));
                obj=reshape(a11*a22*a33-a11*a23*a32-a12*a21*a33...
                            +a12*a23*a31+a13*a21*a32-a13*a22*a31,[]);
                updateFile2table(obj,1);
            elseif osize1(1)==4
                a11=subsref(obj1,struct('type','()','subs',{{1,1}}));
                a12=subsref(obj1,struct('type','()','subs',{{1,2}}));
                a13=subsref(obj1,struct('type','()','subs',{{1,3}}));
                a14=subsref(obj1,struct('type','()','subs',{{1,4}}));
                a21=subsref(obj1,struct('type','()','subs',{{2,1}}));
                a22=subsref(obj1,struct('type','()','subs',{{2,2}}));
                a23=subsref(obj1,struct('type','()','subs',{{2,3}}));
                a24=subsref(obj1,struct('type','()','subs',{{2,4}}));
                a31=subsref(obj1,struct('type','()','subs',{{3,1}}));
                a32=subsref(obj1,struct('type','()','subs',{{3,2}}));
                a33=subsref(obj1,struct('type','()','subs',{{3,3}}));
                a34=subsref(obj1,struct('type','()','subs',{{3,4}}));
                a41=subsref(obj1,struct('type','()','subs',{{4,1}}));
                a42=subsref(obj1,struct('type','()','subs',{{4,2}}));
                a43=subsref(obj1,struct('type','()','subs',{{4,3}}));
                a44=subsref(obj1,struct('type','()','subs',{{4,4}}));
                obj=reshape(a14*a23*a32*a41-a13*a24*a32*a41-a14*a22*a33*a41...
                            +a12*a24*a33*a41+a13*a22*a34*a41-a12*a23*a34*a41...
                            -a14*a23*a31*a42+a13*a24*a31*a42+a14*a21*a33*a42...
                            -a11*a24*a33*a42-a13*a21*a34*a42+a11*a23*a34*a42...
                            +a14*a22*a31*a43-a12*a24*a31*a43-a14*a21*a32*a43...
                            +a11*a24*a32*a43+a12*a21*a34*a43-a11*a22*a34*a43...
                            -a13*a22*a31*a44+a12*a23*a31*a44+a13*a21*a32*a44...
                            -a11*a23*a32*a44-a12*a21*a33*a44+a11*a22*a33*a44,[]);
                updateFile2table(obj,1);
            else
                obj=Tcalculus('det',[],[],obj1.TCindex,{},1);
            end
        end
            
        %function obj=abs(obj1)
        %    obj=Tcalculus('abs',size(obj1),[],obj1.TCindex,{},1);
        %end

        function obj=compose(obj1,varargin)

            for i=1:length(varargin)
                if ~strcmp(class(varargin{i}),'function_handle')
                    varargin{i}
                    error('compose must take a function handles as the %-th argument\n',i);
                end
            end
            fsize=size(varargin{1}(0));
            % remove singletons at the end
            while length(fsize)>0 && fsize(end)==1
                fsize(end)=[];
            end
            obj=Tcalculus('compose',[size(obj1),fsize],varargin,obj1.TCindex,{},1);
        end
        
        function obj=round(obj1)
            obj=compose(obj1,@(x__)round(x__));
            updateFile2table(obj,1);
        end

        function obj=ceil(obj1)
            obj=compose(obj1,@(x__)ceil(x__));
            updateFile2table(obj,1);
        end

        function obj=floor(obj1)
            obj=compose(obj1,@(x__)floor(x__));
            updateFile2table(obj,1);
        end

        function obj=abs(obj1)
            obj=compose(obj1,@(x__)abs(x__),@(x__)sign(x__),@(x__)zeros(size(x__)));
            updateFile2table(obj,1);
        end

        function obj=exp(obj1)
            obj=compose(obj1,@(x__)exp(x__),@(x__)exp(x__),@(x__)exp(x__));
            updateFile2table(obj,1);
        end
        function obj=log(obj1)
            obj=compose(obj1,@(x__)log(x__),@(x__)1./x__,@(x__)-1./x__.^2);
            updateFile2table(obj,1);
        end
        function obj=reciprocal(obj1)
            obj=compose(obj1,@(x__)1./x__,@(x__)-1./x__.^2,@(x__)2./x__.^3);
            updateFile2table(obj,1);
        end
        function obj=cube(obj1)
            obj=compose(obj1,@(x__)x__.^3,@(x__)3*x__.^2,@(x__)6*x__);
            updateFile2table(obj,1);
        end
        function obj=sqr(obj1)
            obj=compose(obj1,@(x__)x__.^2,@(x__)2*x__,@(x__)2*ones(size(x__)));
            updateFile2table(obj,1);
        end
        function obj=sqrt(obj1)
            obj=compose(obj1,@(x__)sqrt(x__),@(x__).5./sqrt(x__),@(x__)-.25./x__.^1.5);
            updateFile2table(obj,1);
        end
        function obj=cos(obj1)
            obj=compose(obj1,@(x__)cos(x__),@(x__)-sin(x__),@(x__)-cos(x__));
            updateFile2table(obj,1);
        end
        function obj=sin(obj1)
            obj=compose(obj1,@(x__)sin(x__),@(x__)cos(x__),@(x__)-sin(x__));
            updateFile2table(obj,1);
        end
        function obj=tan(obj1)
            obj=compose(obj1,@(x__)tan(x__),@(x__)sec(x__).^2,@(x__)2*sec(x__).^2.*tan(x__));
            updateFile2table(obj,1);
        end
        function obj=atan(obj1)
            obj=compose(obj1,@(x__)atan(x__),@(x__)1./(1+x__.^2),@(x__)-2*x__./(1+x__.^2).^2);
            updateFile2table(obj,1);
        end
        function obj=normpdf(obj1)
            obj=compose(obj1,@(x__)exp(-x__.^2/2)./sqrt(2*pi),...
                              @(x__)-x__.*exp(-x__.^2/2)./sqrt(2*pi),...
                              @(x__)(x__.^2-1).*exp(-x__.^2/2)./sqrt(2*pi));
            updateFile2table(obj,1);
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                 Binary operators                         %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj=plus(obj1,obj2,ind1,ind2)
            
            global TCsymbolicExpressions;

            if nargin<3
                ind1=1;
                ind2=1;
            end
            
            if nargin>1 && ~isnan(ind2)
                [obj1,obj2]=toCalculus(obj1,obj2);
                osize1=size(obj1);
                osize2=size(obj2);

                % following not okay because it changes the size of the result
                % if prod(osize1)==1 && prod(osize2)==1
                %     % singleton + singleton
                %     obj1=reshape(obj1,[]);
                %     osize1=[];
                %     obj2=reshape(obj2,[]);
                %     osize2=[];
                % end
                
                % add scalar? multiply by ones of appropriate size
                if isempty(osize1) && ~isempty(osize2)
                    obj1=growScalar(obj1,osize2);
                    osize1=size(obj1);
                elseif isempty(osize2) && ~isempty(osize1)
                    obj2=growScalar(obj2,osize1);
                    osize2=size(obj2);
                end
                
                if  ~myisequal(osize1,osize2)
                    obj1,obj2
                    error('can only add objects of the same size: [%s] ~= [%s]\n',...
                          index2str(osize1),index2str(osize2))
                end
                
                objs={obj1;obj2};
                inds={ind1;ind2};
            else
                osize1=size(obj1);
                objs={obj1};
                inds={ind1};
            end

            ops=cellfun(@(x)x.TCindex,objs);
            
            %% merge nested sums and remove zeros
            for i=length(ops):-1:1
                if strcmp(TCsymbolicExpressions(ops(i)).type,'zeros')
                    ops(i)=[];
                    inds(i)=[];
                elseif strcmp(TCsymbolicExpressions(ops(i)).type,'plus')
                    ops=[ops;TCsymbolicExpressions(ops(i)).operands];
                    inds=[inds;num2cell(cellfun(@(x)inds{i}*x,TCsymbolicExpressions(ops(i)).op_parameters))];
                    ops(i)=[];
                    inds(i)=[];
                end
            end
            
            %% sort with plus before minus & by TCindex
            [~,k]=sortrows([-[inds{:}]',ops(:)]);
            ops=ops(k);
            inds=inds(k);

            if isempty(ops)
                obj=Tzeros(osize1);
            elseif length(ops)==1 && inds{1}>0
                obj=Tcalculus(ops(1));
            else
                obj=Tcalculus('plus',osize1,[],ops,inds,1);
            end
        end

        function obj=minus(obj1,obj2)
            obj=plus(obj1,obj2,1,-1);
            updateFile2table(obj,1);
        end

        function obj=relational(obj1,obj2,type)
            [obj1,obj2]=toCalculus(obj1,obj2);
            obj3=plus(obj1,obj2,1,-1);
            updateFile2table(obj3,2);
            obj=Tcalculus(type,size(obj3),[],obj3.TCindex,{},2);
        end
            
        function obj=eq(obj1,obj2)
            obj=relational(obj1,obj2,'iszero');  
        end

        function obj=ge(obj1,obj2)
            obj=relational(obj1,obj2,'ispositive');  
        end

        function obj=gt(obj1,obj2)
            obj=relational(obj1,obj2,'ispositive');  
        end

        function obj=le(obj2,obj1)
            obj=relational(obj1,obj2,'ispositive');  
        end

        function obj=lt(obj2,obj1)
            obj=relational(obj1,obj2,'ispositive');  
        end

        function obj=mge(obj1,obj2)
            obj=relational(obj1,obj2,'ispositivedefinite');  
        end
            
        function obj=mle(obj1,obj2)
            obj=relational(obj2,obj1,'ispositivedefinite');  
        end
            
        function obj=times(obj1,obj2,times2tprod)

            if nargin<3
                times2tprod=true;
            end
            
            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);

            if  ~myisequal(osize1,osize2)
                obj1,obj2
                error('can only element-wise multiply objects of the same size: [%s] ~= [%s]\n',...
                      index2str(osize1),index2str(osize2))
            end
            
            if times2tprod
                obj=tprod(obj1,1:length(osize1),obj2,1:length(osize2));
                updateFile2table(obj,1);
            else
                if strcmp(type(obj1),'zeros') || strcmp(type(obj2),'zeros')
                    obj=obj1;
                elseif strcmp(type(obj1),'ones') 
                    obj=obj2;
                elseif strcmp(type(obj2),'ones')
                    obj=obj1;
                else
                    obj=Tcalculus('times',osize1,[],[obj1.TCindex;obj2.TCindex],{},1);
                end
            end
        end

        function obj=rdivide(obj1,obj2)

            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);

            % following not okay because it changes the size of the result
            % if prod(osize1)==1 && prod(osize2)==1
            %     % singleton / singleton
            %     obj1=reshape(obj1,[]);
            %     osize1=[];
            %     obj2=reshape(obj2,[]);
            %     osize2=[];
            % end

            if  ~myisequal(osize1,osize2) && ~isempty(osize2)
                obj1,obj2
                error('can only element-wise divide objects of the same size (or divide by scalar): [%s] ~= [%s]\n',...
                      index2str(osize1),index2str(osize2))
            end
            
            if strcmp(type(obj1),'zeros') 
                obj=obj1;
            else
                obj=Tcalculus('rdivide',osize1,[],[obj1.TCindex;obj2.TCindex],{},1);
            end
        end

        function obj=ldivide(obj1,obj2)
            obj=rdivide(obj2,obj1);
            updateFile2table(obj,1);
        end
        
        function obj=mtimes(obj1,obj2,mtimes2tprod)

            if nargin<3
                mtimes2tprod=true;
            end
            
            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);

            if mtimes2tprod
                if length(osize1)==2 && length(osize2)==2
                    % A * B
                    obj=tprod(obj1,[1,-1],obj2,[-1,2]);
                elseif length(osize1)==2 && length(osize2)==1
                    % A * x
                    obj=tprod(obj1,[1,-1],obj2,[-1]);
                elseif length(osize1)==1 && length(osize2)==2 && osize1==osize2(1)
                    % x' * A 
                    obj=tprod(obj1,[-1],obj2,[-1,1]);
                elseif length(osize1)==1 && length(osize2)==2 && osize2(1)==1
                    % x * A 
                    obj=tprod(reshape(obj1,[osize1,1]),[1,-1],obj2,[-1,2]);
                elseif length(osize1)==1 && length(osize2)==1
                    % x .* y (inner product)
                    obj=tprod(obj1,[-1],obj2,[-1]);
                elseif isempty(osize1)  
                    % scalar * A
                    obj=tprod(obj1,[],obj2,1:length(osize2));
                elseif isempty(osize2)  
                    % A * scalar
                    obj=tprod(obj2,[],obj1,1:length(osize1));
                % following not okay because it changes the size of the result
                % elseif prod(osize1)==1 && prod(osize2)==1
                %     % singleton * singleton
                %     obj=reshape(obj1,[])*reshape(obj2,[]);
                else
                    obj1,obj2,error('Tcalculus/mtimes ambigous')
                end
            else
                if length(osize1)==2 && length(osize2)==2
                    % A * B
                    if osize1(2)~=osize2(1)
                        disp(obj1,0,1)
                        disp(obj2,0,1)
                        error('Inner matrix dimensions must agree: %d ~= %d\n',osize1(2),osize2(1))
                    end
                    osize=[osize1(1),osize2(2)];
                elseif length(osize1)==2 && length(osize2)==1
                    % A * x
                    if osize1(2)~=osize2(1)
                        disp(obj1,0,1)
                        disp(obj2,0,1)
                        error('Inner matrix dimensions must agree: %d ~= %d\n',osize1(2),osize2(1))
                    end
                    osize=osize1(1);
                elseif length(osize1)==1 && length(osize2)==2
                    % x' * A 
                    if osize1(1)~=osize2(1)
                        disp(obj1,0,1)
                        disp(obj2,0,1)
                        error('Inner matrix dimensions must agree: %d ~= %d\n',osize1(1),osize2(1))
                    end
                    osize=osize2(2);
                elseif length(osize1)==1 && length(osize2)==1
                    % x .* y
                    if osize1(1)~=osize2(1)
                        disp(obj1,0,1)
                        disp(obj2,0,1)
                        error('Inner matrix dimensions must agree: %d ~= %d\n',osize1(1),osize2(1))
                    end
                    osize=[];
                elseif isempty(osize1)  
                    % scalar * A
                    osize=osize2;
                elseif isempty(osize2)  
                    % A * scalar
                    osize=osize1;
                else
                    obj1,obj2,error('Tcalculus/mtimes ambigous')
                end

                if strcmp(type(obj1),'zeros') || strcmp(type(obj2),'zeros')
                    obj=Tzeros(msize);
                elseif strcmp(type(obj2),'ones') && myisequal(osize2,[1,1])
                    obj=obj1;
                elseif strcmp(type(obj1),'ones') && myisequal(osize1,[1,1])
                    obj=obj2;
                elseif strcmp(type(obj2),'ones') && osize2(1)==1 && osize2(2)>1 
                    obj=repmat(obj1,[1,osize2(2)],0);
                elseif strcmp(type(obj1),'ones') && osize1(2)==1 && osize1(1)>1 
                    obj=repmat(obj2,[osize1],0); % needed for ones(100,1)*rand(1,36)
                elseif (strcmp(type(obj1),'ones') || strcmp(type(obj2),'ones')) &&...
                        ~isempty(osize1) && ~isempty(osize2)
                    obj=osize1(end)*Tones([osize1(1),osize2(end)]);
                else
                    obj=Tcalculus('mtimes',osize,[],[obj1.TCindex;obj2.TCindex],{},1);
                end
            end
            updateFile2table(obj,1);
        end

        function obj=mrdivide(obj1,obj2)

            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);

            if isempty(osize2)
                obj=rdivide(obj1,obj2);
                updateFile2table(obj,1);
                return
            end
            
            if (length(osize2) ~= 2 && length(osize2) ~= 0)
                obj1,obj2
                error('mrdivide takes scalar/1-vector/2-matrix (not [%s]) and scalar/2-matrix (not [%s])\n',...
                      index2str(osize1),index2str(osize2))
            end
            
            if length(osize2)==2 && osize2(2)~=osize1(1)
                obj1,obj2
                error('Inner matrix dimensions must agree: [%s] ~= [%s]\n',...
                      index2str(osize1),index2str(osize2))
            end
            
            msize=osize1;

            if strcmp(type(obj1),'zeros') 
                obj=Tzeros(msize);
                updateFile2table(obj,1);
            else
                obj=Tcalculus('mrdivide',msize,[],[obj1.TCindex;obj2.TCindex],{},1);
            end
        end

        function obj=mldivide(obj1,obj2)

            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);

            if isempty(osize2)
                obj=ldivide(obj1,obj2);
                updateFile2table(obj,1);
                return
            end
            
            if (length(osize1) ~= 2 && length(osize1) ~= 0) || length(osize2) > 2
                obj1,obj2
                error('mldivide takes 2-matrix (not [%s]) and scalar/1-vector/2-matrix (not [%s])\n',...
                      index2str(osize1),index2str(osize2))
            end
            
            if length(osize1)==2 && osize1(2)~=osize2(1)
                obj1,obj2
                error('Inner matrix dimensions must agree: [%s] ~= [%s]\n',...
                      index2str(osize1),index2str(osize2))
            end
            
            msize=osize2;

            if strcmp(type(obj2),'zeros') 
                obj=Tzeros(msize);
                updateFile2table(obj,1);
            else
                obj=Tcalculus('mldivide',msize,[],[obj1.TCindex;obj2.TCindex],{},1);
            end
        end

        function obj=clp(obj1,obj2)
        % Canonical LP: For a given obj1 with all entries >=0 ,
        % computes the scalar 
        %     max { alpha >0 : obj1 + alpha obj2 >= 0 }
        % obj1 and obj2 must have the same size. 
        % Used by ipm Newton solvers to determine the step size.
            
            [obj1,obj2]=toCalculus(obj1,obj2);
            
            if  ~myisequal(size(obj1),size(obj2))
                obj1,obj2
                error('can only element-wise multiply objects of the same size: [%s] ~= [%s]\n',...
                      index2str(size(obj1)),index2str(size(obj2)))
            end
        
            obj=Tcalculus('clp',[],[],[obj1.TCindex;obj2.TCindex],{},1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                 Multi-ary operators                      %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj=interpolate(obj1,obj2,obj3,obj4,method,convert)
            warning('interpolate not tested');
        
            if nargin<6
                convert=false;
            end
        
            [obj1,obj2,obj3,obj4]=toCalculus(obj1,obj2,obj3,obj4);

            osize1=size(obj1);
            osize2=size(obj2);
            osize3=size(obj3);
            osize4=size(obj4);
            K=osize2(end);
            if osize2(end)~=osize3(end)
                obj2,obj3,
                error('inconsistent lengths for table size in Xi and Yi: %d ~= %d\n',...
                      osize2(end),osize3(end))
            end
            if ~myisequal(osize1,osize2(1:end-1))
                obj1,obj2,
                error('inconsistent size between X and Xi: [%s] ~=[%s]\n',...
                      index2str(osize1),index2str(osize2(1:end-1)))
            end
            osize=osize3(1:end-1);
            if ~isempty(osize4)
                obj4,
                error('S must be a scalar ([%s] instead)\n',index2str(osize4))
            end
            
            if convert
                X=Tvariable('DUMMY_',osize1);
                alpha=length(osize1);
                beta=length(osize3)-1;
                D=obj2-tprod(X,1:alpha,Tones(K),alpha+1);
                D2=tprod(D,[-(1:alpha),1],D,[-(1:alpha),1]);
                ED2=exp(-D2./(2*obj4*obj4));
                F=tprod(obj3,[1:beta,-1],ED2,-1);
                switch method
                  case 'ugaussian'
                    obj=substitute(F,X,obj1);
                  case 'ngaussian'
                    f=sum(ED2,1);
                    obj=substitute(F./f,X,obj1);
                  otherwise
                    error('unknown interpolation method ''%s''\n',method);
                end
                updateFile2table(obj,1);
            else
                obj=Tcalculus('interpolate',osize,method,...
                              [obj1.TCindex;obj2.TCindex;obj3.TCindex;obj4.TCindex],{},1);
            end
        end
        
        function obj=Ginterpolate(obj1,obj2,obj3,obj4,method,convert)
            warning('Ginterpolate not tested');
        
            if nargin<6
                convert=false;
            end
        
            [obj1,obj2,obj3,obj4]=toCalculus(obj1,obj2,obj3,obj4);

            osize1=size(obj1);
            osize2=size(obj2);
            osize3=size(obj3);
            osize4=size(obj4);
            K=osize2(end);
            if osize2(end)~=osize3(end)
                obj2,obj3,
                error('inconsistent lengths for table size in Xi and Yi: %d ~= %d\n',...
                      osize2(end),osize3(end))
            end
            if ~myisequal(osize1,osize2(1:end-1))
                obj1,obj2,
                error('inconsistent size between X and Xi: [%s] ~=[%s]\n',...
                      index2str(osize1),index2str(osize2(1:end-1)))
            end
            osize=[osize3(1:end-1),osize1];
            if ~isempty(osize4)
                obj4,
                error('S must be a scalar ([%s] instead)\n',index2str(osize4))
            end
            
            if convert
                X=Tvariable('DUMMY_',osize1);
                alpha=length(osize1);
                beta=length(osize3(1:end-1));
                D=obj2-tprod(X,1:alpha,Tones(K),alpha+1);
                D2=tprod(D,[-(1:alpha),1],D,[-(1:alpha),1]);
                ED2=exp(-D2./(2*obj4*obj4));
                F=tprod(obj3,[1:beta,-1],ED2,-1);
                GF=gradient(F,X);
                switch method
                  case 'ugaussian'
                    obj=substitute(GF,X,obj1);
                  case 'ngaussian'
                    f=sum(ED2,1);
                    Gf=gradient(f,X);
                    obj=substitute(GF./f-tprod(F,1:beta,Gf,beta+1:beta+alpha)./(f*f),X,obj1);
                  otherwise
                    error('unknown interpolation method ''%s''\n',method);
                end
                updateFile2table(obj,1);
            else
                obj=Tcalculus('Ginterpolate',osize,method,...
                              [obj1.TCindex;obj2.TCindex;obj3.TCindex;obj4.TCindex],{},1);
            end
        end
            
            
        function obj=Hinterpolate(obj1,obj2,obj3,obj4,method,convert)
            warning('Hinterpolate not tested');

            if nargin<6
                convert=false;
            end
        
            [obj1,obj2,obj3,obj4]=toCalculus(obj1,obj2,obj3,obj4);

            osize1=size(obj1);
            osize2=size(obj2);
            osize3=size(obj3);
            osize4=size(obj4);
            K=osize2(end);
            if osize2(end)~=osize3(end)
                obj2,obj3,
                error('inconsistent lengths for table size in Xi and Yi: %d ~= %d\n',...
                      osize2(end),osize3(end))
            end
            if ~myisequal(osize1,osize2(1:end-1))
                obj1,obj2,
                error('inconsistent size between X and Xi: [%s] ~=[%s]\n',...
                      index2str(osize1),index2str(osize2(1:end-1)))
            end
            osize=[osize3(1:end-1),osize1,osize1];
            if ~isempty(osize4)
                obj4,
                error('S must be a scalar ([%s] instead)\n',index2str(osize4))
            end
            
            if convert
                X=Tvariable('DUMMY_',osize1);
                alpha=length(osize1);
                beta=length(osize3(1:end-1));
                D=obj2-tprod(X,1:alpha,Tones(K),alpha+1);
                D2=tprod(D,[-(1:alpha),1],D,[-(1:alpha),1]);
                ED2=exp(-D2./(2*obj4*obj4));
                F=tprod(obj3,[1:beta,-1],ED2,-1);
                GF=gradient(F,X);
                HF=gradient(GF,X);
                switch method
                  case 'ugaussian'
                    obj=substitute(HF,X,obj1);
                  case 'ngaussian'
                    f=sum(ED2,1);
                    Gf=gradient(f,X);
                    Hf=gradient(Gf,X);
                    obj=substitute(...
                        HF./f...
                        -(tprod(GF,1:alpha+beta,Gf,alpha+beta+1:beta+2*alpha)...
                          +tprod(GF,[1:beta,beta+alpha+1:beta+2*alpha],Gf,beta+1:beta+alpha)...
                          +tprod(F,1:beta,Hf,beta+1:beta+2*alpha))./(f*f)...
                        +tprod(F,1:beta,Gf,beta+1:beta+alpha,Gf,beta+alpha+1:beta+2*alpha)*(2./(f*f*f)),...
                        X,obj1);
                  otherwise
                    error('unknown interpolation method ''%s''\n',method);
                end
                updateFile2table(obj,1);
            else
                obj=Tcalculus('Hinterpolate',osize,method,...
                              [obj1.TCindex;obj2.TCindex;obj3.TCindex;obj4.TCindex],{},1);
            end
            
        end
            
            
        % function obj=tprod(varargin) -- in separate file

        function obj=tprod_matlab(varargin)
            [tprod_size,sums_size]=checkTprodSizes(varargin{:});

            if any(sums_size==0)
                obj=Tzeros(tprod_size);
                return
            end
            
            %% Compile operands in 'objs' and 'inds' cells
            objs={};
            inds={};
            for i=1:2:nargin
                obji=toCalculus(varargin{i});
                indi=varargin{i+1};
                objs{end+1,1}=obji;
                inds{end+1,1}=indi;
            end              
            ops=cellfun(@(x)x.TCindex,objs);
            obj=Tcalculus('tprod_matlab',tprod_size,sums_size,ops,inds,1);
        end

        function obj=vertcat(varargin)
            obj=cat(1,varargin{:});
            updateFile2table(obj,1);
        end
        
        function obj=horzcat(varargin)
            obj=cat(2,varargin{:});
            updateFile2table(obj,1);
        end
        
        function obj=cat(dim,varargin)
            [varargin{:}]=toCalculus(varargin{:});
            k=cellfun(@(x)prod(x.msize)==0,varargin);
            varargin(k)=[];
            if length(varargin)==0
                % empty cat
                obj=Tconstant([],0);
                updateFile2table(obj,1);
            elseif length(varargin)==1
                % single input cat
                obj=varargin{1};
                updateFile2table(obj,1);
            else        
                isZeros=true;
                isOnes=true;
                k=[];
                % find number of dimensions
                nd=dim;
                for i=1:length(varargin)
                    obji=varargin{i};
                    sizei=size(obji);
                    if length(sizei)>nd
                        nd=length(sizei);
                    end
                end
                for i=1:length(varargin)
                    obji=varargin{i};
                    typei=type(obji);
                    sizei=size(obji);
                    isZeros=isZeros && isequal(typei,'zeros');
                    isOnes=isOnes && isequal(typei,'ones');
                    % add singleton dimensions to reach desired size
                    if length(sizei)<nd
                        sizei=[sizei,ones(1,nd-length(sizei))];
                        obji=reshape(obji,sizei);
                        updateFile2table(obji,1);
                        varargin{i}=obji;
                    end
                    if i==1
                        osize=sizei;
                        if sizei(dim)>0
                            % only concatenate non-empty arguments
                            k=[k,i];
                        end
                    else
                        if ~myisequal(osize([1:dim-1,dim+1:end]),...
                                      sizei([1:dim-1,dim+1:end]))
                            varargin{:}
                            osize
                            sizei
                            error('incompatible size for concatenation along dimension %d\n',dim);
                        end
                        if sizei(dim)>0
                            % only concatenate non-empty arguments
                            osize(dim)=osize(dim)+sizei(dim);
                            k=[k,i];
                        end
                    end
                end % i=1:length(varargin)
                if isZeros
                    obj=Tzeros(osize);
                    updateFile2table(obj,1);
                elseif isOnes
                    obj=Tones(osize);
                    updateFile2table(obj,1);
                else
                    ops=cellfun(@(x)x.TCindex,varargin(k)');
                    obj=Tcalculus('cat',osize,dim,ops,{},1);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   Derivatives                            %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % function grad=gradient(obj,var) - in separate file
        
        function [hess,grad]=hessian(obj,var1,var2)
        % [hess,grad]=hessian(obj,var1,var2)
        % 
        % Computes the gradient of the expression 'obj', with respect to the 
        % variabe 'var1', and also the hessian of the expression 'obj', with
        % respect to the variables 'var1' and 'var2'.
        % When 'var2' is ommitted, it is assumed equal to 'var1'.

            if nargin<3
                var2=var1;
            end
            grad=gradient(obj,var1);
            updateFile2table(grad,1);
            hess=gradient(grad,var2);
            updateFile2table(hess,1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   Substitution                           %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj=substitute(obj,var,obj1)
        % obj=substitute(obj,var,obj1)
        %
        % Searches for any occurences of the variable 'var' in the expression 'obj'
        % and replaces 'var' by the expression 'obj1'. 'var' and 'obj1' should have the 
        % same size.
        %
        % When 'obj' is a cell array, the substitution is applied to every expression
        % in the cell array.
        %
        % When 'var' is a cell array, obj1 should be a 1-dimensional vector that will
        % be broken into pieces, each with the same size as one of the variables.
        % Each variable is the replaced by the corresponding piece
        % of obj1.

            if iscell(obj)
                % Apply substitution to each element of cell array 'obj'
                for k=1:length(obj)
                    obj{k}=substitute(toCalculus(obj{k}),var,obj1);
                end
                return
            end
            
            if iscell(var)
                % Substitute each variable in call array 'var' by a
                % piece of the 1-d tensor obj1
                if length(size(obj1))~=1
                    var
                    obj1
                    error('Input ''var'' is a cell array, but ''ojb1'' is not a 1-vector\n');
                end
                nx=0;
                for k=1:length(var)
                    osize=size(var{k});
                    xk=reshape(subsref(obj1,struct('type','()',...
                                                   'subs',{{nx+1:nx+prod(osize)}})),...
                               osize);
                    nx=nx+prod(osize);
                    obj=substitute(obj,var{k},xk);
                    updateFile2table(obj,1);
                end
                if nx ~= prod(size(obj1))
                    error('Variables and new expression must have the same total sizes\n');
                end
                return
            end
            
            if isa(obj,'Tcalculus')
                obj=substituteRecursive(obj,var,obj1);
            end

        end
        
        function obj=substituteRecursive(obj,var,obj1)
        % obj=substituteRecursive(obj,var,obj1)
        %
        % Recursive function used by substitute()

            if ~isequal(type(var),'variable')
                var
                error('Can only substitute variables (not %s)',var(type));
            end
            
            if ~myisequal(size(var),size(obj1))
                var
                obj1
                error('Variable and new expression must have the same sizes\n');
            end

            global substituteCounter;
            substituteCounter=substituteCounter+1;

            if strcmp(type(obj),'variable')
                if myisequal(parameters(obj),parameters(var)) && ...
                        myisequal(size(obj),size(var))
                    % variable to be replaced
                    obj=obj1;
                end
            else                
                % replace operands
                ops=operands(obj);
                for i=1:length(ops)
                    op=substituteRecursive(Tcalculus(ops(i)),var,obj1);
                    ops(i)=op.TCindex;
                end
                if ~myisequal(ops,operands(obj))
                    % may need to add expression to TC table
                    obj=Tcalculus(type(obj),size(obj),parameters(obj),...
                                  ops,op_parameters(obj),1);
                end
            end
        end            
                                
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Auxiliary functions (private)                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj2=growScalar(obj1,osize)
% obj2=growScalar(obj1,osize)
%
% given a scalar-valued obj1, returns a tensor with size osize,
% with all entries equal to the scalar
    
    if strcmp(type(obj1),'zeros')
        obj2=Tzeros(osize);
    elseif strcmp(type(obj1),'ones')
        obj2=Tones(osize);
    else
        obj2=reshape(obj1,ones(1,length(osize)));
        updateFile2table(obj2,2);
        S.type='()';
        S.subs=cell(1,length(osize));
        for i=1:length(osize)
            S.subs{i}=ones(osize(i),1);
        end
        obj2=subsref(obj2,S);
        updateFile2table(obj2,2);
    end
end



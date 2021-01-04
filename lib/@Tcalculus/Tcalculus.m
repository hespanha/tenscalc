classdef Tcalculus
% Symbolic tensors class
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

    % properties (GetAccess = {?csparse},SetAccess = {})
    % end

    % properties (GetAccess = {},SetAccess = 'private');
    % end

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
        function obj=Tcalculus(type,osize,parameters,operands,op_parameters,file_line,nowarningsamesize,nowarningever)
        % Add entry to table TCsymbolicExpressions of symbolic expressions with fields:
        % type       string
        % osize      row vector
        % parameters [] is no parameters
        %            name     for type='variable'
        %            value    for type='constant'
        %            dim      for type='cat'
        %            cell array of fun for type='compose'
        %            cell array with {dim,sz,subscripts} for type=vec2tensor
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
        %             when called with a double the file and line
        %             number is retrieved from dbstack():
        %             0  - file_line from function that called  addTC2table
        %             1  - file_line from caller of function that
        %                  called  addTC2table
        %             etc
        % hash - for faster search
            if nargin==2
                % create object for entry already in TCsymbolicExpressions
                % if size is given, no need to access table
                obj.TCindex=type;
                obj.osize=osize;
                return
            end

            global TCsymbolicExpressions

            if ~isempty(TCsymbolicExpressions) && ~isstruct(TCsymbolicExpressions)
                error('unexpected global variable TCsymbolicExpressions\n')
            end

            if nargin==1
                % create object for entry already in TCsymbolicExpressions
                if type>length(TCsymbolicExpressions)
                    error('Reference to symbolic variable erased by Tcalculus.clear()\n')
                end
                obj.TCindex=type;
                obj.osize=TCsymbolicExpressions(type).osize;
                return
            end

            if nargin==0
                error('Tcalculus cannot be created without arguments (nargin=%d, nargout=%d)',...
                      nargin,nargout);
            end

            if nargin<7
                nowarningsamesize=false;
            end

            if nargin<8
                nowarningsever=false;
            end

            if ~ischar(file_line)
                if obj.update_file_line
                    st=dbstack();
                    if length(st)<file_line+2
                        file_line='command line';
                    else
                        file_line=sprintf('file ''%s'', function ''%s'', line number %d',...
                                          st(file_line+2).file,st(file_line+2).name,st(file_line+2).line);
                    end
                else
                    file_line='';
                end
            end

            % add new entry to TCsymbolicExpressions

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

            global TCsymbolicExpressionsHash

            % slower
            %ohash=sum(type)+sum(osize)+prod(osize)+sum(operands)+sum(serialize(op_parameters),'all');
            % faster
            ohash=sum(type)+sum(osize)+prod(osize)+sum(operands);
            if ~isempty(TCsymbolicExpressions)
                if isequal(type,'variable')
                    % variable with same name already exists?
                    types={TCsymbolicExpressions(:).type};
                    for i=find(strcmp(types,'variable'))
                        if myisequal(TCsymbolicExpressions(i).parameters,parameters) && ...
                                ~( nowarningever || (nowarningsamesize && isequal(TCsymbolicExpressions(i).osize,osize)))
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
                        %2) XXX_cache, need not match (keep existing)
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
            TCsymbolicExpressions(end).tprod2matlab_cache=zeros(0,1);
            TCsymbolicExpressions(end).substitute_cache=zeros(0,3);
            TCsymbolicExpressionsHash(end+1)=ohash;
            TCindex=length(TCsymbolicExpressions);
            obj.TCindex=TCindex;
            obj.osize=osize;
        end

        function updateFile2table(obj,file_line)
            if ~ischar(file_line)
                if obj.update_file_line
                    st=dbstack();
                    if length(st)<file_line+2
                        file_line='command line';
                    else
                        file_line=sprintf('file ''%s'', function ''%s'', line number %d',...
                                          st(file_line+2).file,st(file_line+2).name,st(file_line+2).line);
                    end
                else
                    file_line='';
                end
            end

            global TCsymbolicExpressions;
            TCsymbolicExpressions(obj.TCindex).file_line=file_line;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                          Display                         %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function disp(obj,tprod2mat,maxDepth)
        % disp - Display Tcalculus formula
        %   disp(x) displays the formula x
            if nargin<2
                tprod2mat=false;
            end

            if nargin<3
                maxDepth=5;
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

        function bool=isTcalculus(obj)
            bool=true;
        end
        
        function type=type(obj)
        % type - type of the (top) symbolic operation
        %
        %   type(x), for a Tcalculus symbolic tensor, returns the type
        %   of the (top) symbolic operation as a character string.

            global TCsymbolicExpressions;
            type=TCsymbolicExpressions(obj.TCindex).type;
        end

        function parameters=parameters(obj)
            global TCsymbolicExpressions;
            parameters=TCsymbolicExpressions(obj.TCindex).parameters;
        end

        function parameters=name(obj)
        % name - name of a symbolic variable
        %
        %   name(x), for a Tcalculus symbolic variable, returns the
        %   name of the variable as a character string.
            global TCsymbolicExpressions;
            if ~strcmp(TCsymbolicExpressions(obj.TCindex).type,'variable')
                disp(obj)
                error('only variables have names');
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

        function tprod2matlab_cache=tprod2matlab_cache(obj)
            global TCsymbolicExpressions;
            tprod2matlab_cache=TCsymbolicExpressions(obj.TCindex).tprod2matlab_cache;
        end

        function add2tprod2matlab_cache(obj,exp)
            global TCsymbolicExpressions;
            TCsymbolicExpressions(obj.TCindex).tprod2matlab_cache(end+1,:)=exp.TCindex;
        end

        function substitute_cache=substitute_cache(obj)
            global TCsymbolicExpressions;
            substitute_cache=TCsymbolicExpressions(obj.TCindex).substitute_cache;
        end

        function add2substitute_cache(obj,var,expIn,expOut)
            global TCsymbolicExpressions;
            TCsymbolicExpressions(obj.TCindex).substitute_cache(end+1,:)=[var.TCindex,expIn.TCindex,expOut.TCindex];
        end

        function [children,depth]=children(obj,maxDepth)
        % children - returns the children of a Tcalculus object
        %
        %   [children,depth]=children(obj), returns all the children
        %   (i.e. operands) of a Tcalculos symbolic object.

        %   [children,depth]=children(obj,maxDepth), returns the children of a
        %   Tcalculos symbolic object up to a depth determine by
        %   'maxDepth'.

            global TCsymbolicExpressions;

            if nargin<2
                maxDepth=inf;
            end

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
        %%%%                       Get sizes                          %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function varargout=size(obj,dim)
        % size - Size of a tensor
        %
        %   m=size(X), for a tensor X with d dimensions, returns the
        %   size of the tensor X as a row vector with d elements
        %
        %   [m1,m2,...,mk]=size(X), for a tensor X with d dimensions,
        %   returns the size of the tensor X in each dimension (up to
        %   k<=d) as separate variables. An error results if X has
        %   less dimensions than the number of outputs (d<k).
        %
        %   m=size(X,dim), for a tensor X with d dimensions, returns
        %   the number of elements of X in the dimensions specified by
        %   dim. dim can be a scalar or a vector, in the latter case, m
        %   has the same size as d.

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

        function len=ndims(obj)
        % ndims - Number of dimensions in a tensor
        %
        %   ndims(X) returns the number of dimensiona in the tensor X,
        %   defined to be length(size(X)).
        %
        % Note that a scalar (which has size = []) has 0 dimensions

            len=length(size(obj));
        end

        function len=length(obj)
        % length - Length of a tensor
        %
        %   length(X) returns the length of the tensor x, defined to
        %   be max(size(X)).
        %
        % Note that a scalar (which has size = []) has [] length

            len=max(size(obj));
        end

        function len=numel(obj)
        % numel - Number of elements in a tensor
        %
        %   numel(X) returns the number of elements of X, i.e.,
        %   prod(size(X)).
        %
        % Note that a scalar (which has size = []) always has 1
        % element
            len=prod(size(obj));
        end

        function msize=msize(obj,dim)
        % msize - Matlab-compatible size
        %
        %   msize(obj,dim), similar to size(obj,dim), except that it
        %   uses matlab's convention that that scalars and n-vectors
        %   are really 1x1 or nx1 matrices.
            msize=size(obj);
            while length(msize)<2
                msize(end+1)=1;
            end
            if nargin==2
                msize=msize(dim);
            end
            if nargin>1
                msize=msize(dim);
            end
            if nargout>1
                if nargout>length(msize)
                    error('too many output arguments for size = [%s]',index2str(msize));
                else
                    varargout=num2cell(msize);
                end
            else
                varargout{1}=msize;
            end
        end

        function bool=isempty(obj)
        % isempty - true for an empty tensor
        %
        %   isempty(X) returns true if X is a tensor with zero
        %   elements.
        %
        % Note that a scalar (which has size = []) always has 1
        % elements, so it is never empty.
            bool=any(size(obj)==0);
        end

        function ind=end(obj,k,n)
        % see help end
            osize=size(obj);
            if n~=length(osize)
                error('mismatch between object length (%d) and indexing length (%d)',length(osize),n);
            end
            ind=osize(k);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   Referencing, rashaping                 %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function n=numArgumentsFromSubscript(obj,s,indexingContent)
        % see help numArgumentsFromSubscript
            if isequal(s.type,'.')
                n=1;
            else
                error('indexing ''%s'' not valid for Tcalculus objects',s.type);
            end
        end

        function obj=subsref(obj1,S)
        % subsref - subcript reference
        %
        %   X(i,j,k,...) returns a subtensor formed by the elements of
        %   X with subscripts vectors i,j,k,... The resulting tensor
        %   has the same number of dimensions as X, with lengths along
        %   each dimension given by length(i), length(j), length(k),
        %   etc.  A colon ":" can be used as a subscript to indicate
        %   all subscripts on that particular dimensiuon.
        %
        %   subsref(X,S) with a subscript reference structure S
        %   behaves as in regular matlab.
        %
        % Attention: subcripts of the cell type {.} and of the field
        % type .field are not valid for Tcalculus objects.

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
                           'subscript length (%d)'],length(osize),length(S.subs));
                end
                trivial=true;
                for i=1:length(S.subs)
                    if ischar(S.subs{i}) && S.subs{i}==':'
                        S.subs{i}=1:osize(i);
                    else
                        if size(S.subs{i},1)>1 && size(S.subs{i},2)>1
                            S.subs{i}
                            error('subsref: indices must have a single dimension (index %d has size [%s] instead)\n',...
                                  i,index2str(size(S.subs{i})));
                        end
                        if islogical(S.subs{i})
                            if length(S.subs{i})==osize(i)
                                S.subs{i}=find(S.subs{i});
                            else
                                error('subsref: subscript in dimension %d using logical array with length %d that does not match tensor size %d\n',i,length(S.subs{i}),osize(i))
                            end
                        end
                        if any(S.subs{i}<1)
                            S.subs{i}
                            error('subsref: subscript in dimension %d smaller than 1\n',i)
                        end
                        if any(S.subs{i}>osize(i))
                            S.subs{i}
                            error('subsref: subscript in dimension %d exceeds tensor dimension (%d)\n',i,osize(i))
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
        % reshape - reshape tensor
        %
        %   reshape(X,n1,n2,...,nd) or reshape(X,[n1,n2,...,nd])
        %   returns a tensor with size [n1,n2,n3,...] whose elements
        %   are taken from X. An error results in the number of
        %   elements of X is not prod([n1,n2,n3,...]).
        %
        % Attention: unlike matlab's regular reshape,
        % Tcalculus/reshape does not support using [] as one of the
        % dimensions.
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

        function obj=repmat(obj1,varargin)
        % repmat - repicate and tile a tensor
        %
        %   B=repmat(A,n1,n2,...,nd) or B=repmat(A,[n1,n2,...,nd]),
        %   given a tensor A with d dimensions, returns a large tensor
        %   B obtained by tiling copies of A n1 times along the 1st
        %   dimension, n2 times along the 2nd dimension, etc. The
        %   resulting array B has the same number of dimensions as A,
        %   but size
        %
        %         [n1*size(A,1), n2*size(A,2), ... ]
        %

            if nargin>2
                sz=[varargin{:}];
            elseif nargin==2
                sz=varargin{1};
            else
                error('repmat must include size multiplier');
            end
            osize1=size(obj1);
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

        function obj=vec2tensor(obj1,sz,subs,dim)
        % vec2tensor - Expands a vector to a sparse tensor
        %
        %   vec2tensor(X,sz,subs), given
        %     * Tcalculus n-vector X (tensor with size [n]) a
        %     * vector sz with d integers
        %     * nxd matrix subs of subscripts
        %   returns
        %     * Tcalculus tensor Y with size sz, with the nonzero
        %       entries taken from X, with
        %          Y(subs(i,:))=X(i) for i=1:n
        %
        %   vec2tensor(X,sz,subs,dim), given
        %     * Tcalculus tensor X with size(X,dim)=n
        %     * vector sz with d integers
        %     * nxd matrix subs of subscripts
        %   returns
        %     * Tcalculus tensor Y with size similar to that of X, but
        %       the dim dimension expanded to the sizes in sz, and the nonzero
        %       entries taken from X, with
        %          Y(...,subs(i,:),...)=X(...,i,...) for i=1:n
        %       where the ... denote indices of the dimensions before and
        %       after dim
            if nargin<4
                dim=1;
            end

            if ndims(obj1)<1
                error('vec2tensor can only be used for tensors with 1 or more dimensions')
            end

            if length(dim)~=1
                disp(dim)
                error('vec2tensorc''s 3rd argument must be a scalar')
            end

            if dim>ndims(obj1)
                error('vec2tensor''s 3rd argument subs (%d) must be no larger than number of tensor dimensions [%s]\n',...
                      dim,index2str(size(obj1)))
            end

            if length(sz)~=size(subs,2)
                error('vec2tensor: length of 2nd argument sz (%d) must match number of columns of 3rd argument subs ([%s])',...
                      length(sz),index2str(subs));
            end

            osize1=size(obj1);
            if size(subs,1)~=osize1(dim)
                error('vec2tensor: number of rows of 3rd argument subs (size(subs)=[%s]) must match size(X[%s],dim=%d)',...
                      index2str(size(subs)),index2str(osize1),dim);
            end

            if size(subs,1)~=size(unique(subs,'rows'),1)
                error('vec2tensor: the number of rows of subs must be unique (%d rows, but only %d unique\n',...
                      size(subs,1),size(unique(subs,'rows'),1));
            end

            for i=1:length(sz)
                if any(subs(:,i)<1)
                    error('subsref: subscript in dimension %d smaller than 1\n',i);
                end
                if any(subs(:,i)>sz(i))
                    error('vec2tensor: subscript in dimension %d exceeds tensor dimension (%d)\n',i,sz(i))
                end
            end

            osize=[osize1(1:dim-1),sz(:)',osize1(dim+1:end)];
            % subscripts are stores in the transposed form, which is more
            % compatible with tenscalc's internal representation
            obj=Tcalculus('vec2tensor',osize,{dim,sz,subs'},obj1.TCindex,{},1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%               Unitary functions/operators                %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=uplus(obj1)
        % + unary plus
        %
        %   +X for tensors is just X
            %obj=plus(obj1,nan,1,nan);
            obj=obj1;
            updateFile2table(obj,1);
        end

        function obj=uminus(obj1)
        % - unary minus
        %
        %   -X returns the symmetric of X
            obj=plus(obj1,nan,-1,nan);
            updateFile2table(obj,1);
        end

        function obj=norm(obj1,p)
        % norm - Vector, matrix, or tensor norm
        %
        %    norm(x) or norm(x,2) returns the 2-norm of the vector x
        %
        %    norm(x,1) returns the 1-norm of the vector or matrix x
        %
        %    norm(x,inf) returns the infinity norm of the vector or matrix x
        %
        %    norm(x,'fro') returns the Frobenius norm of the tensor x
        %
        % Attention:
        %
        % The 1-norm should be avoided in optimization criteria
        % because it is not differentiable at points where the optimum
        % often lies. This can be done by introducing slack variables
        % and constraints. For example, a minimization of the form
        %      min { norm(x,1) : x in R^n , F(x)>= 0 \}
        % can be reformulates as
        %      min { sum(v) : v,x in R^n , -v <=x<= v, F(x)>= 0 \}
        %
        % The infinity norm should also be avoided in optimization
        % criteria for similar reasons. This can be done by
        % introducing slack variables and constraints. For example, a
        % minimization of the form
        %      min { norm(x,inf) : x in R^n , F(x)>= 0 \}
        % can be reformulates as
        %      min { v : v in R ,x in R^n , -v <=x<= v, F(x)>= 0 \}

            if nargin<2
                p=2;
            end

            if strcmp(type(obj1),'zeros')
                obj=Tzeros([]);
                updateFile2table(obj,1);
            end

            switch p
              case 2
                switch ndims(obj1)
                  case 0
                    obj=abs(obj1);
                  case 1
                    obj=sqrt(norm2(obj1));
                    %obj=Tcalculus('normEuclidean',[],[],obj1.TCindex,{},1);
                  otherwise
                    error('2-norm only available for Tcalculus vectors, not [%s]-tensors',index2str(size(obj1)));
                end
              case 1
                switch ndims(obj1)
                  case 0;
                    obj=abs(obj1);
                  case 1
                    obj=sum(abs(obj1),1);
                    %obj=Tcalculus('norm1',[],[],obj1.TCindex,{},1);
                  case 2
                    obj=max(sum(abs(obj1),1),[],1);
                  otherwise
                    error('1-norm only available for Tcalculus vectors and matrices, not [%s]-tensors',index2str(size(obj1)));
                end
              case inf
                switch ndims(obj1)
                  case 0;
                    obj=abs(obj1);
                  case 1
                    obj=max(abs(obj1));
                    %obj=Tcalculus('norminf',[],[],obj1.TCindex,{},1);
                  case 2
                    obj=max(sum(abs(obj1),2),[],1);
                  otherwise
                    error('inf-norm only available for Tcalculus vectors and matrices, not [%s]-tensors',index2str(size(obj1)));
                end
              case 'fro'
                obj=sqrt(norm2(obj1));
              otherwise
                if isnumeric(p) && length(p)==1
                    if ndims(obj1)==1
                        obj=power(sum(abs(obj1).^p,1),1/p);
                    else
                        error('p-norm only available for Tcalculus vectors, not [%s]-tensors',index2str(size(obj1)));
                    end
                else
                    error('second argument for norm(.) must be 1,2,inf,''fro'',or a scalar');
                end
            end

        end

        function obj=norm1(obj1)
        % norm1 - Sum of absolute values of tensor entries
        %
        %   norm1(x) returns the sum of the absolute value of all
        %   entries of the tensor x, which for vectors corresponds to
        %   the 1-norm of x.
        %
        %   norm1(x) is similar to sum(abs(x),'all'), but the former
        %   computes derivatives more efficiently.
        %
        % Attention:
        %
        % norm1 should be avoided in optimization criteria because
        % it is not differentiable at points where the optimum often
        % lies. This can be done by introducing slack variables and
        % constraints. For example, a minimization of the form
        %      min { norm1(x) : x in R^n , F(x)>= 0 \}
        % can be reformulates as
        %      min { sum(v,'all') : v,x in R^n , -v <=x<= v, F(x)>= 0 \}

            if strcmp(type(obj1),'zeros')
                obj=Tzeros([]);
                updateFile2table(obj,1);
            else
                obj=sum(abs(obj1),'all');
                %obj=Tcalculus('norm1',[],[],obj1.TCindex,{},1);
            end
        end

        function obj=norm2(obj1,S)
        % norm2 - Squared quadratic norm
        %
        %   norm2(x) returns the sum of the square of all entries of
        %   the tensor x, which for matrices corresponds to the square
        %   of the Frobenius norm of x.
        %
        %   norm2(x) is similar to sum(x.^2,'all'), but the former
        %   computes derivatives more efficiently.
        %
        %   norm2(x,S) returns the value of the quadratic form
        %   <x,Sx>. This form is only applicable when x is a vector
        %   (tensor with 1 dimension) and S a square matrix
        %   (tensor with 2 dimensions).

            if strcmp(type(obj1),'zeros')
                obj=Tzeros([]);
                updateFile2table(obj,1);
            else
                if nargin<2
                    obj=Tcalculus('norm2',[],[],obj1.TCindex,{},1);
                else
                    mx=size(obj1);
                    ms=size(S);
                    if length(mx)==1 && length(ms)==2 && ms(1)==mx(1) && ms(2)==mx(1)
                        obj=tprod(obj1,-1,S,[-1,-2],obj1,-2);
                    else
                        error('norm2: norm2(x,S) called with x[%s] not a column vector and/or S[%s] not a square matrix of compatible size\n',index2str(size(obj1)),index2str(size(S)));
                    end
                end
            end
        end

        function obj=norminf(obj1)
        % norminf - Maximum absolute value of tensor entries
        %
        %   norminf(x) returns the largest absolute value of all
        %   entries of the tensor x, which for vectors corresponds to
        %   the infinity-norm of x.
        %
        %   norminf(x) is similar to max(abs(x),[],'all'), but the former
        %   computes derivatives more efficiently.
        %
        % Attention:
        %
        % norminf() should be avoided in optimization criteria because
        % it is not differentiable at points where the optimum often
        % lies. This can be done by introducing slack variables and
        % constraints. For example, a minimization of the form
        %      min { norminf(x) : x in R^n , F(x)>= 0 \}
        % can be reformulates as
        %      min { v : v in R,x in R^n , -v <=x<= v, F(x)>= 0 \}
            if strcmp(type(obj1),'zeros')
                obj=Tzeros([]);
                updateFile2table(obj,1);
            else
                obj=Tcalculus('norminf',[],[],obj1.TCindex,{},1);
            end
        end

        function obj=full(obj1)
        % full - convert sparse tensor to full
        %
        %   Y=full(X) converts a sparse tensor X to a full tensor
        %   Y. No entry of the tensor is changed, but subsequent
        %   computations will not take advantage of the information
        %   that some entries are know to be zero.
        %
        %   The main use of this function is to force the results of a
        %   computation to be returned in a linear memory structure
        %   with all the zeros filled in appropriately. This is
        %   generally required to return n-dimensional tensors to
        %   matlab with n>2, since matlab does not support sparse
        %   arrays with more than 2 dimensions.


            if isequal(type(obj1),'full')
                obj=obj1;
            else
                obj=Tcalculus('full',size(obj1),[],obj1.TCindex,{},1);
            end
        end

        function obj=sum(obj1,dimension,sum2tprod)
        % sum - sum of tensor entries
        %
        %   s=sum(X,vecdim) sums the entries of the tenasor X along
        %   the dimensions specified in vecdim, resulting in a tensor
        %   with the same size as X, but with the dimensions in vecdim
        %   removed.
        %
        %   s=sum(X,'all') sums the entries of X along all its
        %   dimensions, resulting in a scalar. Performs similarly to
        %   sum(X,1:ndims(X))

            if nargin<3
                sum2tprod=true;
            end
            if nargin<2
                switch length(size(obj1))
                  case 0
                    dimension=[];
                  case 1
                    dimension=1;
                  otherwise
                    error('Tcalculus.sum: must include dimension to sum');
                end
            end
            if isempty(dimension)
                % no need to sum
                obj1=obj1;
            end
            if myisequal(dimension,'all')
                dimension=1:ndims(obj1);
            end
            if sum2tprod
                ind=zeros(1,length(size(obj1)));
                ind(dimension)=-1:-1:-length(dimension);
                ind(ind==0)=1:length(size(obj1))-length(dimension);
                obj=tprod(obj1,ind);
                updateFile2table(obj,1);
            else
                ssize=size(obj1);
                ssize(dimension)=[];
                obj=Tcalculus('sum',ssize,dimension,obj1.TCindex,{},1);
            end
        end

        function obj=reduce(obj1,dimension,type)
        % used by min,max,all,any
            if myisequal(dimension,'all')
                dimension=1:ndims(obj1);
            end
            osize=size(obj1);
            osize(dimension)=[];
            obj=Tcalculus(type,osize,dimension,obj1.TCindex,{},2);
        end

        function obj=compare(obj1,obj2,type)
        % used by (binary) min and max
            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);
            % scalar? multiply by ones of appropriate size
            if isempty(osize1) && ~isempty(osize2)
                obj1=growScalar(obj1,osize2);
                osize1=size(obj1);
            elseif isempty(osize2) && ~isempty(osize1)
                obj2=growScalar(obj2,osize1);
                osize2=size(obj2);
            end
            osize=size(obj1);
            obj=Tcalculus(type,osize,[],[obj1.TCindex;obj2.TCindex],{},2);
        end


        function obj=min(obj1,obj2,dimension)
        % min - minimum elements of a tensor
        %
        %   M=min(X), for a tensor X with d dimensions, returns a
        %   tensor M with d-1 dimensions whose entries are the
        %   smallest entries of X, taken along the 1st dimension.
        %
        %   M=min(X,Y), for two tensors X and Y with the same size,
        %   returns a tensor M also with the same size, but with
        %   entries taken from X or Y, depending on which entry is
        %   smallest. If X is a scalar, then M has the same size as Y
        %   and its entries are the smallest of the corresponding
        %   entry of Y or the (only) entry of X. Similarly if Y is a
        %   scalar.
        %
        %   M=min(X,[],vecdim) takes the minimum over all entries in
        %   the dimensions specified by vecdim, resulting in a tensor
        %   with the same size as X, but with the dimensions in vecdim
        %   removed. The simplest form min(X) performs similarly to
        %   min(X,[],1).
        %
        %   M=min(X,[],'all') returns a scalar with the smallest entry
        %   of X. Performs similarly to min(X,[],1:ndims(X))

            if nargin==1
                obj=reduce(obj1,1,'min');
            elseif nargin==3 && isempty(obj2)
                obj=reduce(obj1,dimension,'min');
            elseif nargin==2
                %error('min2 not implemented');
                obj=compare(obj1,obj2,'min2');
            end
        end
        function obj=max(obj1,obj2,dimension)
        % max - maximum elements of a tensor
        %
        %   M=max(X), for a tensor X with d dimensions, returns a
        %   tensor M with d-1 dimensions whose entries are the
        %   largest entries of X, taken along the 1st dimension.
        %
        %   M=max(X,Y), for two tensors X and Y with the same size,
        %   returns a tensor M also with the same size, but with
        %   entries taken from X or Y, depending on which entry is
        %   largest. If X is a scalar, then M has the same size as Y
        %   and its entries are the largest of the corresponding entry
        %   of Y or the (only) entry of X. Similarly if Y is a scalar.
        %
        %   M=max(X,[],vecdim) takes the maximum over all entries in
        %   the dimensions specified by vecdim, resulting in a tensor
        %   with the same size as X, but with the dimensions in vecdim
        %   removed. The simplest form max(X) performs similarly to
        %   max(X,[],1).
        %
        %   M=max(X,[],'all') returns a scalar with the largest entry
        %   of X. Performs similarly to max(X,[],1:ndims(X))

            if nargin==1
                obj=reduce(obj1,1,'max');
            elseif nargin==3 && isempty(obj2)
                obj=reduce(obj1,dimension,'max');
            elseif nargin==2
                %error('max2 not implemented');
                obj=compare(obj1,obj2,'max2');
            end
        end

        function obj=all(obj1,dimension)
        % all - Check if all tensor entries are nonzero
        %
        %   b=all(X,vecdim) checks if the entries of the tensor X are
        %   nonzero and performs the Boolean operation 'and' along the
        %   dimensions specified in vecdim, producing the logical
        %   value true if all entries are nonzero. The result is a
        %   tensor with the same size as X, but with the dimensions in
        %   vecdim removed.
        %
        %   b=all(X,'all') operates over all entries of the tensor X,
        %   resulting in a scalar that is equal to true if all entries
        %   are nonzero. Performs similarly to all(X,1:ndims(X))
            obj=reduce(obj1,dimension,'all');
        end
        function obj=any(obj1,dimension)
        % any - Check if any tensor entries are nonzero
        %
        %   b=any(X,vecdim) checks if the entries of the tensor X are
        %   nonzero and performs the Boolean operation 'or' along the
        %   dimensions specified in vecdim, producing the logical
        %   value true if at least one entry is nonzero. The result is
        %   a tensor with the same size as X, but with the dimensions
        %   in vecdim removed.
        %
        %   b=any(X,'all') operates over all entries of the tensor X,
        %   resulting in a scalar that is equal to true if at least
        %   one entry is nonzero. Performs similarly to
        %   any(X,1:ndims(X))
            obj=reduce(obj1,dimension,'any');
        end

        function obj=diag(obj1,k)
        % diag - Diagonal matrices and diagonals of a matrix.
        %
        %   diag(v,k) when v is an n-vector, returns a square matrix
        %   with n+abs(k) rows/columns, with the k-th diagonal equal
        %   to v. k=0 corresponds to the main diagonal, k>0 above the
        %   main diagonal and k<0 below.
        %
        %   diag(v) when v is an n-vector, returns the same as
        %   diag(v,0) and puts v in the main diagonal
        %
        %   diag(A) when A is an n-by-n a matrix, returns an n`-vector
        %   with the main diagonal of A.
        %
        % Attention: diag(A) with A a matrix doe *not* support a
        % second argument specifying a diagonal other than the main
        % diagonal.
            if nargin<2
                k=0;
            end
            osize1=size(obj1);
            if length(osize1)==1
                % vector->matrix
                %obj=tprod(obj1,1,Teye([osize1,osize1]),[1,2]); % only for k=0

                % in principle more efficient since it will never require
                % multiplications, but tprod is also pretty good at
                % handling multiplication by Teye;
                if k>=0
                    obj=vec2tensor(obj1,[osize1+k,osize1+k],[1:osize1;k+1:k+osize1]');
                else
                    obj=vec2tensor(obj1,[osize1-k,osize1-k],[1-k:osize1-k;1:osize1]');
                end
                updateFile2table(obj,1);
            elseif length(osize1)==2 && osize1(1)==osize1(2)
                if nargin>1
                    error('current implementation of diag(x) can only take 2nd argument when x is a vector not with size [%s]',...
                          index2str(osize1))
                end
                % matrix->vector
                obj=tprod(obj1,[1,1]);
                updateFile2table(obj,1);
            else
                obj1
                error('ambigous Tcalculus/diag()');
            end
        end

        function obj=trace(obj1)
        % trace - Sum of diagonal elements.
        %
        %  trace(A) is the sum of the diagonal elements of A, which is
        %  also the sum of the eigenvalues of A.
        %
        %  For more general tensors n-dimensional tensors, one must
        %  have size(A,1)==size(A,2), and trace returns a (n-2)-dimensional
        %  tensor obtained by summing over the same entries
        %  of the 1st two dimensions, as in
        %       Y(i,k,l)=sum_i A(i,i,j,k,l,...)
        %

            osize1=size(obj1);
            if length(osize1)>=2 && osize1(1)==osize1(2)
                obj=tprod(obj1,[-1,-1,1:length(osize1)-2]);
            else
                obj1
                error('ambigous Tcalculus/trace()');
            end
        end

        function obj=transpose(obj1)
        % .' or transpose - transpose of a real-values tensor
        %
        % Attention: TensCalc does not support complex-valued
        % variables and therefore transpose() and ctranspose() return
        % the same values.
            obj=ctranspose(obj1);
            updateFile2table(obj,1);
        end

        function obj=ctranspose(obj1)
        % ' or ctranspose - transpose of a real-values tensor
        %
        % Attention: TensCalc does not support complex-valued
        % variables and therefore transpose() and ctranspose() return
        % the same values.
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
            if length(osize1)<2
                error('Tcalculus: Transpose on scalar (0D-array) or vector (1D-array) not allowed nothing\n');
            end
            if length(osize1)>2
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

        function obj=permute(obj1,order)
            osize1=size(obj1);
            obj=tprod(obj1,order);
        end

        function obj=permute_matlab(obj1,order)
            osize1=size(obj1);
            obj=Tcalculus('permute_matlab',osize1(order),order,obj1.TCindex,{},1);
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
        % ldl - LDL factorization of a symmetric matrix
        %
        %   ldl(A) computes the LDL factorization of a symmetric
        %   matrix (n-by-n tensor).
        %
        %   ldl(A) actually uses "pychologically lower triangular
        %   matrices", i.e., matrices that are diagonal/triangular up
        %   to a permutation, with permutations selected to minimize
        %   the fill-in for sparse matrices and reduce computation
        %   time (see documentation for lu with sparse matrices).
        %
        %   The LDL factorization used requires all elements of the
        %   main diagonal to be nonzero, if this is not the case the
        %   lu() factorization should be used.
        %
        %   All entries of ``A`` above the main diagonal are ignored
        %   and assumed to be equal to the one below the main diagonal,
        %   *without performing any test regarding of whether or not
        %   this is true*.
        %
        % Attention: The output of this function includes the whole
        % factorization as a single entity, in a format that can be
        % passed to functions that require factorizations (such as
        % mldivide, inv, det, logdet, traceinv), but should *not* be
        % used by functions that are not expecting a factorized matrix
        % as an input. See ldl_l() and ldl_d() to get the factors.
            obj=factor(obj1,'ldl',varargin{:});
        end

        function obj=ldl_d(obj1)
        % ldl_d - Diagonal matrix of an LDL factorization
        %
        %   ldl_d(ldl(A)) return the main diagonal of the diagonal
        %   matrix in an LDL factorization
        %
        %   This function can only be applied to a matrix that has
        %   been factorized with ldl().
        %
        %   ldl(A) actually uses "pychologically lower triangular
        %   matrices", i.e., matrices that are diagonal/triangular up
        %   to a permutation, with permutations selected to minimize
        %   the fill-in for sparse matrices and reduce computation
        %   time (see documentation for lu with sparse matrices).
        %
        %   The LDL factorization used requires all elements of the
        %   main diagonal to be nonzero, if this is not the case the
        %   lu() factorization should be used.

        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.

            if ~strcmp(type(obj1),'ldl')
                error('ldl_d can only be called for LDL factorizations (not ''%s'')',type(obj1))
            end
            osize1=size(obj1,1);
            obj=Tcalculus('ldl_d',osize1,[],obj1.TCindex,{},1);
        end

        function obj=ldl_l(obj1)
        % ldl_l - Lower-trangular matrix of an LDL factorization
        %
        %   ldl_l(ldl(A)) return the lower-triangiular matrix in an
        %   LDL factorization
        %
        %   This function can only be applied to a matrix that has
        %   been factorized with ldl().
        %
        %   ldl(A) actually uses "pychologically lower triangular
        %   matrices", i.e., matrices that are diagonal/triangular up
        %   to a permutation, with permutations selected to minimize
        %   the fill-in for sparse matrices and reduce computation
        %   time (see documentation for lu with sparse matrices).
        %
        %   The LDL factorization used requires all elements of the
        %   main diagonal to be nonzero, if this is not the case the
        %   lu() factorization should be used.

        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.
            if ~strcmp(type(obj1),'ldl')
                error('ldl_l can only be called for LDL factorizations (not ''%s'')',type(obj1))
            end
            osize1=size(obj1);
            obj=Tcalculus('ldl_l',osize1,[],obj1.TCindex,{},1);
        end

        function obj=lu_l(obj1)
        % lu_l - Lower-trangular matrix of an LU factorization
        %
        %   lu_l(lu(A)) return the lower-triangular matrix in an LU
        %   factorization (which is guaranteed to be nonsingular)
        %
        %   This function can only be applied to a matrix that has
        %   been factorized with lu().
        %
        %   lu(A) actually uses "pychologically upper/lower triangular
        %   matrices", i.e., matrices that are diagonal/triangular up
        %   to a permutation, with permutations selected to minimize
        %   the fill-in for sparse matrices and reduce computation
        %   time (see documentation for lu with sparse matrices).
        %
        %   The LU factorization used requires all elements of the
        %   main diagonal to be nonzero, if this is not the case the
        %   lu() factorization should be used.

        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.

            if ~strcmp(type(obj1),'lu')
                error('lu_l can only be called for LU factorizations (not ''%s'')',type(obj1))
            end
            osize1=size(obj1);
            obj=Tcalculus('lu_l',osize1,[],obj1.TCindex,{},1);
        end

        function obj=lu_u(obj1)
        % lu_u - Upper-trangular matrix of an LU factorization
        %
        %   lu_u(lu(A)) return the upper-triangiular matrix in an LU
        %   factorization
        %
        %   This function can only be applied to a matrix that has
        %   been factorized with lu().
        %
        %   lu(A) actually uses "pychologically upper/lower triangular
        %   matrices", i.e., matrices that are diagonal/triangular up
        %   to a permutation, with permutations selected to minimize
        %   the fill-in for sparse matrices and reduce computation
        %   time (see documentation for lu with sparse matrices).
        %
        %   The LDL factorization used requires all elements of the
        %   main diagonal to be nonzero, if this is not the case the
        %   lu() factorization should be used.

        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.
            if ~strcmp(type(obj1),'lu')
                error('lu_u can only be called for LU factorizations (not ''%s'')',type(obj1))
            end
            osize1=size(obj1);
            obj=Tcalculus('lu_u',osize1,[],obj1.TCindex,{},1);
        end

        function obj=lu_d(obj1)
        % lu_d - Diagonal matrix of the U matrix in an LU factorization
        %
        %   lu_d(lu(A)) return the main diagonal of the
        %   upper-triangular matrix in an LU factorization. The main
        %   diagonal of the L matrix is equal to 1.
        %
        %   This function can only be applied to a matrix that has
        %   been factorized with lu().

            if ~strcmp(type(obj1),'lu')
                error('lu_d can only be called for LU factorizations (not ''%s'')',type(obj1))
            end
            osize1=size(obj1,1);
            obj=Tcalculus('lu_d',osize1,[],obj1.TCindex,{},1);
        end

        function obj=lu(obj1,varargin)
        % lu - LU factorization
        %
        %    lu(A) computes the LU factorization of a matrix (tensor
        %    with 2 dimensions).
        %
        %    lu(A) actually uses "pychologically upper/lower
        %    triangular matrices", i.e., matrices that are triangular
        %    up to a permutation, with permutations selected to
        %    minimize the fill-in for sparse matrices and reduce
        %    computation time (see documentation for lu with sparse
        %    matrices).
        %
        % Attention: The output of this function includes the whole
        % factorization as a single entity, in a format that can be
        % passed to functions that require factorizations (such as
        % mldivide, inv, det, logdet, traceinv), but should *not* be
        % used by functions that are not expecting a factorized matrix
        % as an input.  See lu_l(), lu_u(), and lu_d() to get the factors.
            obj=factor(obj1,'lu',varargin{:});
        end
        function obj=lu_sym(obj1,varargin)
        % lu - LU factorization of a symmetric matrix
        %
        %    lu(A) computes the LU factorization of a symmetric matrix
        %    (tensor with 2 dimensions).
        %
        %    lu(A) actually uses "pychologically lower triangular
        %    matrices", i.e., matrices that are triangular up to a
        %    permutation, with permutations selected to minimize the
        %    fill-in for sparse matrices and reduce computation time
        %    (see documentation for lu with sparse matrices).
        %
        %    All entries of ``A`` above the main diagonal are ignored
        %    and assumed to be equal to the one below the main diagonal,
        %    *without performing any test regarding of whether or not
        %    this is true*.
        %
        % Attention: The output of this function includes the whole
        % factorization as a single entity, in a format that can be
        % passed to functions that require factorizations (such as
        % mldivide, inv, det, logdet, traceinv), but should *not* be
        % used by functions that are not expecting a factorized matrix
        % as an input.

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

        % function obj=det(obj1)
        %     %obj1=toCalculus(obj1);
        %     osize1=size(obj1);
        %     if length(osize1)~=2
        %         error('det is only defined for 2D matrices');
        %     end
        %     if osize1(1)~=osize1(2)
        %         error('det is only defined for square matrices');
        %     end
        %     if osize1(1)==1
        %         obj=reshape(obj1,[]);
        %     elseif osize1(1)==2
        %         a11=subsref(obj1,struct('type','()','subs',{{1,1}}));
        %         a12=subsref(obj1,struct('type','()','subs',{{1,2}}));
        %         a21=subsref(obj1,struct('type','()','subs',{{2,1}}));
        %         a22=subsref(obj1,struct('type','()','subs',{{2,2}}));
        %         obj=reshape(a11*a22-a12*a21,[]);
        %         updateFile2table(obj,1);
        %     elseif osize1(1)==3
        %         a11=subsref(obj1,struct('type','()','subs',{{1,1}}));
        %         a12=subsref(obj1,struct('type','()','subs',{{1,2}}));
        %         a13=subsref(obj1,struct('type','()','subs',{{1,3}}));
        %         a21=subsref(obj1,struct('type','()','subs',{{2,1}}));
        %         a22=subsref(obj1,struct('type','()','subs',{{2,2}}));
        %         a23=subsref(obj1,struct('type','()','subs',{{2,3}}));
        %         a31=subsref(obj1,struct('type','()','subs',{{3,1}}));
        %         a32=subsref(obj1,struct('type','()','subs',{{3,2}}));
        %         a33=subsref(obj1,struct('type','()','subs',{{3,3}}));
        %         obj=reshape(a11*a22*a33-a11*a23*a32-a12*a21*a33...
        %                     +a12*a23*a31+a13*a21*a32-a13*a22*a31,[]);
        %         updateFile2table(obj,1);
        %     elseif osize1(1)==4
        %         a11=subsref(obj1,struct('type','()','subs',{{1,1}}));
        %         a12=subsref(obj1,struct('type','()','subs',{{1,2}}));
        %         a13=subsref(obj1,struct('type','()','subs',{{1,3}}));
        %         a14=subsref(obj1,struct('type','()','subs',{{1,4}}));
        %         a21=subsref(obj1,struct('type','()','subs',{{2,1}}));
        %         a22=subsref(obj1,struct('type','()','subs',{{2,2}}));
        %         a23=subsref(obj1,struct('type','()','subs',{{2,3}}));
        %         a24=subsref(obj1,struct('type','()','subs',{{2,4}}));
        %         a31=subsref(obj1,struct('type','()','subs',{{3,1}}));
        %         a32=subsref(obj1,struct('type','()','subs',{{3,2}}));
        %         a33=subsref(obj1,struct('type','()','subs',{{3,3}}));
        %         a34=subsref(obj1,struct('type','()','subs',{{3,4}}));
        %         a41=subsref(obj1,struct('type','()','subs',{{4,1}}));
        %         a42=subsref(obj1,struct('type','()','subs',{{4,2}}));
        %         a43=subsref(obj1,struct('type','()','subs',{{4,3}}));
        %         a44=subsref(obj1,struct('type','()','subs',{{4,4}}));
        %         obj=reshape(a14*a23*a32*a41-a13*a24*a32*a41-a14*a22*a33*a41...
        %                     +a12*a24*a33*a41+a13*a22*a34*a41-a12*a23*a34*a41...
        %                     -a14*a23*a31*a42+a13*a24*a31*a42+a14*a21*a33*a42...
        %                     -a11*a24*a33*a42-a13*a21*a34*a42+a11*a23*a34*a42...
        %                     +a14*a22*a31*a43-a12*a24*a31*a43-a14*a21*a32*a43...
        %                     +a11*a24*a32*a43+a12*a21*a34*a43-a11*a22*a34*a43...
        %                     -a13*a22*a31*a44+a12*a23*a31*a44+a13*a21*a32*a44...
        %                     -a11*a23*a32*a44-a12*a21*a33*a44+a11*a22*a33*a44,[]);
        %         updateFile2table(obj,1);
        %     else
        %         obj=Tcalculus('det',[],[],obj1.TCindex,{},1);
        %     end
        % end

        function obj=det(obj1)
        % det - Determinant of a factorized matrix
        %
        %   det(lu(A)) or det(ldl(A)) return the determinant of the
        %   square matrix A.
        %
        %   This function can only be applied to a matrix that has
        %   been factorized by lu() or ldl().
        %
        %   The ldl() factorization should only be used for symmetric
        %   matrices with nonzero elements in the main diagonal.
        %
        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.
            if ~strcmp(type(obj1),'ldl') && ~strcmp(type(obj1),'lu') && ~strcmp(type(obj1),'lu_sym')
                error('det can only be called for LDL/LU factorizations (not ''%s'')',type(obj1))
            end
            obj1=toCalculus(obj1);
            osize1=size(obj1);
            if length(osize1)~=2
                error('logdet is only defined for 2D matrices');
            end
            if osize1(1)~=osize1(2)
                error('logdet is only defined for square matrices');
            end
            obj=Tcalculus('det',[],[],obj1.TCindex,{},1);
        end

        function obj=logdet(obj1)
        % logdet - Natural logarithm of the determinant of a factorized matrix
        %          (with positive determinant)
        %
        %   logdet(lu(A)) or logdet(ldl(A)) return the natural
        %   logarithm of the determinant of the square matrix A.
        %
        %   This function results in the same value as log(det(A)),
        %   but in the context of optimizations, it is more efficient
        %   to use ``logdet(A)`` because the latter simplifies the
        %   computation derivatives.
        %
        %   logdet() can only be applied to a matrix that has been
        %   factorized by lu() or ldl().
        %
        %   The ldl() factorization should only be used for symmetric
        %   matrices with nonzero elements in the main diagonal.
        %
        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.
        %
        %   |tenscalc| assumes that det(A)>0 and actually returns
        %   log(abs(det(A))), but ignores the abs() operation when
        %   computing derivatives of logdet.

            if ~strcmp(type(obj1),'ldl') && ~strcmp(type(obj1),'lu') && ~strcmp(type(obj1),'lu_sym')
                error('logdet can only be called for LDL/LU factorizations (not ''%s'')',type(obj1))
            end
            obj1=toCalculus(obj1);
            osize1=size(obj1);
            if length(osize1)~=2
                error('logdet is only defined for 2D matrices');
            end
            if osize1(1)~=osize1(2)
                error('logdet is only defined for square matrices');
            end
            obj=Tcalculus('logdet',[],[],obj1.TCindex,{},1);
        end

        function obj=traceinv(obj1)
        % traceinv - Trace of the inverse of a factorized matrix
        %
        %   traceinv(lu(A)) or traceinv(ldl(A)) return the natural
        %   logarithm of the determinant of the square matrix A.
        %
        %   This function results in the same value as trace(inv(A)),
        %   but in the context of optimizations, it is more efficient
        %   to use ``traceinv(A)`` because the latter simplifies the
        %   computation of derivatives.
        %
        %   traceinv() can only be applied to a matrix that has been
        %   factorized by lu() or ldl().
        %
        %   The ldl() factorization should only be used for symmetric
        %   matrices with nonzero elements in the main diagonal.
        %
        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.
            if ~strcmp(type(obj1),'ldl') && ~strcmp(type(obj1),'lu') && ~strcmp(type(obj1),'lu_sym')
                error('traceinv can only be called for LDL/LU factorizations (not ''%s'')',type(obj1))
            end
            obj1=toCalculus(obj1);
            osize1=size(obj1);
            if length(osize1)~=2
                error('traceinv is only defined for 2D matrices');
            end
            if osize1(1)~=osize1(2)
                error('traceinv is only defined for square matrices');
            end
            obj=Tcalculus('traceinv',[],[],obj1.TCindex,{},1);
        end

        function obj=inv(obj1)
        % inv - Inverse of a factorized matrix
        %
        %   inv(lu(A)) or inv(ldl(A)) return the inverse of the square matrix A.
        %
        %   This function can only be applied to a matrix that has
        %   been factorized by lu() or ldl().
        %
        %   The ldl() factorization should only be used for symmetric
        %   matrices with nonzero elements in the main diagonal.
        %
        %   When an the LDL factorization is used, all entries of A
        %   above the main diagonal are ignored and assumed to be
        %   equal to the one below the main diagonal, without
        %   performing any test regarding of whether or not this is
        %   true.
            if ~strcmp(type(obj1),'ldl') && ~strcmp(type(obj1),'lu') && ~strcmp(type(obj1),'lu_sym')
                error('inv can only be called for LDL/LU factorizations (not ''%s'')',type(obj1))
            end
            obj1=toCalculus(obj1);
            osize1=size(obj1);
            if length(osize1)~=2
                error('inv is only defined for 2D matrices');
            end
            if osize1(1)~=osize1(2)
                error('inv is only defined for square matrices');
            end
            obj=Tcalculus('inv',osize1,[],obj1.TCindex,{},1);
        end

        function obj=componentwise(obj1,fun,cfun,varargin)
        % y=componentwise(x,fun,cfun,derivative1,derivative2)
        %
        % Returns a tensor with the same size of ''x'' that results
        % from applying the function ''fun'' to each entry of ''x'':
        %       y_{ijk} = fun(x_{ijk})
        % ''fun'' should be a function handle.
        %
        % The string ''cfun'' provides code that can be used to
        % evaluate the function when '%s' is replaced the entry:
        %
        % The optional argument ''derivative'' is a handle to a matlab
        % function that returns a Tcalculus tensor with the 1st
        % derivative. Typically, this handle includes a call to
        % componentwise, as in the example
        %   componentwise(x,@sin,'sin(%s)',...
        %                         @(x)componentwise(x,@cos,'cos(%s)',...
        %                                   @(x)-componentwise(x,@sin,'sin(%s)')));
        % but it may not, as in the example
        %   componentwise(x,@(x)2*x,'2*%s',@(x)2*Tones(size(x)),@(x)Tzeros(size(x)));

            if ~isa(fun,'function_handle')
                error('componentwise 2nd argument must be a function handle, class ''%s'' instead\n',class(fun))
            end

            if ~ischar(cfun)
                error('componentwise 3rd argument must be a character string, class ''%s'' instead\n',class(cfun))
            end

            for i=1:length(varargin)
                if ~isa(varargin{i},'function_handle')
                    error('componentwise %dth argument must be a function handle, class ''%s'' instead\n',3+i,class(varargin{i}));
                end
            end
            obj=Tcalculus('componentwise',size(obj1),{fun,cfun,varargin{:}},obj1.TCindex,{},1);
        end

        function obj=exp(obj1)
        % Exponential of tensor entries.
            obj=componentwise(obj1,@exp,'exp(%s)',@exp);
            updateFile2table(obj,1);
        end
        function obj=log(obj1)
        % Natural logarithm of tensor entries.
            obj=componentwise(obj1,@log,'log(%s)',@(x)1./x);
            updateFile2table(obj,1);
        end
        function obj=sin(obj1)
        % Sine of tensor entries in radians.
            obj=componentwise(obj1,@sin,'sin(%s)',@cos);
            updateFile2table(obj,1);
        end
        function obj=cos(obj1)
        % Cosine of tensor entries in radians.
            obj=componentwise(obj1,@cos,'cos(%s)',@(x)-sin(x));
            updateFile2table(obj,1);
        end
        function obj=tan(obj1)
        % Tangent of tensor entries in radians.
            obj=componentwise(obj1,@tan,'tan(%s)',@sec2);
            updateFile2table(obj,1);
        end
        function obj=sec2(obj1)
            obj=componentwise(obj1,@(x)sec(x).^2,'1/(cos(%s)*cos(%s))',@(x)2*sec2(x).*tan(x));
            updateFile2table(obj,1);
        end
        function obj=atan(obj1)
        % Inverse tangent in radians of tensor entries.
            obj=componentwise(obj1,@atan,'atan(%s)',@datan);
            updateFile2table(obj,1);
        end
        function obj=datan(obj1)
            obj=componentwise(obj1,@(x)1./(1+x.^2),'1/(1+%s*%s)',@ddatan);
            updateFile2table(obj,1);
        end
        function obj=ddatan(obj1)
            obj=componentwise(obj1,@(x)-2*x./(1+x.^2).^2,'-2*%s/((1+%s*%s)*(1+%s*%s))');
            updateFile2table(obj,1);
        end
        function obj=sqr(obj1)
        % sqr - Square of tensor entries
        %
        %    sqr(X) returns the element-wise square of X, i.e., X.^2 or X.*X
            obj=componentwise(obj1,@(x)x.^2,'%s*%s',@(x)2*x);
            updateFile2table(obj,1);
        end
        function obj=cube(obj1)
        % cube - Cube of tensor entries
        %
        %    cube(X) returns the element-wise cube of X, i.e., X.^3 or X.*X.*X
            obj=componentwise(obj1,@(x)x.^3,'%s*%s*%s',@(x)3*x.*x);
            updateFile2table(obj,1);
        end
        function obj=power(obj1,p)
        % power - Power of tensor entries
        %
        %    X.^Y or power(x,y) returns the element-wise X raised to
        %    the power Y.
        %
        %    The power Y must be a regular numeric scalar, *not* a
        %    Tcalculus symbolic expression
        %
        % Attention: unlike in matlab's regular power(), the power Y must
        % be a scalar.
        %
            if length(p)>1
                error('Tcalculus/power only supports scalar exponents');
            end
            if ~isnumeric(p)
                error('Tcalculus/power only supports constant exponents');
            end
            if p==1
                obj=obj1;
            else
                obj=componentwise(obj1,eval(sprintf('@(x)x.^%.16g',p)),...
                                  sprintf('pow(%%s,%.16g)',p),...
                                  eval(sprintf('@(x)%.16g*power(x,%.16g)',p,p-1)));
                updateFile2table(obj,1);
            end
        end
        function obj=sqrt(obj1)
        % cube - Square root of tensor entries
            obj=componentwise(obj1,@sqrt,'sqrt(%s)',@(x).5./sqrt(x)); % reuses sqrt computation
            %obj=componentwise(obj1,@sqrt,'sqrt(%s)',@(x).5*power(x,-.5));
            %obj=componentwise(obj1,@sqrt,'sqrt(%s)',@dsqrt);
            updateFile2table(obj,1);
        end
        function obj=dsqrt(obj1)
            obj=componentwise(obj1,@(x).5./sqrt(x),'.5*pow(%s,-.5)',@ddsqrt);
            updateFile2table(obj,1);
        end
        function obj=ddsqrt(obj1)
            obj=componentwise(obj1,@(x)-.25./x.^1.5,'-.25*pow(%s,-1.5)');
            updateFile2table(obj,1);
        end
        function obj=round(obj1)
        % round - Rounds tensor entries to the nearest integer
        %
        %    round(X) returns as tensor with the same size as X, with
        %    every entry of X rounded to the nearest integer.
        %
        % Attention: unlike in regular matlab, this function does
        % *not* support a second argument specifying a desired number
        % of digits for rounding to a decimal.
            obj=componentwise(obj1,@round,'round(%s)',@(x)Tzeros(size(x)));
            updateFile2table(obj,1);
        end
        function obj=ceil(obj1)
        % ceil - Rounds tensor entries to the nearest integer towards +inf
        %
        %    ceil(X) returns as tensor with the same size as X, with
        %    every entry of X rounded to the nearest integer towards
        %    +inf.
            obj=componentwise(obj1,@ceil,'ceil(%s)',@(x)Tzeros(size(x)));
            updateFile2table(obj,1);
        end
        function obj=floor(obj1)
        % floor - Rounds tensor entries to the nearest integer towards -inf
        %
        %    floor(X) returns as tensor with the same size as X, with
        %    every entry of X rounded to the nearest integer towards
        %    -inf.
            obj=componentwise(obj1,@floor,'floor(%s)',@(x)Tzeros(size(x)));
            updateFile2table(obj,1);
        end
        function obj=abs(obj1)
        % abs - Absolute value of tensor entries
            obj=componentwise(obj1,@abs,'fabs(%s)',@sign);
            updateFile2table(obj,1);
        end
        function obj=sign(obj1)
        % sign - Signum function.
        %
        %    sign(X) returns as tensor with the same size as X, with
        %    an entry equal to 1, 0, or -1, depending on whether the
        %    corresponding entry of X is positive, zero, or negative,
        %    respectively.
            obj=componentwise(obj1,@sign,'(%s>0)?1:((%s<0)?-1:0)',@(x)Tzeros(size(x)));
            updateFile2table(obj,1);
        end
        function obj=heaviside(obj1)
        % heaviside - Step or heaviside function applied to the tensor entries.
        %
        %    heaviside(X) returns as tensor with the same size as X, with
        %    an entry equal to 1, .5, or 0, depending on whether the
        %    corresponding entry of X is positive, zero, or negative,
        %    respectively.
            obj=componentwise(obj1,@heaviside,'(%s>0)?1:((%s<0)?0:.5)',@(x)Tzeros(size(x)));
            updateFile2table(obj,1);
        end
        function obj=relu(obj1)
        % relu - Rectified linear unit activation function.
        %
        %    relu(X) returns a tensor with the same size as X, with
        %    each entry equal the the corresponding entry of X or 0,
        %    depending on whether the entry is positive or not. Same
        %    as max(X,0).
            obj=componentwise(obj1,@(x)max(x,0),'(%s>0)?%s:0',@heaviside);
            updateFile2table(obj,1);
        end
        function obj=srelu(obj1)
        % srelu - Soft rectified linear unit activation function.
        %
        %    srelu(X) returns a tensor with the same size as X, with
        %    the "soft" rectified linear unit activation function of
        %    the entries of X. Same as log(1+exp(x)) or
        %    x+log(1+exp(-x))
            obj=componentwise(obj1,@(x)log(1+exp(x)),'log(1+exp(%s))',@sheaviside);
            updateFile2table(obj,1);
        end
        function obj=sheaviside(obj1)
        % sheaviside - Soft heaviside function.
        %
        %    sheaviside(X) returns a tensor with the same size as X,
        %    with the "soft" heaviside function of the entries of
        %    X. Same as 1./(1+exp(-x))
            obj=componentwise(obj1,@(x)1./(1+exp(-x)),'1/(1+exp(-%s))',@dsheaviside);
            updateFile2table(obj,1);
        end
        function obj=dsheaviside(obj1)
        % dsheaviside - Derivative of the soft heaviside function.
        %
        %    dsheaviside(X) returns a tensor with the same size as X,
        %    with the derivative of the "soft" heaviside function of
        %    the entries of X. Same as 1./(2+exp(x)+exp(-x))
            obj=componentwise(obj1,@(x)1./(2+exp(x)+exp(-x)),'1/(2+exp(%s)+exp(-%s))');
            updateFile2table(obj,1);
        end
        function obj=normpdf(obj1)
        % normpdf - Normal probability density function (pdf).
        %
        %    normpdf(X) returns the pdf of the standard normal
        %    distribution evaluated at the values in X.
        %
        % Attention: unlike in ragular matlab, this function does
        % *not* support second and third arguments specifying a mean
        % and standard deviation different than 0 and 1,
        % respectively.
            obj=componentwise(obj1,@normpdf,'exp(-%s*%s/2)/sqrt(2*M_PI)',@(x)-x.*normpdf(x));
            updateFile2table(obj,1);
        end
        function obj=lngamma(obj1)
            obj=componentwise(obj1,@(x)log(gamma(x)),'lngamma(%s)',@digamma);
            updateFile2table(obj,1);
        end

        function obj=compose(obj1,varargin)
        % y=compose(x,fun,fun1,fun2)
        %
        % Returns a tensor that results from applying the function
        % ''fun'' to each entry of ''obj'', e.g.,:
        %       y_{ijk} = fun(x_{ijk})
        % ''fun'' should be a function handle. When the result of ''fun''
        % is itself a tensor, the dimensions of fun are appended to the
        % end of the dimensions of x, e.g.,:
        %       y_{ijklm} = fun(x_{ijk})_{lmn}
        %
        % The optional argument ''fun1'' is a handle to a matlab
        % function that returns a Tcalculus tensors with the 1st
        % derivatives.
        %
        % The optional argument ''fun2'' is a handle to a matlab
        % function that returns a Tcalculus tensors with the 2nd
        % derivatives.
        %
        % Attention: compose can handle nonscalar functions, which
        % componentwise cannot.

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

        % function obj=sqr(obj1)
        %     obj=compose(obj1,@(x__)x__.^2,@(x__)2*x__,@(x__)2*ones(size(x__)));
        %     updateFile2table(obj,1);
        % end
        % function obj=sqrt(obj1)
        %     obj=compose(obj1,@(x__)sqrt(x__),@(x__).5./sqrt(x__),@(x__)-.25./x__.^1.5);
        %     updateFile2table(obj,1);
        % end
        % function obj=cos(obj1)
        %     obj=compose(obj1,@(x__)cos(x__),@(x__)-sin(x__),@(x__)-cos(x__));
        %     updateFile2table(obj,1);
        % end
        % function obj=sin(obj1)
        %     obj=compose(obj1,@(x__)sin(x__),@(x__)cos(x__),@(x__)-sin(x__));
        %     updateFile2table(obj,1);
        % end
        % function obj=tan(obj1)
        %     obj=compose(obj1,@(x__)tan(x__),@(x__)sec(x__).^2,@(x__)2*sec(x__).^2.*tan(x__));
        %     updateFile2table(obj,1);
        % end
        % function obj=atan(obj1)
        %     obj=compose(obj1,@(x__)atan(x__),@(x__)1./(1+x__.^2),@(x__)-2*x__./(1+x__.^2).^2);
        %     updateFile2table(obj,1);
        % end
        % function obj=exp(obj1)
        %     obj=compose(obj1,@(x__)exp(x__),@(x__)exp(x__),@(x__)exp(x__));
        %     updateFile2table(obj,1);
        % end
        % function obj=log(obj1)
        %     obj=compose(obj1,@(x__)log(x__),@(x__)1./x__,@(x__)-1./x__.^2);
        %     updateFile2table(obj,1);
        % end
        % function obj=reciprocal(obj1)
        %     obj=compose(obj1,@(x__)1./x__,@(x__)-1./x__.^2,@(x__)2./x__.^3);
        %     updateFile2table(obj,1);
        % end
        % function obj=cube(obj1)
        %     obj=compose(obj1,@(x__)x__.^3,@(x__)3*x__.^2,@(x__)6*x__);
        %     updateFile2table(obj,1);
        % end
        % function obj=normpdf(obj1)
        %     obj=compose(obj1,@(x__)exp(-x__.^2/2)./sqrt(2*pi),...
        %                       @(x__)-x__.*exp(-x__.^2/2)./sqrt(2*pi),...
        %                       @(x__)(x__.^2-1).*exp(-x__.^2/2)./sqrt(2*pi));
        %     updateFile2table(obj,1);
        % end
        % function obj=round(obj1)
        %     obj=compose(obj1,@(x__)round(x__));
        %     updateFile2table(obj,1);
        % end

        % function obj=ceil(obj1)
        %     obj=compose(obj1,@(x__)ceil(x__));
        %     updateFile2table(obj,1);
        % end

        % function obj=floor(obj1)
        %     obj=compose(obj1,@(x__)floor(x__));
        %     updateFile2table(obj,1);
        % end
        % function obj=abs(obj1)
        %     obj=compose(obj1,@(x__)abs(x__),@(x__)sign(x__),@(x__)zeros(size(x__)));
        %     updateFile2table(obj,1);
        % end
        % function obj=relu(obj1)
        %     obj=compose(obj1,@(x__)relu(x__),@(x__)heaviside(x__),@(x__)zeros(size(x__)));
        %     updateFile2table(obj,1);
        % end
        % function obj=heaviside(obj1)
        %     obj=compose(obj1,@(x__)heaviside(x__),@(x__)zeros(size(x__)),@(x__)zeros(size(x__)));
        %     updateFile2table(obj,1);
        % end
        % function obj=srelu(obj1)
        %     obj=compose(obj1,@(x__)log(1+exp(x__)),...
        %                 @(x__)1./(1+exp(-x__)),@(x__)1./(2+exp(-x__)+exp(x__)));
        %     updateFile2table(obj,1);
        % end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                 Binary operators                         %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=plus(obj1,obj2,ind1,ind2)
        % + or plus
        %
        %   X+Y or plus(X,Y) add entry-by-entry tensors of the same
        %   size or a scalar (tensor with 0 dimensions) with a tensor
        %   of arbitrary size.
        %
        % Attention: unlike matlab's regular plus(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.

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
        % - or minus - Minus
        %
        %   X-Y or minus(X,Y) subtract entry-by-entry tensors of the
        %   same size or a scalar (tensor with 0 dimensions) with a
        %   tensor of arbitrary size.
        %
        % Attention: unlike matlab's regular minus(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.

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
        % == or eq -  Equal
        %
        %   X==Y performs an entry-by-entry equality comparison of
        %   tensors of the same size or a scalar (tensor with 0
        %   dimensions) with a tensor of arbitrary size.
        %
        % Attention: unlike matlab's regular eq(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.
            obj=relational(obj1,obj2,'iszero');
        end

        function obj=ge(obj1,obj2)
        % >= or ge - Greater than or equal to
        %
        %   X<=Y performs an entry-by-entry greater than or equal to
        %   comparison of tensors of the same size or a scalar (tensor
        %   with 0 dimensions) with a tensor of arbitrary size.
        %
        %   From the perspective of a constrained optimization
        %   numerical solver, due to finite numerical precision, X>=Y
        %   and X>Y represent the same constraint.
        %
        % Attention: unlike matlab's regular ge(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.
            obj=relational(obj1,obj2,'ispositive');
        end

        function obj=gt(obj1,obj2)
        % > or gt - Greater than
        %
        %   X<=Y performs an entry-by-entry greater than comparison of
        %   tensors of the same size or a scalar (tensor with 0
        %   dimensions) with a tensor of arbitrary size.
        %
        %   From the perspective of a constrained optimization
        %   numerical solver, due to finite numerical precision, X>=Y
        %   and X>Y represent the same constraint.
        %
        % Attention: unlike matlab's regular gt(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.
            obj=relational(obj1,obj2,'ispositive');
        end

        function obj=le(obj2,obj1)
        % <= or le - Smaller than or equal to
        %
        %   X<=Y performs an entry-by-entry smaller than or equal to
        %   comparison of tensors of the same size or a scalar (tensor
        %   with 0 dimensions) with a tensor of arbitrary size.
        %
        %   From the perspective of a constrained optimization
        %   numerical solver, due to finite numerical precision, X<=Y
        %   and X<Y represent the same constraint.
        %
        % Attention: unlike matlab's regular le(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.
            obj=relational(obj1,obj2,'ispositive');
        end

        function obj=lt(obj2,obj1)
        % < or lt - Smaller than
        %
        %   X<=Y performs an entry-by-entry smaller than comparison of
        %   tensors of the same size or a scalar (tensor with 0
        %   dimensions) with a tensor of arbitrary size.
        %
        %   From the perspective of a constrained optimization
        %   numerical solver, due to finite numerical precision, X<=Y
        %   and X<Y represent the same constraint.
        %
        % Attention: unlike matlab's regular lt(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.
            obj=relational(obj1,obj2,'ispositive');
        end

        function obj=mge(obj1,obj2)
            obj=relational(obj1,obj2,'ispositivedefinite');
        end

        function obj=mle(obj1,obj2)
            obj=relational(obj2,obj1,'ispositivedefinite');
        end

        function obj=times(obj1,obj2,times2tprod)
        % .* or times - Element-wise tensor multiply
        %
        %   X.*Y or times(X,Y) multiply entry-by-entry tensors of the
        %   same size or a scalar (tensor with 0 dimensions) with a
        %   tensor of arbitrary size.
        %
        % Attention: unlike matlab's regular times(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.

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

        function obj=mtimes(obj1,obj2,mtimes2tprod)
        % * or mtimes - Matrix multiply
        %
        %   X*Y is the matrix product of X and Y, which adapts to the
        %   size of the operands as follows:
        %   * regular matrix multiplication when X and Y are both
        %     matrices (tensors with 2 dimensions) with
        %     size(X,2)==size(Y,1)
        %   * matrix by column-vector multiplication when X is a matrix
        %     (tensor with 2 dimensions) and Y a vector (tensor with 1
        %     dimension) with size(X,2)==size(Y,1)
        %   * row vector by matrix multiplication when X is a vector
        %     (tensor with 1 dimension) and Y a matrix (tensors with 2
        %     dimensions) with size(X,1)==size(Y,1)
        %   * inner product when X and Y are both vectors (tensors
        %     with 1 dimension) with size(X,1)==size(Y,1)
        %   * entry-by-entry multiplication with either X or Y are
        %     scalars (tensor with 0 dimensions).
        %
        % Attention: Depending on the sizes of the parameters, this
        % operation may behave quite differently from MATLAB©’s matrix
        % multiplication mtimes().

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
                    obj=repmat(obj1,[1,osize2(2)]);
                elseif strcmp(type(obj1),'ones') && osize1(2)==1 && osize1(1)>1
                    obj=repmat(obj2,[osize1]); % needed for ones(100,1)*rand(1,36)
                elseif (strcmp(type(obj1),'ones') || strcmp(type(obj2),'ones')) &&...
                        ~isempty(osize1) && ~isempty(osize2)
                    obj=osize1(end)*Tones([osize1(1),osize2(end)]);
                else
                    obj=Tcalculus('mtimes',osize,[],[obj1.TCindex;obj2.TCindex],{},1);
                end
            end
            updateFile2table(obj,1);
        end

        function obj=rdivide(obj1,obj2)
        % ./ or rdivide - Right tensor divide
        %
        %   X./Y or rdividde(X,Y) right-divide entry-by-entry tensors
        %   of the same size or a scalar (tensor with 0 dimensions)
        %   with a tensor of arbitrary size.
        %
        % Attention: unlike matlab's regular times(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.

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

            if  ~myisequal(osize1,osize2) && ~isempty(osize1) && ~isempty(osize2)
                obj1,obj2
                error('can only element-wise divide objects of the same size, divide by scalar, or divide scalar: [%s] ~= [%s]\n',...
                      index2str(osize1),index2str(osize2))
            end

            if strcmp(type(obj1),'zeros')
                obj=obj1;
            else
                if isempty(osize1)
                    obj=Tcalculus('rdivide',osize2,[],[obj1.TCindex;obj2.TCindex],{},1);
                elseif isempty(osize2)
                    obj=Tcalculus('rdivide',osize1,[],[obj1.TCindex;obj2.TCindex],{},1);
                else
                    obj=Tcalculus('rdivide',osize1,[],[obj1.TCindex;obj2.TCindex],{},1);
                end
            end
        end

        function obj=ldivide(obj1,obj2)
        % .\ or ldivide - Left tensor divide
        %
        %   X.\Y or ldividde(X,Y) left-divide entry-by-entry tensors
        %   of the same size or a scalar (tensor with 0 dimensions)
        %   with a tensor of arbitrary size.
        %
        % Attention: unlike matlab's regular times(), expansion upon
        % singleton dimensions is not performed automatically to match
        % the tensors' sizes.

            obj=rdivide(obj2,obj1);
            updateFile2table(obj,1);
        end

        function obj=mrdivide(obj1,obj2)
        % Code generation not (yet) implemented
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
        % \ or mldivide - Left matrix divide.
        %
        %   A\B is the matrix division of A into B, which is rougly
        %   the same as inv(A)*B
        %
        %   A must be a square matrix, but B can be a more genarl tensor,
        %   as long as its 1st dimension matches the size of A.
            [obj1,obj2]=toCalculus(obj1,obj2);
            osize1=size(obj1);
            osize2=size(obj2);

            if isempty(osize2)
                obj=ldivide(obj1,obj2);
                updateFile2table(obj,1);
                return
            end

            if length(osize1) ~= 0 && (length(osize1) ~= 2 || osize1(1)~=osize1(2) ) %|| length(osize2) > 2
                obj1,obj2
                error('1st argument of mldivide must be a square matrix (not [%s]))\n',...
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

        % function obj=tprod(varargin) - in separate file

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
        % cat - Concatenate tensors.
        %
        %   Y=cat(dim,A,B,C,...) concatenates the tensors
        %   A,B,C,... along the dimension dim. All the input tensors
        %   A,B,C,... should have similar sizes, except for the
        %   dimension dim along which concatenation will take
        %   place. The resulting tensor Y will also have a similar
        %   size, except for the dimension dim, for which we will have
        %        size(Y,dim)=size(A,dim)+size(B,dim)+...

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


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   Derivatives                            %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % function grad=gradient(obj,var) - in separate file

        function [hess,grad]=hessian(obj,var1,var2)
        % hessian - Hessian of a tensor-values symbolic expression with
        %          respect to a pair of tensor-valued symbolic variables
        %
        %    [h,g]=hessian(f,x), returns a tensor g with the partial
        %    derivatives of the entries of f with respect to the
        %    entries of the variable x, and a tensor h with the
        %    partial second derivatives of the entries of f with
        %    respect to the entries of the variable x
        %
        %      When f is a tensor with size
        %        [n1,n2,...,nN]
        %      and x a tensor-valued variable (created with Tvariable) with size
        %        [m1,m2,...,mM]
        %      then
        %        h=hessian(f,x)
        %      results in a tensor h with size
        %        [n1,n2,...,nN,m1,m2,...,mM,m1,m2,...,mM]
        %      with
        %        h(i1,i2,...,iN,j1,j2,...,jM,k1,k2,....,kM)
        %             =d2 f(i1,i2,...,iN) / d x(j1,j2,...,jM)d x(k1,k2,...,kM)
        %
        %    [h,g]=hessian(f,x,y), returns a tensor g with the partial
        %    derivatives of the entries of f with respect to the
        %    entries of the variable x, and a tensor with the partial
        %    second derivatives of the entries of f with respect to
        %    the entries of the variables x and y
        %
        %      When f is a tensor with size
        %        [n1,n2,...,nN]
        %      and x a tensor-valued variable (created with Tvariable) with size
        %        [m1,m2,...,mM]
        %      and x a tensor-valued variable (created with Tvariable) with size
        %        [l1,l2,...,lL]
        %      then
        %        h=hessian(f,x,y)
        %      results in a tensor y with size
        %        [n1,n2,...,nN,m1,m2,...,mM,l1,l2,...,lL]
        %      with
        %        y(i1,i2,...,iN,j1,j2,...,jM,k1,k2,....,kL)
        %             =d2 f(i1,i2,...,iN) / d x(j1,j2,...,jM) d y(k1,k2,...,kL)
        %
        % Copyright 2012-2017 Joao Hespanha

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
        % When 'obj' is a structure, the substitution is applied to every elements
        % of the structure
        %
        % When 'var' is a cell array, obj1 should be a 1-dimensional vector that will
        % be broken into pieces, each with the same size as one of the variables.
        % Each variable is the replaced by the corresponding piece
        % of obj1.

            if isnumeric(obj)
                return
            end

            if iscell(obj)
                % Apply substitution to each element of cell array 'obj'
                for k=1:length(obj)
                    obj{k}=substitute(toCalculus(obj{k}),var,obj1);
                end
                return
            end

            if isstruct(obj)
                % Apply substitution to each element of structure 'obj'
                fn=fieldnames(obj);
                for k=1:length(fn)
                    obj.(fn{k})=substitute(toCalculus(obj.(fn{k})),var,obj1);
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

            if ~isa(var,'Tcalculus')
                error('2nd argument of substitute must be a Tcalculus object');
            end

            if isa(obj,'Tcalculus')
                obj=substituteRecursive(obj,var,obj1);
            end

        end

        function obj=substituteRecursive(obj,var,obj1)
        % obj=substituteRecursive(obj,var,obj1)
        %
        % Recursive function used by substitute()

            %% Check if already in the cache
            obj_substitute_cache=substitute_cache(obj); % [variable.TCindex,derivative.TCindex;...]
            k=find(obj_substitute_cache(:,1)==var.TCindex & obj_substitute_cache(:,2)==obj1.TCindex);
            if length(k)==1
                %fprintf('.');
                obj=Tcalculus(obj_substitute_cache(k,3));
                return
            elseif length(k)>1
                error('substitute_cache has repeated entries\n');
            end
            originalObj=obj;

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

            %% Add to the cache
            add2substitute_cache(originalObj,var,obj1,obj);
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

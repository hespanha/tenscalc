function [str_obj]=spy(obj,indent)
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

    if nargin<2
        indent=0;
    end

    str_children='';
    for i=1:length(obj.objs)
        str_children=[str_children,spy@Tcalculus(obj.objs{i},indent+1)];
    end

    switch obj.type

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%                        Tokens                            %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case {'zeros'}
        str_obj=sprintf('zeros(%s)',index2str(obj.size));
      case {'ones'}
        str_obj=sprintf('ones(%s)',index2str(obj.size));
      case {'eye'}
        str_obj=sprintf('eye(%s)',index2str(obj.size));

      case {'constant'}
        str_obj=sprintf('constant(%s)',index2str(obj.size));
      case {'variable'}
        str_obj=obj.name;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%                   Referencing, rashaping                 %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case {'subsref'}
        switch obj.inds.type
          case '()',
            str_obj=sprintf('subsref(%s|%s|,?)',obj.objs{1}.type,index2str(size(obj.objs{1})));
          otherwise,
            error('subsref of type ''%s'' not implemented\n',obj.inds.type);
        end

      case {'reshape'}
        str_obj=sprintf('reshape(%s|%s|,[%s])',obj.objs{1}.type,index2str(size(obj.objs{1})),index2str(size(obj)));

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%             Unitary Functions/Operators                  %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case {'ctranspose','transpose','iszero','ispositive','ispositivedefinite','norm2','chol','ldl','lu'}
        str_obj=sprintf('%s(%s|%s|)',obj.type,obj.objs{1}.type,index2str(size(obj.objs{1})));
      case {'sum'}
        str_obj=sprintf('sum(.,%d)',obj.inds);
      case {'diag'}
        if length(size(obj))==1
            % matrix->vector
            str_obj=sprintf('spdiags(.,0)',obj.objs{1}.type,index2str(size(obj.objs{1})));
        else
            % vector->matrix
            str_obj=sprintf('spdiags(%s|%s|,0,%s)',obj.objs{1}.type,index2str(size(obj.objs{1})),index2str(size(obj)));
        end
      case {'repmat'}
        str_obj=sprintf('repmat(%s|%s|,[%s])',obj.objs{1}.type,index2str(size(obj.objs{1})),index2str(obj.inds));
      case {'gradient'}
        str_obj=sprintf('gradient(%s|%s|,%s)',obj.objs{1}.type,index2str(size(obj.objs{1})),obj.objs{2});
      case {'compose'}
        str_fun=char(obj.fun);
        str_children=spy@Tcalculus(obj.objs{1},indent+1);
        str_obj=sprintf('feval(%s,%s|%s|)',str_fun,obj.objs{1}.type,index2str(size(obj.objs{1})));

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%                 Binary operators                         %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case {'pptrs','clp'}
        str_obj=sprintf('%s(%s|%s|,%s|%s|)',obj.type,obj.objs{1}.type,index2str(size(obj.objs{1})),obj.type,obj.objs{2}.type,index2str(size(obj.objs{2})));

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%                 Multi-ary operators                      %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case {'plus'}
        str_obj='plus(';
        for i=1:length(obj.objs)
            if i>1
                str_obj=[str_obj,','];
            end
            if obj.inds{i}>0
                str_obj=sprintf('%s+%s|%s|',str_obj,obj.objs{i}.type,index2str(size(obj.objs{i})));
            else
                str_obj=sprintf('%s-%s|%s|',str_obj,obj.objs{i}.type,index2str(size(obj.objs{i})));
            end
        end
        str_obj=sprintf('%s)',str_obj);
      case {'mtimes'}
        str_obj='';
        for i=1:length(obj.objs)
            if i>1
                str_obj=[str_obj,'*'];
            end
            str_obj=sprintf('%s%s|%s|',str_obj,obj.type,obj.objs{i}.type,index2str(size(obj.objs{i})));
        end

      case {'times'}
        str_obj='';
        for i=1:length(obj.objs)
            if i>1
                str_obj=[str_obj,'.*'];
            end
            str_obj=sprintf('%s%s|%s|',str_obj,obj.type,obj.objs{i}.type,index2str(size(obj.objs{i})));
        end

      case {'rdivide'}
        str_obj='';
        for i=1:length(obj.objs)
            if i>1
                str_obj=[str_obj,'./'];
            end
            str_obj=sprintf('%s%s|%s|',str_obj,obj.type,obj.objs{i}.type,index2str(size(obj.objs{i})));
        end

      case {'mldivide'}
        str_obj='';
        for i=1:length(obj.objs)
            if i>1
                str_obj=[str_obj,'\'];
            end
            str_obj=sprintf('%s%s|%s|',str_obj,obj.type,obj.objs{i}.type,index2str(size(obj.objs{i})));
        end

      case {'tprod'}
        str_obj='tprod(';
        for i=1:length(obj.objs)
            if i>1
                str_obj=[str_obj,','];
            end
            str_obj=sprintf('%s%s|%s|,[%s]',str_obj,obj.objs{i}.type,index2str(size(obj.objs{i})),index2str(obj.inds{i}));
        end
        str_obj=sprintf('%s)',str_obj);

      case {'tprod_matlab'}
        str_obj='tprod_matlab(';
        for i=1:length(obj.objs)
            if i>1
                str_obj=[str_obj,','];
            end
            str_obj=sprintf('%s%s|%s|,[%s]',str_obj,obj.objs{i}.type,index2str(size(obj.objs{i})),index2str(obj.inds{i}));
        end
        str_obj=sprintf('%s)',str_obj);

      case {'cat'}
        str_obj=sprintf('cat(%d',obj.inds);
        for i=1:length(obj.objs)
            str_obj=sprintf('%s,%s|%s|',str_obj,obj.objs{i}.type,index2str(size(obj.objs{i})));
        end
        str_obj=[str_obj,')'];

      otherwise
        error('str(%s) incomplete\n',obj.type)
    end

    str_obj=sprintf('%s [%s]: %s\n%s',repmat('   ',1,indent),index2str(size(obj)),str_obj,str_children);
end

function writeMatlabInstructions(obj,fid,ks)
% Type of code produced:
%    ''Matlab'' - all computations done pure matlab code.
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

    for i=1:length(ks)
        k=ks(i);
        [type,parameters,operands]=getInstruction(int64(k));
        operands=obj.memoryLocations(operands);
        parameters=getArrayFromByteStream(uint8(parameters));
        osize=parameters{1};
        parameters=parameters{2};
        switch (type)
          case obj.Itypes.I_set
          case obj.Itypes.I_load
            if issparse(parameters)
                fprintf(fid,'\t\tobj.m%d=%s; %% [%s]\n',...
                        obj.memoryLocations(k),mat2str_compact(parameters),index2str(osize));
                
            else
                fprintf(fid,'\t\tobj.m%d=[%s]; %% [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),index2str(osize));
            end
                
          case obj.Itypes.I_Mones
            while length(parameters)<2
                parameters(end+1)=1;
            end
            fprintf(fid,'\t\tobj.m%d=ones(%s); %% [%s]\n',...
                    obj.memoryLocations(k),mymat2str(parameters),index2str(osize));
                
          case obj.Itypes.I_Meye
            if length(parameters)==2
                fprintf(fid,'\t\tobj.m%d=speye(%s); %% [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),index2str(osize));
            else
                fprintf(fid,'\t\tobj.m%d=myeye([%s]); %% [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),index2str(osize));
            end
        
          case obj.Itypes.I_Mzeros
            while length(parameters)<2
                parameters(end+1)=1;
            end
            if length(parameters)==2
                fprintf(fid,'\t\tobj.m%d=sparse(%d,%d); %% [%s]\n',...
                        obj.memoryLocations(k),parameters(1),parameters(2),index2str(osize));
            else
                fprintf(fid,'\t\tobj.m%d=zeros(%s); %% [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),index2str(osize));
            end
                
          case obj.Itypes.I_Mrepmat
            while length(parameters)<2
                parameters(end+1)=1;
            end
            if length(parameters~=2)
                fprintf(fid,'\t\tobj.m%d=repmat(full(obj.m%d),[%s]); %% [%s]\n',...
                        obj.memoryLocations(k),operands,index2str(parameters),index2str(osize));
            else
                fprintf(fid,'\t\tobj.m%d=repmat(obj.m%d,[%s]); %% [%s]\n',...
                        obj.memoryLocations(k),operands,index2str(parameters),index2str(osize));
            end
            
          case obj.Itypes.I_Mreshape
            while length(parameters)<2
                parameters(end+1)=1;
            end
            if length(parameters~=2)
                fprintf(fid,'\t\tobj.m%d=reshape(full(obj.m%d),[%s]); %% [%s]\n',...
                        obj.memoryLocations(k),operands,index2str(parameters),index2str(osize));
            else
                fprintf(fid,'\t\tobj.m%d=reshape(obj.m%d,[%s]); %% [%s]\n',...
                        obj.memoryLocations(k),operands,index2str(parameters),index2str(osize));
            end
                    
          case obj.Itypes.I_Mmin
            fprintf(fid,'\t\tobj.m%d=min(obj.m%d,[],%s);\n',...
                    obj.memoryLocations(k),operands,mymat2str(parameters));
            
          case obj.Itypes.I_Msum
            fprintf(fid,'\t\tobj.m%d=sum(obj.m%d,%d);\n',...
                    obj.memoryLocations(k),operands,parameters(1));
            for i=2:length(parameters)
                fprintf(fid,'\t\tobj.m%d=obj.m%d+sum(obj.m%d,%d);\n',...
                        obj.memoryLocations(k),obj.memoryLocations(k),operands,parameters(i));
            end
            sz=osize;while (length(sz)<2);sz(end+1)=1;end
            fprintf(fid,'\t\tobj.m%d=reshape(obj.m%d,[%s]); %% [%s]\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),mymat2str(sz),index2str(osize));
            
          case obj.Itypes.I_Mclp
            fprintf(fid,'\t\tobj.m%d=clp(obj.m%d,obj.m%d); %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mcat
            fprintf(fid,'\t\tobj.m%d=cat(%s',obj.memoryLocations(k),mymat2str(parameters));
            fprintf(fid,',obj.m%d',operands);
            fprintf(fid,'); %% [%s]\n',index2str(osize));
            
          case obj.Itypes.I_Mnorm2
            fprintf(fid,'\t\tobj.m%d=sum(obj.m%d(:).^2); %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mnorm1
            fprintf(fid,'\t\tobj.m%d=sum(abs(obj.m%d(:))); %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mnorminf
            fprintf(fid,'\t\tobj.m%d=max(abs(obj.m%d(:))); %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mfull
            fprintf(fid,'\t\tobj.m%d=full(obj.m%d); %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mdiag
            fprintf(fid,'\t\tobj.m%d=spdiags(obj.m%d,0,%d,%d); %% [%s]\n',...
                    obj.memoryLocations(k),operands,osize(1),osize(2),index2str(osize));
            
          case obj.Itypes.I_Mcompose
            functions={'@(x__)log(x__)';'@(x__)1./x__';'@(x__)-1./x__.^2';
                       '@(x__)2./x__.^3';
                       '@(x__)exp(x__)';
                       '@(x__)sin(x__)';'@(x__)-sin(x__)';
                       '@(x__)cos(x__)';'@(x__)-cos(x__)';
                       '@(x__)abs(x__)';'@(x__)sign(x__)';
                       '@(x__)atan(x__)';'@(x__)1./(1+x__.^2)';'@(x__)-2*x__./(1+x__.^2).^2';
                       '@(x__)sqrt(x__)';'@(x__).5./sqrt(x__)';'@(x__)-.25./x__.^1.5';
                       '@(x__)x__.^2';'@(x__)2*x__';'@(x__)2*ones(size(x__))';
                       '@(x__)x__.^3';'@(x__)3*x__.^2';'@(x__)6*x__';
                       };
            q=find(strcmp(parameters,functions));
            if isempty(q)
                error('computeScalarInstructions: compose object not implemented for function ''%s''\n',parameters);
            end
            call=regexprep(functions{q}(7:end),'x__',sprintf('obj.m%d',operands));
            fprintf(fid,'\t\tobj.m%d=%s; %% [%s]\n',obj.memoryLocations(k),call,index2str(osize));
            
          case obj.Itypes.I_Msubsref
            switch parameters.type 
              case '()',
                str='(';
                for i=1:length(parameters.subs)
                    if i ~= 1
                        str=[str,','];
                    end
                    if ischar(parameters.subs{i}) && parameters.subs{i}==':'
                        str=sprintf('%s:',str);
                    else
                        str=sprintf('%s%s',str,mat2str_compact(parameters.subs{i}(:)));
                    end
                end
                str=sprintf('%s)',str);
              otherwise,
                error('subsref of type ''%s'' not implemented\n',parameters.type);
            end
            fprintf(fid,'\t\tobj.m%d=obj.m%d%s; %% [%s]\n',...
                    obj.memoryLocations(k),operands,str,index2str(osize));
            
          case obj.Itypes.I_Mctranspose
            fprintf(fid,'\t\tobj.m%d=obj.m%d''; %% [%s]\n',obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mmtimes
            fprintf(fid,'\t\tobj.m%d=obj.m%d*obj.m%d; %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mtimes
            fprintf(fid,'\t\tobj.m%d=obj.m%d.*obj.m%d; %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mrdivide
            fprintf(fid,'\t\tobj.m%d=obj.m%d./obj.m%d; %% [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(osize));
            
          case {obj.Itypes.I_Mtprod,obj.Itypes.I_Mtprod_matlab}
            fprintf(fid,'\t\tobj.m%d=mytprod',obj.memoryLocations(k));
            sep='(';
            parameters=parameters(2:end); %% only keep indices 
            for i=1:length(parameters)
                fprintf(fid,'%cobj.m%d,[%s]',sep,operands(i),index2str(parameters{i}));
                sep=',';
            end
            fprintf(fid,'); %% [%s]\n',index2str(osize));
            
          case obj.Itypes.I_Mplus
            fprintf(fid,'\t\tobj.m%d=',obj.memoryLocations(k));
            for i=1:length(parameters)
                if parameters(i)>0
                    fprintf(fid,'+obj.m%d',operands(i));
                else
                    fprintf(fid,'-obj.m%d',operands(i));
                end
            end
            fprintf(fid,'; %% [%s]\n',index2str(osize));
            
            %% LU factorization, et al.
            % WW=randn(6,6);WW=sparse(WW);b=randn(6,1);
            % dxl=WW\b
            % [L,U,p,q,r]=lu(WW,'vector');s=inv(r);  %% I_Mlu
            % m=s*b;m=L\m(p);               %% I_Mmldivide_l1
            % m=U\m;m(q)=m                  %% I_Mmldivide_u
            
          case obj.Itypes.I_Mlu
            fprintf(fid,'\t\tobj.m%d=struct();\n',obj.memoryLocations(k));
            fprintf(fid,'\t\t[obj.m%d.L,obj.m%d.U,obj.m%d.p,obj.m%d.q,obj.m%d.r]=lu(sparse(obj.m%d),''vector'');\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),...
                    obj.memoryLocations(k),obj.memoryLocations(k),obj.memoryLocations(k),operands);
            fprintf(fid,'\t\tobj.m%d.s=inv(obj.m%d.r); %% [%s]\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),index2str(osize));
            
          case obj.Itypes.I_Mmldivide_l1
            fprintf(fid,'\t\tobj.m%d=obj.m%d.s*obj.m%d;\n',...
                    obj.memoryLocations(k),operands(1),operands(2));            %%% problem for ldl factorization
            fprintf(fid,'\t\tobj.m%d=obj.m%d.L\\(obj.m%d(obj.m%d.p)); %% [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),operands(1),index2str(osize));
            
          case obj.Itypes.I_Mmldivide_u
            fprintf(fid,'\t\tobj.m%d=obj.m%d.U\\(obj.m%d);\n',...
                    obj.memoryLocations(k),operands(1),operands(2));
            fprintf(fid,'\t\tobj.m%d(obj.m%d.q)=obj.m%d; %% [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),index2str(osize));
            
            %% LDL factorization, et al.
            % WW=randn(6,6);WW=sparse(WW+WW');b=randn(6,1);
            % dxl=WW\b
            % [L,D,p,s]=ldl(WW,'vector');  %% I_Mldl
            % m=s*b;m=L\m(p);              %% I_Mmldivide_l1
            % m=D\m;                       %% I_Mmldivide_d
            % m=L'\m;m(p)=m;m=s*m;         %% I_Mmldivide_u1
            
          case obj.Itypes.I_Mldl
            fprintf(fid,'\t\tobj.m%d=struct();\n',obj.memoryLocations(k));
            fprintf(fid,'\t\t[obj.m%d.L,obj.m%d.D,obj.m%d.p,obj.m%d.s]=ldl(sparse(obj.m%d),''vector''); %% [%s]\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),obj.memoryLocations(k),obj.memoryLocations(k),operands,index2str(osize));
            
          case obj.Itypes.I_Mmldivide_u1
            fprintf(fid,'\t\tobj.m%d=obj.m%d.L''\\(obj.m%d);\n',...
                    obj.memoryLocations(k),operands(1),operands(2));
            fprintf(fid,'\t\tobj.m%d(obj.m%d.p)=obj.m%d;\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k));
            fprintf(fid,'\t\tobj.m%d=obj.m%d.s*obj.m%d; %% [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),index2str(osize));
                  
          case obj.Itypes.I_Mmldivide_d
            fprintf(fid,'\t\tobj.m%d=obj.m%d.D\\(obj.m%d); %% [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(2),index2str(osize));
                  
      otherwise
        instructionTypes
        error('writeMatlabInstructions: instruction %d not implemented (see instructionTypes)\n',type);
    end
end

end
    
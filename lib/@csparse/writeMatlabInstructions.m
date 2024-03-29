function writeMatlabInstructions(obj,fid,ks)
% Type of code produced:
%    ''Matlab'' - all computations done pure matlab code.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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
                fprintf(fid,'\t\tobj.m%d=%s; %% op %d: [%s]\n',...
                        obj.memoryLocations(k),mat2str_compact(parameters),k,index2str(osize));

            else
                fprintf(fid,'\t\tobj.m%d=[%s]; %% op %d: [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),k,index2str(osize));
            end

          case obj.Itypes.I_Mones
            while length(parameters)<2
                parameters(end+1)=1;
            end
            fprintf(fid,'\t\tobj.m%d=ones(%s); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),mymat2str(parameters),k,index2str(osize));

          case obj.Itypes.I_Meye
            if length(parameters)==2
                fprintf(fid,'\t\tobj.m%d=speye(%s); %% op %d: [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),k,index2str(osize));
            else
                fprintf(fid,'\t\tobj.m%d=myeye([%s]); %% op %d: [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),k,index2str(osize));
            end

          case obj.Itypes.I_Mzeros
            while length(parameters)<2
                parameters(end+1)=1;
            end
            if length(parameters)==2
                fprintf(fid,'\t\tobj.m%d=sparse(%d,%d); %% op %d: [%s]\n',...
                        obj.memoryLocations(k),parameters(1),parameters(2),k,index2str(osize));
            else
                fprintf(fid,'\t\tobj.m%d=zeros(%s); %% op %d: [%s]\n',...
                        obj.memoryLocations(k),mymat2str(parameters),k,index2str(osize));
            end

          case obj.Itypes.I_Mrepmat
            while length(parameters)<2
                parameters(end+1)=1;
            end
            if false && length(parameters>2)
                fprintf(fid,'\t\tobj.m%d=repmat(full(obj.m%d),[%s]); %% op %d: [%s]\n',...
                        obj.memoryLocations(k),operands,index2str(parameters),k,index2str(osize));
            else
                fprintf(fid,'\t\tobj.m%d=repmat(obj.m%d,[%s]); %% op %d: [%s]\n',...
                        obj.memoryLocations(k),operands,index2str(parameters),k,index2str(osize));
            end

          case obj.Itypes.I_Mreshape
            while length(parameters)<2
                parameters(end+1)=1;
            end
            if length(parameters)>2 % ND sparse arrays are not supported
                fprintf(fid,'\t\tobj.m%d=reshape(full(obj.m%d),[%s]); %% op %d: [%s]\n',...
                        obj.memoryLocations(k),operands,index2str(parameters),k,index2str(osize));
            else
                %[t,~,o]=getInstruction(int64(operands))
                if prod(osize)>100
                    fprintf(fid,'\t\tobj.m%d=sparse(reshape(obj.m%d,[%s])); %% op %d: [%s]\n',...
                            obj.memoryLocations(k),operands,index2str(parameters),k,index2str(osize));
                else
                    fprintf(fid,'\t\tobj.m%d=reshape(obj.m%d,[%s]); %% op %d: [%s]\n',...
                            obj.memoryLocations(k),operands,index2str(parameters),k,index2str(osize));
                end
            end

          case obj.Itypes.I_Mvec2tensor
            fprintf(fid,'\t\tobj.m%d=vec2tensor(obj.m%d,[%s],[%s],%d); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,...
                    index2str(parameters{2}),index2str(parameters{3}'),parameters{1},k,index2str(osize));

          case obj.Itypes.I_Mpermute_matlab
            fprintf(fid,'\t\tobj.m%d=permute(obj.m%d,[%s]); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,index2str(parameters),k,index2str(osize));

          case obj.Itypes.I_Mmin
            fprintf(fid,'\t\tobj.m%d=min(obj.m%d,[],%s);\n',...
                    obj.memoryLocations(k),operands,mymat2str(parameters));

          case obj.Itypes.I_Mmax
            fprintf(fid,'\t\tobj.m%d=max(obj.m%d,[],%s);\n',...
                    obj.memoryLocations(k),operands,mymat2str(parameters));

          case obj.Itypes.I_Mmin2
            fprintf(fid,'\t\tobj.m%d=min(obj.m%d,obj.m%d);\n',...
                    obj.memoryLocations(k),operands);

          case obj.Itypes.I_Mmax2
            fprintf(fid,'\t\tobj.m%d=max(obj.m%d,obj.m%d);\n',...
                    obj.memoryLocations(k),operands);

          case obj.Itypes.I_Msum
            fprintf(fid,'\t\tobj.m%d=sum(obj.m%d,%d);\n',...
                    obj.memoryLocations(k),operands,parameters(1));
            for i=2:length(parameters)
                fprintf(fid,'\t\tobj.m%d=obj.m%d+sum(obj.m%d,%d);\n',...
                        obj.memoryLocations(k),obj.memoryLocations(k),operands,parameters(i));
            end
            sz=osize;while (length(sz)<2);sz(end+1)=1;end
            fprintf(fid,'\t\tobj.m%d=reshape(obj.m%d,[%s]); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),mymat2str(sz),k,index2str(osize));

          case obj.Itypes.I_Mclp
            fprintf(fid,'\t\tobj.m%d=clp(obj.m%d,obj.m%d); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mcat
            fprintf(fid,'\t\tobj.m%d=cat(%s',obj.memoryLocations(k),mymat2str(parameters));
            fprintf(fid,',obj.m%d',operands);
            fprintf(fid,'); %% op %d: [%s]\n',k,index2str(osize));

          case obj.Itypes.I_Mnorm2
            fprintf(fid,'\t\tobj.m%d=sum(obj.m%d(:).^2); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mnorm1
            fprintf(fid,'\t\tobj.m%d=sum(abs(obj.m%d(:))); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mnorminf
            fprintf(fid,'\t\tobj.m%d=max(abs(obj.m%d(:))); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mfull
            fprintf(fid,'\t\tobj.m%d=full(obj.m%d); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mdiag
            fprintf(fid,'\t\tobj.m%d=spdiags(obj.m%d,0,%d,%d); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,osize(1),osize(2),k,index2str(osize));

          case obj.Itypes.I_Mcomponentwise
            mfun=char(parameters{1});
            s=regexp(mfun,'^@\(([^)]+)\)(.+)$','tokens');
            if isempty(s)
                call=sprintf('%s(obj.m%d)',mfun,operands);
            else
                call=regexprep(s{1}{2},['\<',s{1}{1},'\>'],sprintf('obj.m%d',operands));
            end
            fprintf(fid,'\t\tobj.m%d=%s; %% op %d: [%s]\n',obj.memoryLocations(k),call,k,index2str(osize));

          case obj.Itypes.I_Mcompose
            functions={'@(x__)log(x__)';'@(x__)1./x__';'@(x__)-1./x__.^2';
                       '@(x__)2./x__.^3';
                       '@(x__)exp(x__)';
                       '@(x__)sin(x__)';'@(x__)-sin(x__)';
                       '@(x__)cos(x__)';'@(x__)-cos(x__)';
                       '@(x__)round(x__)';'@(x__)ceil(x__)';'@(x__)floor(x__)';
                       '@(x__)abs(x__)';'@(x__)sign(x__)';
                       '@(x__)atan(x__)';'@(x__)1./(1+x__.^2)';'@(x__)-2*x__./(1+x__.^2).^2';
                       '@(x__)sqrt(x__)';'@(x__).5./sqrt(x__)';'@(x__)-.25./x__.^1.5';
                       '@(x__)x__.^2';'@(x__)2*x__';'@(x__)2*ones(size(x__))';
                       '@(x__)x__.^3';'@(x__)3*x__.^2';'@(x__)6*x__';
                       '@(x__)log(1+x__)/log(2)';'@(x__)1./(1+x__)/log(2)';'@(x__)-1./(1+x__).^2';
                       '@(x__)log(1+exp(x__))';'@(x__)1./(1+exp(-x__))';...
                                               '@(x__)1./(2+exp(-x__)+exp(x__))';
                       '@(x__)relu(x__)';'@(x__)heaviside(x__)';
                      };
            q=find(strcmp(parameters,functions));
            if isempty(q)
                error('computeScalarInstructions: compose object not implemented for function ''%s''\n',parameters);
            end
            call=regexprep(functions{q}(7:end),'x__',sprintf('obj.m%d',operands));
            fprintf(fid,'\t\tobj.m%d=%s; %% op %d: [%s]\n',obj.memoryLocations(k),call,k,index2str(osize));

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
                        str=sprintf('%s%s',str,mat2str_compact(parameters.subs{i}));
                    end
                end
                if length(parameters.subs)<2
                    str=sprintf('%s,1)',str);                    
                else
                    str=sprintf('%s)',str);
                end
              case '(:)',
                % compressed subscript as stored in computeMatlabInstructions
                str='(';
                for i=1:length(parameters.subs)
                    if i ~= 1
                        str=[str,','];
                    end
                    if ischar(parameters.subs{i}) && parameters.subs{i}==':'
                        str=sprintf('%s:',str);
                    else
                        str=sprintf('%s%d:%d',str,parameters.subs{i});
                    end
                end
                if length(parameters.subs)<2
                    str=sprintf('%s,1)',str);                    
                else
                    str=sprintf('%s)',str);
                end
              otherwise,
                error('subsref of type ''%s'' not implemented\n',parameters.type);
            end
            fprintf(fid,'\t\tobj.m%d=obj.m%d%s; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,str,k,index2str(osize));

          case obj.Itypes.I_Mctranspose
            fprintf(fid,'\t\tobj.m%d=obj.m%d''; %% op %d: [%s]\n',obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mmtimes
            fprintf(fid,'\t\tobj.m%d=obj.m%d*obj.m%d; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mtimes
            fprintf(fid,'\t\tobj.m%d=obj.m%d.*obj.m%d; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case obj.Itypes.I_Mrdivide
            fprintf(fid,'\t\tobj.m%d=obj.m%d./obj.m%d; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands,k,index2str(osize));

          case {obj.Itypes.I_Mtprod,obj.Itypes.I_Mtprod_matlab}
            fprintf(fid,'\t\tobj.m%d=mytprod',obj.memoryLocations(k));
            sep='(';
            parameters=parameters(2:end); %% only keep indices
            for i=1:length(parameters)
                fprintf(fid,'%cobj.m%d,[%s]',sep,operands(i),index2str(parameters{i}));
                sep=',';
            end
            fprintf(fid,'); %% op %d: [%s]\n',k,index2str(osize));

          case obj.Itypes.I_Mplus
            fprintf(fid,'\t\tobj.m%d=',obj.memoryLocations(k));
            for i=1:length(parameters)
                if parameters(i)>0
                    fprintf(fid,'+obj.m%d',operands(i));
                else
                    fprintf(fid,'-obj.m%d',operands(i));
                end
            end
            fprintf(fid,'; %% op %d: [%s]\n',k,index2str(osize));

            %% UMFpack LU factorization
            % [L,U,P,Q,D] = umfpack (A)
            % L*U=P*(D\WW)*Q
            %    WW * dxl = b
            %    P*inv(D)*WW*Q*inv(Q)*dxl = P*inv(D)*b
            %    L*U*inv(Q)*dxl = P*inv(D)*b
            %    dxl=Q*inv(U)*inv(L)*P*inv(D)*b
            %    dxl=Q*(U\(L\(P*(D\b))))

          case obj.Itypes.I_Mlu_sym
            fprintf(fid,'\t\tobj.m%d=struct();\n',obj.memoryLocations(k));
            fprintf(fid,'\t\t[obj.m%d.L,obj.m%d.U,obj.m%d.P,obj.m%d.Q,obj.m%d.D]=umfpack(sparse(obj.m%d));\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),...
                    obj.memoryLocations(k),obj.memoryLocations(k),obj.memoryLocations(k),operands);

          case obj.Itypes.I_Mmldivide_lu
            fprintf(fid,'\t\tobj.m%d=obj.m%d.Q*(obj.m%d.U\\(obj.m%d.L\\(obj.m%d.P*(obj.m%d.D\\obj.m%d))));\n',...
                    obj.memoryLocations(k),operands(1),operands(1),operands(1),operands(1),operands(1),operands(2));            %%% problem for ldl factorization

            %% LU factorization, et al.
            % WW=randn(6,6);WW=sparse(WW);b=randn(6,1);
            %
            % [L,U,p,q,d]=lu(WW,'vector');s=inv(d);  %% I_Mlu
            % L*U=d(:,p)\WW(:,q);  % per MATLAB's help
            %
            %    WW * dxl = b
            %    WW(:,q) * dxl(q) = b
            %    inv(d(:,p)) * WW(:,q) * dxl(q) = inv(d(:,p)) * b
            %        d * s = I <=> d(:,p) * s(p,:) = I => inv(d(:,p))=s(p,:)
            %    s(p,:) * WW(:,q) * dxl(q) = s(p,:) * b
            %    L * U * dxl(q) = s(p,:) * b
            %
            % m=s*b;m=L\m(p,:);               %% I_Mmldivide_l1
            % m=U\m;m(q,:)=m                  %% I_Mmldivide_u
            %

          case obj.Itypes.I_Mlu
            fprintf(fid,'\t\tobj.m%d=struct();\n',obj.memoryLocations(k));
            fprintf(fid,'\t\t[obj.m%d.L,obj.m%d.U,obj.m%d.p,obj.m%d.q,obj.m%d.d]=lu(sparse(obj.m%d),''vector'');\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),...
                    obj.memoryLocations(k),obj.memoryLocations(k),obj.memoryLocations(k),operands);
            fprintf(fid,'\t\tobj.m%d.s=inv(obj.m%d.d); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),k,index2str(osize));

          case obj.Itypes.I_Mmldivide_l1
            fprintf(fid,'\t\tobj.m%d=obj.m%d.s*obj.m%d;\n',...
                    obj.memoryLocations(k),operands(1),operands(2));            %%% problem for ldl factorization
            fprintf(fid,'\t\tobj.m%d=obj.m%d.L\\(obj.m%d(obj.m%d.p,:)); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),operands(1),k,index2str(osize));

          case obj.Itypes.I_Mmldivide_u
            fprintf(fid,'\t\tobj.m%d=obj.m%d.U\\(obj.m%d);\n',...
                    obj.memoryLocations(k),operands(1),operands(2));
            fprintf(fid,'\t\tobj.m%d(obj.m%d.q,:)=obj.m%d; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),k,index2str(osize));

            %% LU factorization
            % A = D*P'*L*U*Q'
          case obj.Itypes.I_Mlu_u
            fprintf(fid,'\t\tobj.m%d(:,obj.m%d.q)=obj.m%d.U; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(1),k,index2str(osize));

          case obj.Itypes.I_Mlu_l
            fprintf(fid,'\t\tobj.m%d(obj.m%d.p,:)=obj.m%d.L; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(1),k,index2str(osize));
            fprintf(fid,'\t\tobj.m%d=obj.m%d.d*obj.m%d; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),k,index2str(osize));

          case obj.Itypes.I_Mlu_d
            fprintf(fid,'\t\tobj.m%d=diag(obj.m%d.U); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),k,index2str(osize));

            %% LDL factorization, et al.
            % WW=randn(6,6);WW=sparse(WW+WW');b=randn(6,1);
            % dxl=WW\b
            % [L,D,p,s]=ldl(WW,'vector');  %% I_Mldl
            % m=s*b;m=L\m(p,:);              %% I_Mmldivide_l1
            % m=D\m;                       %% I_Mmldivide_d
            % m=L'\m;m(p,:)=m;m=s*m;         %% I_Mmldivide_u1

          case obj.Itypes.I_Mldl
            fprintf(fid,'\t\tobj.m%d=struct();\n',obj.memoryLocations(k));
            % THRESH as the pivot tolerance in the algorithm. THRESH
            % must be a double scalar lying in the interval [0,
            % 0.5]. The default value for THRESH is 0.01. Using
            % smaller values of THRESH may give faster factorization
            % times and fewer entries, but may also result in a less
            % stable factorization.
            fprintf(fid,'\t\t[obj.m%d.L,obj.m%d.D,obj.m%d.p,obj.m%d.s]=ldl(sparse(obj.m%d),%g,''vector''); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),obj.memoryLocations(k),obj.memoryLocations(k),obj.memoryLocations(k),...
                    operands,obj.LDLthreshold,k,index2str(osize));

          case obj.Itypes.I_Mmldivide_u1
            fprintf(fid,'\t\tobj.m%d=obj.m%d.L''\\(obj.m%d);\n',...
                    obj.memoryLocations(k),operands(1),operands(2));
            fprintf(fid,'\t\tobj.m%d(obj.m%d.p,:)=obj.m%d;\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k));
            fprintf(fid,'\t\tobj.m%d=obj.m%d.s*obj.m%d; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),k,index2str(osize));

          case obj.Itypes.I_Mmldivide_d
            fprintf(fid,'\t\tobj.m%d=obj.m%d.D\\(obj.m%d); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(2),k,index2str(osize));

            %% LDL factorization
            % A = inv(s)*P*L*D*L'*P'*inv(s)
          case obj.Itypes.I_Mldl_d
            fprintf(fid,'\t\tobj.m%d=diag(obj.m%d.D); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),k,index2str(osize));

          case obj.Itypes.I_Mldl_l
            fprintf(fid,'\t\tobj.m%d(obj.m%d.p,:)=obj.m%d.L; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(1),k,index2str(osize));
            fprintf(fid,'\t\tobj.m%d=inv(obj.m%d.s)*obj.m%d; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),obj.memoryLocations(k),k,index2str(osize));

            %% det/logdet
            % [L,D,p,s]=ldl(A,'vector');                  %% I_Mldl
            % log(det(A))=log(det(D)/det(s)^2);           %% I_logdet_ldl
            %            =sum(log(diag(D)))-sum(log(diag(s.^2))) % if diag(D)>0 & D diagonal

            % [L,U,p,q,d]=lu(A,'vector');s=inv(d);        %% I_Mlu
            % P(p,p)=eye(size(A));Q(q,q)=eye(size(A));
            % log(det(A))=log(prod(diag(U))*prod(diag(d)))*det(P)*det(Q);  %% I_logdet_lu
            %            =sum(log(diag(U).*diag(d)))+log(det(P)*det(Q))
          case obj.Itypes.I_Mlogdet_ldl
            fprintf(fid,'\t\t if isdiag(obj.m%d.D)\n',operands(1));
            % only if D is really diagonal
            fprintf(fid,'\t\t\tobj.m%d=sum(log(abs(diag(obj.m%d.D))))-sum(log(diag(obj.m%d.s).^2)); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(1),k,index2str(osize));
            fprintf(fid,'else\n');
            % general case
            fprintf(fid,'\t\t\tobj.m%d=log(abs(det(obj.m%d.D)))-sum(log(diag(obj.m%d.s).^2)); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(1),k,index2str(osize));
            fprintf(fid,'end\n');
          case obj.Itypes.I_Mdet_ldl
            fprintf(fid,'\t\tobj.m%d=det(obj.m%d.D)/det(obj.m%d.s)^2; %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(1),k,index2str(osize));

          case obj.Itypes.I_Mlogdet_lu
            fprintf(fid,'\t\tobj.m%d.P(obj.m%d.p,1:length(obj.m%d.p))=speye(length(obj.m%d.p)); %% op %d: [%s]\n',...
                    operands(1),operands(1),operands(1),operands(1),k,index2str(osize));
            fprintf(fid,'\t\tobj.m%d.Q(obj.m%d.q,1:length(obj.m%d.q))=speye(length(obj.m%d.q)); %% op %d: [%s]\n',...
                    operands(1),operands(1),operands(1),operands(1),k,index2str(osize));
            if 0
                % only if U is really diagonal with all positive entries
fprintf(fid,'\t\tobj.m%d=sum(log(diag(obj.m%d.U).*diag(obj.m%d.d)))+log(det(obj.m%d.P)*det(obj.m%d.Q)); %% op %d: [%s]\n',...
                         obj.memoryLocations(k),operands(1),operands(1),operands(1),operands(1),k,index2str(osize));
            else
fprintf(fid,'\t\tobj.m%d=sum(log(abs(diag(obj.m%d.U).*diag(obj.m%d.d))))+log(abs(det(obj.m%d.P)*det(obj.m%d.Q))); %% op %d: [%s]\n',...
                         obj.memoryLocations(k),operands(1),operands(1),operands(1),operands(1),k,index2str(osize));
            end
          case obj.Itypes.I_Mdet_lu
            fprintf(fid,'\t\tobj.m%d.P(obj.m%d.p,1:length(obj.m%d.p))=speye(length(obj.m%d.p)); %% op %d: [%s]\n',...
                    operands(1),operands(1),operands(1),operands(1),k,index2str(osize));
            fprintf(fid,'\t\tobj.m%d.Q(obj.m%d.q,1:length(obj.m%d.q))=speye(length(obj.m%d.q)); %% op %d: [%s]\n',...
                    operands(1),operands(1),operands(1),operands(1),k,index2str(osize));
            fprintf(fid,'\t\tobj.m%d=prod(diag(obj.m%d.U).*diag(obj.m%d.d))*det(obj.m%d.P)*det(obj.m%d.Q); %% op %d: [%s]\n',...
                    obj.memoryLocations(k),operands(1),operands(1),operands(1),operands(1),k,index2str(osize));

          otherwise
            instructionTypes
            error('writeMatlabInstructions: instruction %d not implemented (see instructionTypes)\n',type);
        end
    end
end

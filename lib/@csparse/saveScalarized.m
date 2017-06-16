function saveScalarized(obj,filename)
% saveScalarized(obj,filename)
%
% Saves the csparse object as a compuational graph. See computation
% graphs for the file format.
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

    t0=clock;

    fprintf('saveScalarize...');

    intType='int32';

    function const=addConstant(value,name,description)
    % Adds constant to array of constants to eventually by added to
    % .cgc constants file.
    % Returns the record number within .cgc constants file.
        const=length(constants);  % 0-index
        constants(const+1,1)=struct('value',value,'name',name,'description',description);
    end


    function var=addConstantVariable(value,name,description)
    % First adds constant to the array of constants to eventually by
    % added to .cgc constants file and then adds a function node
    % with the corresponding constant.
    % Returns the variable index.

        nVariables=nVariables+1;
        var=nVariables;
        const=addConstant(value,name,description);
        % write constant function type, constant index, output variable
        fwrite(fg,[const_funct,const,var-1],intType); % 0-based
        %[const_funct,const,var]

        % symbols
        if ~isempty(name) || ~isempty(description)
            fwrite(fs,[symfun,nFuncts,length(name),length(description)],intType); % 0-based
            fwrite(fs,[name,description],'char');
            fwrite(fs,[symvar,var-1,length(name),length(description)],intType); % 0-based
            fwrite(fs,[name,description],'char');
            nSym=nSym+2;
        end
        nFuncts=nFuncts+1;
    end

    %% get types
    CG=CGregistration();
    magics={CG.magic(:).name};

    functs={CG.function(:).name};
    functs=regexprep(functs,'^S','I_');  % replace leading S by I_
    ftypes=int64([CG.function(:).type]);

    consts={CG.constant(:).name};
    ctypes=int64([CG.constant(:).type]);

    ios={CG.io(:).name};
    iotypes=int64([CG.io(:).type]);

    symbols={CG.symbol(:).name};
    symtypes=int64([CG.symbol(:).type]);
    symcon=symtypes(strcmp('constant',symbols));
    symfun=symtypes(strcmp('function',symbols));
    symvar=symtypes(strcmp('variable',symbols));
    symio=symtypes(strcmp('io',symbols));

    const_funct=ftypes(strcmp('constant',functs));
    s2ft_funct=ftypes(strcmp('scalars2ftensor',functs));
    ft2s_funct=ftypes(strcmp('ftensor2scalars',functs));

    %% get types for obj.Itypes
    fnames=fieldnames(obj.Itypes);
    iTypes=nan(length(fnames),1);
    for i=1:length(fnames)
        if fnames{i}(3)=='M'
            continue;
        elseif strcmp(fnames{i},'I_set')
            iTypes(getfield(obj.Itypes,fnames{i}))=NaN;
        elseif strcmp(fnames{i},'I_load')
            iTypes(getfield(obj.Itypes,fnames{i}))=const_funct;
        else
            k=find(strcmp(fnames{i},functs));
            if isempty(k)
                functs
                error('csparse/saveScalarized: type ''%s'' not registered\n',fnames{i});
            end
            iTypes(getfield(obj.Itypes,fnames{i}))=ftypes(k);
        end
    end
    %iTypes
    
    %% compute scalar instructions
    fprintf('  computeScalarInstructions... ');t0=clock;
    checkVariableSets(obj);
    computeScalarInstructions(obj,1:height(obj.vectorizedOperations));
    fprintf('done (%.3f sec)\n',etime(clock(),t0));

    %% open files & write magic
    [~,~,ext]=fileparts(filename);
    if ~isempty(ext)
        error('csparse/saveScalarized: filename "%s" must not include extension\n',filename);
    end
    fg=fopen(sprintf('%s.cg',filename),'w');
    fc=fopen(sprintf('%s.cgc',filename),'w');
    fio=fopen(sprintf('%s.cgio',filename),'w');
    fs=fopen(sprintf('%s.cgs',filename),'w');
    switch intType
      case 'int64'
        magic=CG.magic(strcmp('int64l',magics)).type;
      case 'int32'
        magic=CG.magic(strcmp('int32l',magics)).type;
      otherwise
        error('unimplemented intType "%s"\n',intType);
    end
    fwrite(fg,magic,'uchar');
    fwrite(fc,magic,'uchar');
    fwrite(fio,magic,'uchar');
    fwrite(fs,magic,'uchar');

    %% initialize constants
    nVariables=int64(instructionsTableHeight());
    constants=struct('value',{},'name',{},'description',{});

    nFuncts=0;
    nSym=0;
    %% Save functions graph (and the symbolic names of function and variable nodes)
    for thisExp=1:instructionsTableHeight()
        if mod(thisExp,1000)==0
            fprintf('%d/%d\n',thisExp,instructionsTableHeight());
        end
        [type,parameters,operands]=getInstruction(int64(thisExp)); % 1-based operands

        switch (type)
          case {obj.Itypes.I_plus_minus_dot,...
                obj.Itypes.I_inv,...
                obj.Itypes.I_minus_inv_sqr,...
                obj.Itypes.I_div,...
                obj.Itypes.I_minus_dot,...
                obj.Itypes.I_minus_dot_div,...
                obj.Itypes.I_plus_minus_dot,...
                obj.Itypes.I_plus_sqr,...
                obj.Itypes.I_plus_minus_dot_div,...
                obj.Itypes.I_min,...
                obj.Itypes.I_min0,...
                obj.Itypes.I_max,...
                obj.Itypes.I_max0,...
                obj.Itypes.I_max_abs,...
                obj.Itypes.I_clp,...
                obj.Itypes.I_exp,...
                obj.Itypes.I_log,...
                obj.Itypes.I_cos,...
                obj.Itypes.I_minus_cos,...
                obj.Itypes.I_sin,...
                obj.Itypes.I_minus_sin,...
                obj.Itypes.I_abs,...
                obj.Itypes.I_sqrt,...
                obj.Itypes.I_Dsqrt,...
                obj.Itypes.I_DDsqrt,...
                obj.Itypes.I_sqr,...
                obj.Itypes.I_2times,...
                obj.Itypes.I_2,...
                obj.Itypes.I_atan,...
                obj.Itypes.I_Datan,...
                obj.Itypes.I_DDatan...
                }
            % function type, n inputs, 1 output
            fwrite(fg,[iTypes(type),length(operands),1,operands-1,thisExp-1],intType);  % 0-based;
            %[iTypes(type),length(operands),operands-1,thisExp-1]
          
          case {obj.Itypes.I_sum,...
                obj.Itypes.I_sumprod...
                }
            var=addConstantVariable(int64(parameters(:)),'','sum');
            % function type, 1+n inputs, 1 output
            fwrite(fg,[iTypes(type),1+length(operands),1,var-1,operands-1,thisExp-1],intType);  % 0-based;
            %[iTypes(type),1+length(operands),[var-1,operands-1],thisExp-1]
            
          case obj.Itypes.I_set
            % do nothing
            nFuncts=nFuncts-1;
          
          case obj.Itypes.I_load
            var=addConstant(double(parameters),'','load');
            % write constant function type, constant index, output variable
            fwrite(fg,[const_funct,var,thisExp-1],intType);  % 0-based;
            %[const_funct,var-1,thisExp-1]
            
          otherwise
            error('csparse/saveScalarized: missing type %d\n',type);
        end
        nFuncts=nFuncts+1;
    end

    %% Save I/O functions (and their symbolic names)
    ioget=iotypes(strcmp('get',ios));
    ioset=iotypes(strcmp('set',ios));
    iocopy=iotypes(strcmp('copy',ios));
    nIOs=0;
    for i=1:length(obj.gets)
        % add functions to collect data into full tensors
        newvars=[nVariables+1:nVariables+length(obj.gets(i).source)];
        nVariables=newvars(end);
        for j=1:length(obj.gets(i).source)
            subscripts=getOne(obj.vectorizedOperations,'subscripts',obj.gets(i).source(j));
            instructions=getOne(obj.vectorizedOperations,'instructions',obj.gets(i).source(j));
            var=addConstantVariable(int64(subscripts),'','get');
            % scalar-to-full tensor function type, 1+n inputs, 1 output
            fwrite(fg,[s2ft_funct,1+length(instructions),1,var-1,instructions'-1,newvars(j)-1],intType);  % 0-based;
            nFuncts=nFuncts+1;
        end        

        % io type, # gets, variables
        fwrite(fio,[ioget,length(obj.gets(i).source),newvars],intType);  % 0-based;
        fwrite(fs,[symio,nIOs,length(obj.gets(i).functionName),0],intType);
        fwrite(fs,obj.gets(i).functionName,'char');
        nIOs=nIOs+1;
        nSym=nSym+1;
    end
    for i=1:length(obj.sets)
        % add functions to scalarized full tensors
        newvars=[nVariables+1:nVariables+length(obj.sets(i).destination)];
        nVariables=newvars(end);
        for j=1:length(obj.sets(i).destination)
            subscripts=getOne(obj.vectorizedOperations,'subscripts',obj.sets(i).destination(j));
            instructions=getOne(obj.vectorizedOperations,'instructions',obj.sets(i).destination(j));
            % get subscripts in right order
            [subscripts,k]=sortrows(subscripts',size(subscripts,1):-1:1);
            subscripts=subscripts';
            instructions=instructions(k);
            % full tensor-to-scalar function type, 1 inputs, n outputs
            fwrite(fg,[ft2s_funct,1,length(instructions),newvars(j)-1,instructions'-1],intType);  % 0-based;
            nFuncts=nFuncts+1;
        end        

        % io type, # sets, variables
        fwrite(fio,[ioset,length(obj.sets(i).destination),newvars],intType);  % 0-based;
        fwrite(fs,[symio,nIOs,length(obj.sets(i).functionName),0],intType);
        fwrite(fs,obj.sets(i).functionName,'char');
        nIOs=nIOs+1;
        nSym=nSym+1;
    end
    for i=1:length(obj.copies)
        instructionsDestination=[];
        instructionsSource=[];
        for k=1:length(obj.copies(i).destination)
            subscriptsDestination=getOne(obj.vectorizedOperations,...
                                         'subscripts',obj.copies(i).destination(k));
            iD=getOne(obj.vectorizedOperations,'instructions',obj.copies(i).destination(k));
            instructionsDestination=[instructionsDestination;iD];
            
            subscriptsSource=getOne(obj.vectorizedOperations,'subscripts',obj.copies(i).source(k));
            iS=getOne(obj.vectorizedOperations,'instructions',obj.copies(i).source(k));
            instructionsSource=[instructionsSource;iS];
            
            if ~isequal(subscriptsDestination,subscriptsSource)
                subscriptsDestination,subscriptsSource
                error('saveScalarized: source and destination for Copy command with different sparsity structures\n');
            end
        end
        
        % io type, # copies, source variables, destination variables
        fwrite(fio,[iocopy,length(instructionsSource),...
                    instructionsSource'-1,instructionsDestination'-1],intType);  % 0-based;
        fwrite(fs,[symio,nIOs,length(obj.copies(i).functionName),0],intType);
        fwrite(fs,obj.copies(i).functionName,'char');
        nIOs=nIOs+1;
        nSym=nSym+1;
    end

    %% Save constants (and their symbolic names)
    cfstring=ctypes(strcmp('string',consts));
    cfvint=ctypes(strcmp('int64l_full_vector',consts));
    cftint=ctypes(strcmp('int64l_full_tensor',consts));
    cfdouble=ctypes(strcmp('double_full_tensor',consts));
    for thisC=1:length(constants)
        value=constants(thisC).value;
        sz=size(value);
        nd=length(sz);
        cl=class(value);
        switch cl
          case 'char'
            % constant type, length in byte of data, data
            fwrite(fc,[cfstring,prod(sz)],intType);
            fwrite(fc,value(:)','char');
          case 'int64'
            if size(value,2)==1
                % constant type, length in byte of data, data
                fwrite(fc,[cfvint,8*prod(sz)],intType);
                fwrite(fc,value(:),'int64');
            else
                % constant type, length in byte of data, size, entries
                fwrite(fc,[cftint,8*(1+nd+prod(sz))],intType);
                fwrite(fc,[nd,sz,value(:)'],'int64');
            end
          case 'double'
            % constant type, length in byte of data, # dim, size, data
            fwrite(fc,[cfdouble,8*(1+nd+prod(sz))],intType);
            fwrite(fc,[nd,sz],'int64');
            fwrite(fc,value(:)','double');
          otherwise
            constants(thisC)
            error('csparse/saveScalarized: missing constant class "%s"\n',cl);
        end
        name=constants(thisC).name;
        description=constants(thisC).description;
        if ~isempty(name) || ~isempty(description)
            fwrite(fs,[symcon,thisC-1,length(name),length(description)],intType); % 0-based;
            fwrite(fs,[name,description],'char');
            nSym=nSym+1;
        end
    end

    fclose(fg);
    fclose(fc);
    fclose(fio);
    fclose(fs);

    fprintf('done %d functions, %d ios, %d constants, %d symbols (%.2fms)\n',...
            nFuncts,nIOs,length(constants),nSym,1000*etime(clock,t0));
end
    
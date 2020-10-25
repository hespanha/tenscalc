function saveVectorized(obj,filename)
% saveVectorized(obj,filename)
%
% Saves the csparse object as a computational graph. See computation
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

    fprintf('saveVectorize...\n');

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
    functs=regexprep(functs,'^V','');  % remove leaving "V" to match Tcalculus types
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

    %% open files & write magic
    [~,~,ext]=fileparts(filename);
    if ~isempty(ext)
        error('csparse/saveVectorized: filename "%s" must not include extension\n',filename);
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
    nVariables=int64(height(obj.vectorizedOperations));
    constants=struct('value',{},'name',{},'description',{});

    %% find aliases
    aliases=1:height(obj.vectorizedOperations);
    for thisExp=1:height(obj.vectorizedOperations)
        type=getOne(obj.vectorizedOperations,'type',thisExp);
        operands=getOne(obj.vectorizedOperations,'operands',thisExp)';
        if strcmp(type,'variable') && ~isempty(operands)
            aliases(thisExp)=operands;
        end
    end

    nFuncts=0;
    nSym=0;
    %% Save functions graph (and the symbolic names of function and variable nodes)
    for thisExp=1:height(obj.vectorizedOperations)
        type=getOne(obj.vectorizedOperations,'type',thisExp);
        operands=aliases(getOne(obj.vectorizedOperations,'operands',thisExp));
        name=getOne(obj.vectorizedOperations,'name',thisExp);
        description=getOne(obj.vectorizedOperations,'description',thisExp);

        if ~strcmp(type,'variable')
            k=find(strcmp(type,functs));
            if isempty(k)
                error('csparse/saveVectorized: function "%s" not registered\n',type);
            end
        end
        switch(type)
          case 'variable'
            fwrite(fs,[symvar,thisExp-1,length(name),length(description)],intType); % 0-based
            fwrite(fs,[name,description],'char');
            nSym=nSym+1;
            if ~isempty(operands)
                fprintf('  alias: thisExp=%d, alias to=%d, name=%s\n',thisExp-1,operands-1,name);
            end
            % do nothing

          case 'constant'
            value=getOne(obj.vectorizedOperations,'parameters',thisExp);
            var=addConstant(double(value),name,description);
            % write constant function type, constant index, output variable
            fwrite(fg,[const_funct,var,thisExp-1],intType);  % 0-based;

          case {'zeros','ones','eye'} %% dimension
            osize=getOne(obj.vectorizedOperations,'osize',thisExp);
            var=addConstantVariable(int64(osize),'',[type,' osize']);
            % function type, 1 input (osize), 1 output
            fwrite(fg,[ftypes(k),1,1,var-1,thisExp-1],intType);  % 0-based;

          case {'norm2','norminf','full','diag','transpose','ctranspose','lu','chol'} %% 1 tensor
            % function type, 1 input, 1 output
            fwrite(fg,[ftypes(k),1,1,operands-1,thisExp-1],intType);  % 0-based;

          case {'mldivide_l1','mldivide_u','clp','rdivide'} %% 2 tensors
            % function type, 2 inputs (LU matrix,vector), 1 output
            fwrite(fg,[ftypes(k),2,1,operands-1,thisExp-1],intType);  % 0-based;

          case 'subsref' %% list of indices + tensor
            S=getOne(obj.vectorizedOperations,'parameters',thisExp)';
            if ~strcmp(S.type,'()')
                error('csparse/saveVectorized: unkown subsref type "%s" not registered\n',S.type);
            end
            var=cell(length(S.subs),1);
            for i=1:length(S.subs)
                var{i}=addConstantVariable(S.subs{i},'',[type,' subs'])-1; % 0-based
            end
            % function type, n+1 inputs (index1,index2,...,obj), 1 output
            fwrite(fg,[ftypes(k),1+length(var),1,var{:},operands-1,thisExp-1],intType);  % 0-based;

          case {'all','min','max','sum','repmat'} %% dimension + tensor
            dim=getOne(obj.vectorizedOperations,'parameters',thisExp)';
            var=addConstantVariable(int64(dim),'',[type,' dim']);
            % function type, 2 inputs (dim,obj), 1 output
            fwrite(fg,[ftypes(k),2,1,var-1,operands-1,thisExp-1],intType);  % 0-based;

          case {'reshape'} %% size + tensor
            osize=getOne(obj.vectorizedOperations,'osize',thisExp);
            var=addConstantVariable(int64(osize),'',[type,' osize']);
            % function type, 2 inputs (size,obj), 1 output
            fwrite(fg,[ftypes(k),2,1,var-1,operands-1,thisExp-1],intType);  % 0-based;

          case {'cat'} %% dimension + list of tensors
            dim=getOne(obj.vectorizedOperations,'parameters',thisExp)';
            var=addConstantVariable(int64(dim),'',[type,' var']);
            % function type, 1+n inputs (dim,obj), 1 output
            fwrite(fg,[ftypes(k),1+length(operands),1,var-1,operands-1,thisExp-1],intType);  % 0-based;

          case 'plus' %% signs + list of tensors
            signs=getOne(obj.vectorizedOperations,'parameters',thisExp)';
            var=addConstantVariable(int64(signs),'',[type,' signs']);
            % function type, 1+n inputs (signs, terms), 1 output
            fwrite(fg,[ftypes(k),1+length(signs),1,var-1,operands-1,thisExp-1],intType);  % 0-based;

          case 'tprod' %% tensors and indices
            osize=getOne(obj.vectorizedOperations,'osize',thisExp);
            indices=getOne(obj.vectorizedOperations,'parameters',thisExp)';
            var1=addConstantVariable(int64(osize),'',[type,' osize']);
            var2=addConstantVariable(int64(indices{1}),'',[type,' sumsizes']);
            var=cell(length(operands),1);
            for i=1:length(operands)
                var{i}=addConstantVariable(int64(indices{i+1}),'',[type,' indices'])-1; % 0-based
            end
            % function type, 2+2*n inputs, 1 output, size, sumsizes, tensor1, indices1, ..., 1 output
            fwrite(fg,[ftypes(k),2*length(indices),1,var1-1,var2-1],intType);  % 0-based;
            for i=1:length(operands)
                fwrite(fg,[operands(i)-1,var{i}],intType);
            end
            fwrite(fg,thisExp-1,intType);  % 0-based;

          case 'compose' %% fun + tensor
            fun=getOne(obj.vectorizedOperations,'parameters',thisExp)';
            var=addConstantVariable(fun,'',[type,' fun']);
            % function type, 2 inputs (fun, terms), 1 output
            fwrite(fg,[ftypes(k),2,1,var-1,operands-1,thisExp-1],intType);  % 0-based;

          otherwise
            error('csparse/saveVectorized: missing type "%s"\n',type);
        end

        if ~strcmp(type,'variable')
            fwrite(fs,[symfun,nFuncts,length(name),length(description)],intType); % 0-based
            fwrite(fs,[name,description],'char');
            fwrite(fs,[symvar,thisExp-1,length(name),length(description)],intType); % 0-based
            fwrite(fs,[name,description],'char');
            nFuncts=nFuncts+1;
            nSym=nSym+2;
        end
    end

    %% Save I/O functions (and their symbolic names)
    ioget=iotypes(strcmp('get',ios));
    ioset=iotypes(strcmp('set',ios));
    iocopy=iotypes(strcmp('copy',ios));
    nIOs=0;
    for i=1:length(obj.gets)
        % io type, # gets, variables
        fwrite(fio,[ioget,length(obj.gets(i).source),...
                    obj.gets(i).source'-1],intType);  % 0-based;
        fwrite(fs,[symio,nIOs,length(obj.gets(i).functionName),0],intType);
        fwrite(fs,obj.gets(i).functionName,'char');
        nIOs=nIOs+1;
        nSym=nSym+1;
    end
    for i=1:length(obj.sets)
        % io type, # sets, variables
        fwrite(fio,[ioset,length(obj.sets(i).destination),...
                    obj.sets(i).destination'-1],intType);  % 0-based;
        fwrite(fs,[symio,nIOs,length(obj.sets(i).functionName),0],intType);
        fwrite(fs,obj.sets(i).functionName,'char');
        nIOs=nIOs+1;
        nSym=nSym+1;
    end
    for i=1:length(obj.copies)
        % io type, # copies, source variables, destination variables
        fwrite(fio,[iocopy,length(obj.copies(i).source),...
                    obj.copies(i).source'-1,obj.copies(i).destination'-1],intType);  % 0-based;
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
            error('csparse/saveVectorized: missing constant class "%s"\n',cl);
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

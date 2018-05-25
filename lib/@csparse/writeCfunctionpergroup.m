function writeCfunctionpergroup(obj,codeType,Cfunction,Hfunction,logFile,folder,profiling)
% writeCfunctionpergroup(obj,codeType,Cfunction,Hfunction,logFile,folder,profiling)
%    Write C code to implements all the gets, sets, and copies in
%    a csparse object.
% Inputs:
% obj - csparse object
% codeType - Type of code produced:'
%              ''C'' - all computations done pure C code.'
%                      Impact on non-optimized compilation:'
%                        . medium compilation times'
%                        . largest code size'
%                        . slowest run times'
%                      Impact on optimized code (ideal for final code):'
%                        Gives the most freedom to the compiler for optimization'
%                        . slowest compile optimization times'
%                        + fastest run times'
%                        + smallest code sizes'
%              ''C+asmSB'' - little C code, with most of the computations done'
%                      by small blocks of inlined assembly code'
%                      Impact on non-optimized compilation:'
%                        . medium compilation times'
%                        . medium code size'
%                        . medium run times'
%                      Impact on optimized code:'
%                        Most of the compiler optimization is restricted to re-ordering'
%                        and/or inlining the small blocks of asm code'
%                        . medium compile optimization times,'
%                        . medium run times'
%                        . medium code sizes'
%              ''C+asmLB'' - little C code, with most of the computations done'
%                      by large blocks of inlined assembly code'
%                      Impact on non-optimized compilation (ideal for testing):'
%                        + fastest compilation times'
%                        + smallest code size'
%                        + fastest run times'
%                      Impact on optimized compilation:'
%                        Most of the compiler optimization is restricted to re-ordering'
%                        and/or inlining the large blocks of asm code'
%                        . fastest compile optimization times'
%                        . slowest run times'
%                        . largest optimized code sizes (due to inlining large blocks)'
% Cfunction - file where the C code should be written
% Hfunction - file where the C function headers should be written
% logFile   - file where statistics information should be written
% folder    - folder where all files will be written
% profiling - when non-zero adds profiling code
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

verboseLevel=0;

allFunctionSameFile=true;

m_define=false;  % use of #define seems to lead to larger and
                 % slower code (both for -O0 and -O1)

m_storageclass='register'; % 'register' does not seem to have any
                           % effect on -O0 or -O1 for the apple
                           % compiler, but may help for other compilers

if verboseLevel>0
    fprintf('  writeCfunctionpergroup ');
    t0=clock;
end

[~,basename]=fileparts(Cfunction);

computeFunctions=sprintf('%s_group',basename);
callComputeFunctions=sprintf('  if (groupStatus[%%d]==0) {groupStatus[%%d]=1;%s%%d();}\n',computeFunctions);
NcallComputeFunctions=3;
if profiling
    callComputeFunctions=[callComputeFunctions,'  countCallGroup[%d]++;\n'];
    NcallComputeFunctions=NcallComputeFunctions+1;
    [iTypes,pTypes]=instructionTypes();
end
if obj.debug>0
    callComputeFunctions=[callComputeFunctions,'  else { fprintf(stderr,"(skipping group %d)\\n"); }\n'];
    NcallComputeFunctions=NcallComputeFunctions+1;
end

debugFunction=sprintf('%s_debug',basename);
callDebugFunction=sprintf('  %s();\n',debugFunction);

if ~isempty(Hfunction)
    fih=fopen(fullfile(folder,Hfunction),'w');
    if fih<0
        error('writeCfunctionpergroup: unable to open ''%s''\n',...
              fullfile(folder,Hfunction))
    end
    includeFile(fih,'GPL.c');
    fprintf(fih,'#include <stdint.h> /* needed by uint64_t */\n');
    fprintf(fih,'#include <matrix.h> /* needed by IPM solvers */\n');
end

if ~isempty(logFile)
    fil=fopen(fullfile(folder,logFile),'w');
    if fil<0
        error('writeCfunctionpergroup: unable to open ''%s''\n',...
              fullfile(folder,logFile))
    end
end
fid=fopen(fullfile(folder,Cfunction),'w');
if fid<0
    error('writeCfunctionpergroup: unable to open ''%s''\n',...
          fullfile(folder,Cfunction))
end
includeFile(fid,'GPL.c');
fprintf(fid,'\n/***********************************************/\n');
fprintf(fid,'/* Start of code generated by compile2C/writeC */\n');
fprintf(fid,'/***********************************************/\n\n');
fprintf(fid,'#include <string.h> /* needed by memcpy */\n');
fprintf(fid,'#include <stdint.h> /* needed by uint64_t */\n');
fprintf(fid,'#include <inttypes.h> /* needed by PRId64 */\n');
fprintf(fid,'#include <stdio.h>  /* needed by printf */\n');
fprintf(fid,'#include <fcntl.h>  /* needed by open */\n');
fprintf(fid,'#include <float.h>  /* needed by DBL_MAX */\n');
fprintf(fid,'#include <math.h>   /* needed by fmin, fabs, etc. */\n');
fprintf(fid,'#include <time.h>   /* needed by clock() */\n\n');

fprintf(fid,'#ifdef __APPLE__\n');
fprintf(fid,'#include <unistd.h> /* needed by open */\n');
fprintf(fid,'/* To get nano-sec time in OSX */\n');
fprintf(fid,'#include <mach/mach.h>        /* needed by OSX mach_absolute_time() */\n');
fprintf(fid,'#include <mach/mach_time.h>   /* needed by OSX mach_absolute_time() */\n\n');
fprintf(fid,'#define clock_t uint64_t\n');
fprintf(fid,'#define clock mach_absolute_time\n');
fprintf(fid,'#undef  CLOCKS_PER_SEC\n');
fprintf(fid,'#define CLOCKS_PER_SEC (uint64_t)1000000000\n');
fprintf(fid,'#elif __linux__\n');
fprintf(fid,'#include <unistd.h> /* needed by open */\n');
fprintf(fid,'#define clock_t uint64_t\n');
fprintf(fid,'#define clock clock_ns\n');
fprintf(fid,'#undef  CLOCKS_PER_SEC\n');
fprintf(fid,'#define CLOCKS_PER_SEC (uint64_t)1000000000\n');
fprintf(fid,'clock_t clock_ns() { struct timespec tp;  clock_gettime(CLOCK_REALTIME,&tp); return CLOCKS_PER_SEC*tp.tv_sec+tp.tv_nsec;}\n');
fprintf(fid,'#elif _WIN32\n');
fprintf(fid,'#include <Windows.h>\n');
fprintf(fid,'#include <tchar.h>\n');
fprintf(fid,'#include <stdio.h>\n');
fprintf(fid,'#include <stdint.h>\n');
fprintf(fid,'#define clock_t uint64_t\n');
fprintf(fid,'#define clock clock_ns\n');
fprintf(fid,'#undef  CLOCKS_PER_SEC\n');
fprintf(fid,'#define CLOCKS_PER_SEC (uint64_t)1000000000\n');
fprintf(fid,'clock_t clock_ns() {LARGE_INTEGER li; LARGE_INTEGER freq; QueryPerformanceCounter(&li); QueryPerformanceFrequency(&freq); return (li.QuadPart*CLOCKS_PER_SEC/freq.QuadPart);}\n');
fprintf(fid,'#endif\n');

%groups=getMulti(obj.instructions,'group');
groups=obj.instructionsGroup;

nGroups=max(groups);

if verboseLevel>0
    fprintf('WARNING: writeC: each instruction is being assigned a specific memory address but one should reuse memory\n');
end

%% status of dependency groups
fprintf(fid,'/* Status of the dependency groups */\n');
fprintf(fid,'int groupStatus[]={%s};\n\n',index2str(zeros(1,nGroups)));
if profiling
    fprintf(fid,'/* Profiling */\n');
    fprintf(fid,'long    countCallGroup[]={%s};\n',index2str(zeros(1,nGroups)));
    fprintf(fid,'long    countExecuteGroup[]={%s};\n',index2str(zeros(1,nGroups)));
    fprintf(fid,'clock_t timeExecuteGroup[]={%s};\n',index2str(zeros(1,nGroups)));
    names=sprintf('"%s",',obj.sets(:).functionName);
    fprintf(fid,'char   *setNames[]={%s};\n',names(1:end-1));
    fprintf(fid,'long    countCallSet[]={%s};\n',index2str(zeros(1,length(obj.sets))));
    fprintf(fid,'clock_t timeExecuteSet[]={%s};\n',index2str(zeros(1,length(obj.sets))));
    names=sprintf('"%s",',obj.gets(:).functionName);
    fprintf(fid,'char    *getNames[]={%s};\n',names(1:end-1));
    fprintf(fid,'long    countCallGet[]={%s};\n',index2str(zeros(1,length(obj.gets))));
    fprintf(fid,'clock_t timeExecuteGet[]={%s};\n',index2str(zeros(1,length(obj.gets))));
    names=sprintf('"%s",',obj.copies(:).functionName);
    fprintf(fid,'char    *copyNames[]={%s};\n',names(1:end-1));
    fprintf(fid,'long    countCallCopy[]={%s};\n',index2str(zeros(1,length(obj.copies))));
    fprintf(fid,'clock_t timeExecuteCopy[]={%s};\n',index2str(zeros(1,length(obj.copies))));
    fprintf(fid,'uint64_t countFlops[]={%s};\n',index2str(zeros(1,length(fields(pTypes)))));
    names=fields(pTypes);
    names=sprintf('"%s",',names{:});
    fprintf(fid,'char    *flopsNames[]={%s};\n',names(1:end-1));
    fprintf(fid,'/* Copied from writeCprofiling.c */\n');
    wh=which('csparse');
    f=fopen(fullfile(fileparts(wh),'writeCprofiling.c'),'r');
    s=fread(f);
    fprintf(fid,'%s',s);
    fclose(f);
    fprintf(fid,'/* end copy from writeCprofiling.c */\n\n');
end

%% Storage areas
fprintf(fid,'#define SCRATCHBOOK_TYPE %s\n',obj.scratchbookType);
fprintf(fid,'SCRATCHBOOK_TYPE dbl_max=DBL_MAX;\n');
%% Write declarations for atomic variables
if length(obj.atomicVariables)>0
    fprintf(fid,'#include <umfpack.h>\n');
    fprintf(fid,'double *null = (double *) NULL ;\n');
    % Symbolic & numeric
    fprintf(fid,'void *Symbolic[%d]={NULL',length(obj.atomicVariables));
    fprintf(fid,[repmat(',NULL',1,length(obj.atomicVariables)-1),'};\n']);
    fprintf(fid,'void *Numeric[%d]={NULL',length(obj.atomicVariables));
    fprintf(fid,[repmat(',NULL',1,length(obj.atomicVariables)-1),'};\n']);
    % Ax
    for i=1:length(obj.atomicVariables)
        fprintf(fid,'double Ax%d[%d];\n',i-1,length(obj.atomicVariables(i).Ai));
    end
    fprintf(fid,'double *Ax[%d]={Ax0',length(obj.atomicVariables));
    if length(obj.atomicVariables)>1
        fprintf(fid,',Ax%d',1:length(obj.atomicVariables)-1);
    end
    fprintf(fid,'};\n');
    % Ap
    for i=1:length(obj.atomicVariables)
        fprintf(fid,'long Ap%d[%d]={%d',i-1,length(obj.atomicVariables(i).Ap),obj.atomicVariables(i).Ap(1));
        if length(obj.atomicVariables(i).Ap)>1
            for j=2:length(obj.atomicVariables(i).Ap)
                fprintf(fid,',%d',obj.atomicVariables(i).Ap(j));
                if mod(j,50)==0
                    fprintf(fid,'\n');
                end
            end
        end
        fprintf(fid,'};\n');
    end
    fprintf(fid,'long *Ap[%d]={Ap0',length(obj.atomicVariables));
    if length(obj.atomicVariables)>1
        fprintf(fid,',Ap%d',1:length(obj.atomicVariables)-1);
    end
    fprintf(fid,'};\n');
    % Ai
    for i=1:length(obj.atomicVariables)
        fprintf(fid,'long Ai%d[%d]={%d',i-1,length(obj.atomicVariables(i).Ai),obj.atomicVariables(i).Ai(1));
        if length(obj.atomicVariables(i).Ai)>1
            for j=2:length(obj.atomicVariables(i).Ai)
                fprintf(fid,',%d',obj.atomicVariables(i).Ai(j));
                if mod(j,50)==0
                    fprintf(fid,'\n');
                end
            end
        end
        fprintf(fid,'};\n');
    end
    fprintf(fid,'long *Ai[%d]={Ai0',length(obj.atomicVariables));
    if length(obj.atomicVariables)>1
        fprintf(fid,',Ai%d',1:length(obj.atomicVariables)-1);
    end
    fprintf(fid,'};\n');
    % b
    for i=1:length(obj.atomicVariables)
        fprintf(fid,'double b%d[%d];\n',i-1,obj.atomicVariables(i).n);
    end
    fprintf(fid,'double *b[%d]={b0',length(obj.atomicVariables));
    if length(obj.atomicVariables)>1
        fprintf(fid,',b%d',1:length(obj.atomicVariables)-1);
    end
    fprintf(fid,'};\n');
    % x
    for i=1:length(obj.atomicVariables)
        fprintf(fid,'double x%d[%d];\n',i-1,obj.atomicVariables(i).n);
    end
    fprintf(fid,'double *x[%d]={x0',length(obj.atomicVariables));
    if length(obj.atomicVariables)>1
        fprintf(fid,',x%d',1:length(obj.atomicVariables)-1);
    end
    fprintf(fid,'};\n');
end
%% Dynamic library
fprintf(fid,'#ifdef DYNAMIC_LIBRARY\n');
fprintf(fid,'#ifdef __APPLE__\n');
fprintf(fid,'#define EXPORT __attribute__((visibility("default")))\n');
fprintf(fid,'#elif __linux__\n');
fprintf(fid,'#define EXPORT __attribute__((visibility("default")))\n');
fprintf(fid,'#elif _WIN32\n');
fprintf(fid,'#define EXPORT __declspec(dllexport)\n');
fprintf(fid,'#endif\n');
fprintf(fid,'#include <stdlib.h> /* needed for malloc */\n');
fprintf(fid,'/* Storage area */\n');
fprintf(fid,'SCRATCHBOOK_TYPE *scratchbook=NULL;\n');
fprintf(fid,'/* Auxiliary values */\n');
fprintf(fid,'#ifdef __APPLE__\n');
fprintf(fid,'/* Initializer */\n');
fprintf(fid,'__attribute__((constructor))\n');
fprintf(fid,'static void initializer(void) {\n');
fprintf(fid,'  scratchbook=malloc(sizeof(*scratchbook)*(size_t)%d);\n',max(obj.memoryLocations));
fprintf(fid,'  printf("%%s: loaded dynamic library & allocated memory for scrapbook (%%"PRIXPTR")\\n", __FILE__,(uintptr_t)scratchbook);\n');
fprintf(fid,'}\n');
fprintf(fid,'/* Finalizer */\n');
fprintf(fid,'__attribute__((destructor))\n');
fprintf(fid,'static void finalizer(void) {\n');
for i=1:length(obj.atomicVariables)
    %% Somehow this memory is being freed before getting here.
    % fprintf(fid,'  if (Numeric[%d]) umfpack_di_free_numeric(&Numeric[%d]);\n',i-1,i-1);
end
if profiling
    fprintf(fid,'  profilingView("%s.profile");\n',fullfile(folder,Cfunction));
end
fprintf(fid,'  free(scratchbook);\n');
fprintf(fid,'  printf("%%s: freed scrapbook, unloading dynamic library\\n", __FILE__);\n');
fprintf(fid,'}\n');
fprintf(fid,'#elif __linux__\n');
fprintf(fid,'/* Initializer */\n');
fprintf(fid,'__attribute__((constructor))\n');
fprintf(fid,'static void initializer(void) {\n');
fprintf(fid,'  scratchbook=malloc(sizeof(*scratchbook)*(size_t)%d);\n',max(obj.memoryLocations));
fprintf(fid,'  printf("%%s: loaded dynamic library & allocated memory for scrapbook (%%"PRIXPTR")\\n", __FILE__,(uintptr_t)scratchbook);\n');
fprintf(fid,'}\n');
fprintf(fid,'/* Finalizer */\n');
fprintf(fid,'__attribute__((destructor))\n');
fprintf(fid,'static void finalizer(void) {\n');
for i=1:length(obj.atomicVariables)
    %% Somehow this memory is being freed before getting here.
    % fprintf(fid,'  if (Numeric[%d]) umfpack_di_free_numeric(&Numeric[%d]);\n',i-1,i-1);
end
if profiling
    fprintf(fid,'  profilingView("%s.profile");\n',fullfile(folder,Cfunction));
end
fprintf(fid,'  free(scratchbook);\n');
fprintf(fid,'  printf("%%s: freed scrapbook, unloading dynamic library\\n", __FILE__);\n');
fprintf(fid,'}\n');
fprintf(fid,'#elif _WIN32\n');
fprintf(fid,'#include <windows.h>\n');
fprintf(fid,'BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)\n');
fprintf(fid,'{\n');
fprintf(fid,'    if (fdwReason == DLL_PROCESS_ATTACH) {\n');
fprintf(fid,'       scratchbook=malloc(sizeof(*scratchbook)*(size_t)%d);\n',max(obj.memoryLocations));
fprintf(fid,'       printf("%%s: loaded dynamic library & allocated memory for scrapbook (%%"PRIXPTR")\\n", __FILE__,(uintptr_t)scratchbook);\n');
fprintf(fid,'       return TRUE; }\n');
fprintf(fid,'    else if (fdwReason == DLL_PROCESS_DETACH) {\n');
for i=1:length(obj.atomicVariables)
    %% Somehow this memory is being freed before getting here.
    % fprintf(fid,'       if (Numeric[%d]) umfpack_di_free_numeric(&Numeric[%d]);\n',i-1,i-1);
end
if profiling
    fprintf(fid,'  profilingView("%s.profile");\n',fullfile(folder,Cfunction));
end
fprintf(fid,'       free(scratchbook);\n');
fprintf(fid,'       printf("%%s: freed scrapbook, unloading dynamic library\\n", __FILE__);\n');
fprintf(fid,'       return TRUE; }\n');
fprintf(fid,'}\n');
fprintf(fid,'#endif\n');

fprintf(fid,'#else\n');
%% Static library
fprintf(fid,'#define EXPORT \n');
fprintf(fid,'SCRATCHBOOK_TYPE scratchbook[%d];\n',max(obj.memoryLocations));
fprintf(fid,'#endif\n');

%% Create debug function to print memory in the screen
if obj.debug>1
    fprintf(fid,'EXPORT void %s(void) {\n',debugFunction);
    if m_define
        fprintf(fid,'#define m scratchbook\n');
    else
        fprintf(fid,'  %s SCRATCHBOOK_TYPE *m=scratchbook;\n',m_storageclass);
    end
    for i=1:height(obj.vectorizedOperations)
        name=getOne(obj.vectorizedOperations,'name',i);
        osize=getOne(obj.vectorizedOperations,'osize',i);
        subscripts=getOne(obj.vectorizedOperations,'subscripts',i);
        instructions=getOne(obj.vectorizedOperations,'instructions',i);
        fprintf(fid,'  fprintf(stderr,"  %-20s [%-5s]:");\n',name,index2str(osize));
        for j=1:size(subscripts,2)
            fprintf(fid,'  fprintf(stderr," [%-7s] = %%10.3lg,",m[%d]);\n',...
                    index2str(subscripts(:,j))',obj.memoryLocations(instructions(j))-1);
        end
        fprintf(fid,'  fprintf(stderr,"\\n");\n');
        
    end
    fprintf(fid,'}\n');
end

if ~isempty(logFile)

    if 1
        %% Write description of the csparse object 
        %% very time consuming for obj.str(1)
        fprintf(fil,'* Sparse object:\n');
        fprintf(fil,'%s',obj.str(0));
    end
    
    %% Write Dependency groups to the log file
    fprintf(fil,'\n* Dependency groups:\n');
    fprintf(fil,'%35s','group');
    fprintf(fil,'%5d',0:nGroups-1);
    fprintf(fil,'\n');
    fprintf(fil,'%s\n',repmat('-',1,35));
    for i=1:size(obj.dependencyGroups,2)
        fprintf(fil,'%-35s',obj.dependencyGroupColName{i});
        fprintf(fil,'%5d',full(obj.dependencyGroups(:,i)));
        fprintf(fil,'\n');
    end
    fprintf(fil,'\n');
end

if verboseLevel>0
    fprintf('writeC: creating computation functions\n');
end

%% Create functions to do the computations
for i=1:nGroups

    if allFunctionSameFile
        fig=fid;
    else
        name=sprintf('%s_GROUP%d.c',regexprep(Cfunction,'\.c',''),i);
        fprintf(fid,'extern void %s%d();\n',computeFunctions,i-1);
        fig=fopen(fullfile(folder,name),'w');
        fprintf(fig,'#include <math.h>   /* needed by fmin, fabs, etc. */\n');
        fprintf(fid,'extern SCRATCHBOOK_TYPE scratchbook[%d];\n',max(obj.memoryLocations));
        fprintf(fig,'#ifdef DYNAMIC_LIBRARY\n');
        fprintf(fig,'#define EXPORT __attribute__((visibility("default")))\n');
        fprintf(fig,'#else\n');
        fprintf(fig,'#define EXPORT \n');
        fprintf(fig,'#endif\n');
    end
    
    if ~isempty(Hfunction)
        fprintf(fih,'extern void %s%d();\n',computeFunctions,i-1);
    end
    fprintf(fig,'EXPORT void %s%d() {\n',computeFunctions,i-1);
    k=find(groups==i);
    if obj.debug>0
        fprintf(fig,'  fprintf(stderr,"computing group %d\\n");\n',i-1);
    end
    if profiling
        fprintf(fig,'  countExecuteGroup[%d]++;\n',i-1);
        fprintf(fig,'  clock_t t0=clock();\n',i-1);
    end
    switch (codeType)
      case 'C',
        % Most flexibility for the compiler to optimize
        if m_define
            fprintf(fid,'#define m scratchbook\n');
        else
            fprintf(fid,'  %s SCRATCHBOOK_TYPE *m=scratchbook;\n',m_storageclass);
        end
        if 0
            writeCinstructions(obj,fig,k);    
        else
            % not pretty, but apparently C cannot use matlab's fig
            countFlops=writeCinstructionsC(int64(k),int64(obj.memoryLocations'));
            f=fopen('tmp_toremove.c','r');
            str=fread(f,inf);
            fclose(f);
            delete('tmp_toremove.c');
            fwrite(fig,str);

            for q=1:length(countFlops)
                if countFlops(q)>0
                    if profiling
                        fprintf(fig,'   countFlops[%d] += %d;\n',q-1,countFlops(q));
                    else
                        fprintf(fig,'   // countFlops[%d] += %d;\n',q-1,countFlops(q));
                    end
                end
            end
        end
      case 'C+asmSB',
        if ~strcmp(obj.scratchbookType,'double')
            error('writeC: scratchbookType %s not implemented for codeType %s\n',scratchbookType,codeType);
        end
        writeCasmSBinstructions(obj,fig,k); % Compiler can mostly optimize by inlining (small blocks)
      case 'C+asmLB',
        if ~strcmp(obj.scratchbookType,'double')
            error('writeC: scratchbookType %s not implemented for codeType %s\n',scratchbookType,codeType);
        end
        if 0
            writeCasmLBinstructions(obj,fig,k);  % Compiler can mostly optimize by inlining (large blocks)
        else
            % not pretty, but apparently C cannot use matlab's fig
            writeAsmInstructionsC(int64(k),int64(obj.memoryLocations'));  
            f=fopen('tmp_toremove.c','r');
            str=fread(f,inf);
            fclose(f);
            delete('tmp_toremove.c');
            fwrite(fig,str);
        end
      otherwise
        error('writeC: unknown value of codeType (''%s'')\n',codeType);
    end
    if profiling
        fprintf(fig,'  timeExecuteGroup[%d] += clock()-t0;\n',i-1);
    end
    if obj.debug>1
        fprintf(fig,'  %s\n',callDebugFunction);
    end
    fprintf(fig,'}\n');

    %% one function per file?
    if ~allFunctionSameFile
        %% one function per file 
        fclose(fig);
    end

end

if verboseLevel>0
    fprintf('writeC: creating methods functions\n');
end

%% Create functions for sets
if verboseLevel>0
    fprintf('writeC: creating sets functions\n');
end
for i=1:length(obj.sets)
    if verboseLevel>0
        fprintf('   void %s(double *input);\n',obj.sets(i).functionName);
    end
    if ~isempty(Hfunction)
        fprintf(fih,'extern void %s(double *);\n',obj.sets(i).functionName);
    end
    fprintf(fid,'EXPORT void %s(double *input) {\n',obj.sets(i).functionName);
    if profiling
        fprintf(fid,'  countCallSet[%d]++;\n',i-1);
        fprintf(fid,'  clock_t t0=clock();\n',i-1);
    end
    if m_define
        fprintf(fid,'#define m scratchbook\n');
    else
        fprintf(fid,'  %s SCRATCHBOOK_TYPE *m=scratchbook;\n',m_storageclass);
    end
    if obj.debug>1
        fprintf(fid,'        %s\n',callDebugFunction);
    end
    if obj.debug>0
        fprintf(fid,'  fprintf(stderr,"running %s\\n");\n',obj.sets(i).functionName);
    
    end
    
    subscripts=getOne(obj.vectorizedOperations,'subscripts',obj.sets(i).destination);
    instructions=getOne(obj.vectorizedOperations,'instructions',obj.sets(i).destination);
    instructions=obj.memoryLocations(instructions);
    osize=getOne(obj.vectorizedOperations,'osize',obj.sets(i).destination);
    
    if size(subscripts,2)~=prod(osize)
        subscripts,size(subscripts,2),prod(osize)
        error('setting non-full variable %s\n',obj.gets(i).functionName);
    end

    % get subscripts in right order
    [subscripts,k]=sortrows(subscripts',size(subscripts,1):-1:1);
    subscripts=subscripts';
    instructions=instructions(k);
    
    % copy values
    if all(instructions(2:end)-instructions(1:end-1)==1) && strcmp(obj.scratchbookType,'double')
        % consecutive source
        fprintf(fid,'  memcpy(m+%d,input,(size_t)%d*sizeof(double));\n',...
                instructions(1)-1,length(instructions));
    else
        fprintf(fid,'  m[%d]=*(input++);\n',instructions-1);
    end
    % set group status
    if isempty(obj.sets(i).childrenGroups)
        %fprintf('writeC: %s has no effect\n',obj.sets(i).functionName);
    else
        if ~isempty(obj.sets(i).childrenGroups)
            fprintf(fid,'  groupStatus[%d]=0;\n',obj.sets(i).childrenGroups-1);
        end
    end
    if profiling
        fprintf(fid,'  timeExecuteSet[%d] += clock()-t0;\n',i-1);
    end
    fprintf(fid,'}\n');
end

%% Create functions for gets
if verboseLevel>0
    fprintf('writeC: creating gets functions\n');
end
for i=1:length(obj.gets)
    nSources=length(obj.gets(i).source);
    % write function header
    fprintf(fid,'EXPORT void %s(',obj.gets(i).functionName);
    if ~isempty(Hfunction)
        fprintf(fih,'extern void %s(',obj.gets(i).functionName);
    end
    for j=1:nSources
        if length(osize)==2 & size(subscripts,2)~=prod(osize)
            % sparse
            if ~isempty(Hfunction)
                fprintf(fih,'int64_t **y%d_p,int64_t **y%d_i,double **y_v%d',j,j,j);
            end
            fprintf(fid,'double *y%d',j);
        else
            % full
            if ~isempty(Hfunction)
                fprintf(fih,'double *y%d',j);
            end
            fprintf(fid,'double *y%d',j);
        end
        if j<nSources
            if ~isempty(Hfunction)
                fprintf(fih,',');
            end
            fprintf(fid,',');
        end
    end
    if ~isempty(Hfunction)
        fprintf(fih,');\n');
    end
    fprintf(fid,') {\n');
    % write function body
    if profiling
        fprintf(fid,'  countCallGet[%d]++;\n',i-1);
        fprintf(fid,'  clock_t t0=clock();\n',i-1);
    end
    if m_define
        fprintf(fid,'#define m scratchbook\n');
    else
        fprintf(fid,'  %s SCRATCHBOOK_TYPE *m=scratchbook;\n',m_storageclass);
    end
    if obj.debug>1
        fprintf(fid,'        %s\n',callDebugFunction);
    end
    if obj.debug>0
        fprintf(fid,'  fprintf(stderr,"running %s\\n");\n',obj.gets(i).functionName);
    end
    
    % perform needed computations
    if ~isempty(obj.gets(i).parentGroups)
        fprintf(fid,callComputeFunctions,...
                repmat(obj.gets(i).parentGroups(:)'-1,NcallComputeFunctions,1));
    end
    
    for j=1:nSources
        atomic=getOne(obj.vectorizedOperations,'atomic',obj.gets(i).source(j));
        if (atomic)
            error('get not supported for atomic operator');
        end
        subscripts=getOne(obj.vectorizedOperations,'subscripts',obj.gets(i).source(j));
        instructions=getOne(obj.vectorizedOperations,'instructions',obj.gets(i).source(j));
        instructions=obj.memoryLocations(instructions)';
        osize=getOne(obj.vectorizedOperations,'osize',obj.gets(i).source(j));
        
        if length(osize)==2 & size(subscripts,2)~=prod(osize)
            %% sparse matrices returned in compressed-column form

            error('return sparse matrices is not complete, make matrix full');
            
            % compute sparse arrays
            %[subscripts',instructions]
            [ji,k]=sortrows(subscripts(end:-1:1,:)');
            % http://www.mathworks.com/help/matlab/apiref/mxsetir.html
            ir=ji(:,2)-1; 
            % http://www.mathworks.com/help/matlab/apiref/mxsetjc.html
            jc=zeros(osize(2)+1,1);
            this=0;
            for l=1:osize(2)
                this=this+sum(ji(:,1)==l);
                jc(l+1)=this;
            end
            %jc
            instrX=instrX(k);
            
            pr=instructions(k)-1; % 0-based indexing

            %% INCOMPLETE
            
            fprintf(fid,'  int64_t p[%d]={0,',osize(2)+1);
            for l=1:osize(2)
                fprintf(fid,'%d,',min(find(ij(:,1)>=l))-1); %0-based indexing
            end
            fprintf(fid,'%d};\n  (*y%d_p)=&p;\n',osize(2),j);
            fprintf(fid,'  int64_t i[%d]={%s};\n',index2str(ir));
            fprintf(fid,'%d};\n  (*y%d_i)=&i;\n',ij(end,2),j);
            fprintf(fid,'  if (*y%d_v==0) y%d_v=malloc((size_t)sizeof(double)*%d);\n',j,j,length(instructions));
            if all(instructions(2:end)-instructions(1:end-1)==1) && strcmp(obj.scratchbookType,'double')
                % consecutive destination
                fprintf(fid,'  memcpy(*y%d_v,m+%d,(size_t)%d*sizeof(double));\n',...
                        j,instructions(1)-1,length(instructions));
            else
                fprintf(fid,'{ double *y=*y%d_v;\n',j);
                fprintf(fid,'  *(y++)=m[%d];\n',instructions-1);
                fprintf(fid,'  }\n');
            end

        else
            % for # dimensions other than 2, must be full
            if size(subscripts,2)~=prod(osize)
                fprintf('\nNonzero entries %d/%d:',size(subscripts,2),prod(osize));
                disp(subscripts);
                if strcmp(obj.gets(i).functionName,'getGrad__')
                    % error specific to optimization
                    error('some optimization variables do not affect the criteria nor the constraints (look for zeros in gradient) ');
                else
                    error('getting non-full variable %s(%d)\n',obj.gets(i).functionName,j);
                end
            end

            if ~isempty(instructions)
                % get subscripts in right order
                [subscripts,k]=sortrows(subscripts',size(subscripts,1):-1:1);
                subscripts=subscripts';
                instructions=instructions(k);
        
                if all(instructions(2:end)-instructions(1:end-1)==1) &&...
                        strcmp(obj.scratchbookType,'double')
                    % consecutive destination
                    fprintf(fid,'  memcpy(y%d,m+%d,(size_t)%d*sizeof(double));\n',...
                            j,instructions(1)-1,length(instructions));
                else
                    fprintf(fid,'  *(y%d++)=m[%d];\n',...
                            [repmat(j,1,length(instructions));instructions'-1]);
                end
            end
        end
    end
    if profiling
        fprintf(fid,'  timeExecuteGet[%d] += clock()-t0;\n',i-1);
    end
    fprintf(fid,'}\n');
end

%% Create functions for copies
if verboseLevel>0
    fprintf('writeC: creating copies functions\n');
end
for i=1:length(obj.copies)
    if verboseLevel>0
        fprintf('   void %s();\n',obj.copies(i).functionName);
    end
    if ~isempty(Hfunction)
        fprintf(fih,'extern void %s();\n',obj.copies(i).functionName);
    end
    fprintf(fid,'EXPORT void %s() {\n',obj.copies(i).functionName);
    if profiling
        fprintf(fid,'  countCallCopy[%d]++;\n',i-1);
        fprintf(fid,'  clock_t t0=clock();\n',i-1);
    end
    if m_define
        fprintf(fid,'#define m scratchbook\n');
    else
        fprintf(fid,'  %s SCRATCHBOOK_TYPE *m=scratchbook;\n',m_storageclass);
    end
    if obj.debug>1
        fprintf(fid,'        %s\n',callDebugFunction);
    end
    if obj.debug>0
        fprintf(fid,'  fprintf(stderr,"running %s\\n");\n',obj.copies(i).functionName);
    end
    
    % perform needed computations
    if ~isempty(obj.copies(i).parentGroups)
        fprintf(fid,callComputeFunctions,...
                repmat(obj.copies(i).parentGroups(:)'-1,NcallComputeFunctions,1));
    end
    
    % copy commands
    instructionsDestination=[];
    instructionsSource=[];
    for k=1:length(obj.copies(i).destination)
        atomic=getOne(obj.vectorizedOperations,'atomic',obj.copies(i).destination(k));
        if (atomic)
            error('copies not supported for atomic operator');
        end
        subscriptsDestination=getOne(obj.vectorizedOperations,...
                                     'subscripts',obj.copies(i).destination(k));
        iD=getOne(obj.vectorizedOperations,'instructions',obj.copies(i).destination(k));
        instructionsDestination=[instructionsDestination;iD];
        
        atomic=getOne(obj.vectorizedOperations,'atomic',obj.copies(i).source(k));
        if (atomic)
            error('copies not supported for atomic operator');
        end
        subscriptsSource=getOne(obj.vectorizedOperations,'subscripts',obj.copies(i).source(k));
        iS=getOne(obj.vectorizedOperations,'instructions',obj.copies(i).source(k));
        instructionsSource=[instructionsSource;iS];
        
        if ~isequal(subscriptsDestination,subscriptsSource)
            subscriptsDestination,subscriptsSource
            error('writeC: source and destination for Copy command with different sparsity structures\n');
        end
    end
    instructionsDestination=obj.memoryLocations(instructionsDestination);
    instructionsSource=obj.memoryLocations(instructionsSource);
    
    if length(instructionsDestination)>1 &&...
            all(instructionsDestination(2:end)-instructionsDestination(1:end-1)==1) && ...
            all(instructionsSource(2:end)-instructionsSource(1:end-1)==1)  && ...
            strcmp(obj.scratchbookType,'double')
        % consecutive source and destination
        fprintf(fid,'  memcpy(m+%d,m+%d,(size_t)%d*sizeof(*m));\n',...
                instructionsDestination(1)-1,instructionsSource(1)-1,...
                length(instructionsDestination));
    elseif ~isempty(instructionsDestination)
        fprintf(fid,'  m[%d]=m[%d];\n',[instructionsDestination(:)';instructionsSource(:)']-1);
    end

    % set group status
    if ~isempty(obj.copies(i).childrenGroups)
        fprintf(fid,'  groupStatus[%d]=0;\n',obj.copies(i).childrenGroups-1);
    end

    if profiling
        fprintf(fid,'  timeExecuteCopy[%d] += clock()-t0;\n',i-1);
    end
    fprintf(fid,'}\n');
end

%% Create functions for saves
if verboseLevel>0
    fprintf('writeC: creating saves functions\n');
end
for i=1:length(obj.saves)
    % create function
    if verboseLevel>0
        fprintf('   void %s(magic=%d);\n',obj.saves(i).functionName,obj.saves(i).magic);
    end
    if ~isempty(Hfunction)
        fprintf(fih,'extern void %s(char *filename);\n',obj.saves(i).functionName);
    end
    fprintf(fid,'EXPORT void %s(char *filename) {\n',obj.saves(i).functionName);
    if m_define
        fprintf(fid,'#define m scratchbook\n');
    else
        fprintf(fid,'  %s SCRATCHBOOK_TYPE *m=scratchbook;\n',m_storageclass);
    end
    fprintf(fid,'#if _WIN32\n');
    fprintf(fid,'  int fd=open(filename,O_WRONLY|O_CREAT|O_TRUNC);\n');
    fprintf(fid,'#else\n');
    fprintf(fid,'  int fd=open(filename,O_WRONLY|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR);\n');
    fprintf(fid,'#endif\n');
    fprintf(fid,'  int64_t magic=%d;\n',obj.saves(i).magic);
    fprintf(fid,'  printf("Saving \\"%%s\\" (fd=%%d,magic=%%"PRId64")\\n",filename,fd,magic);\n');
    fprintf(fid,'  write(fd,&magic,sizeof(int64_t));\n');

    atomic=getOne(obj.vectorizedOperations,'atomic',obj.saves(i).source);
    if (atomic)
        error('saves not supported for atomic operator');
    end
    subscripts=getOne(obj.vectorizedOperations,'subscripts',obj.saves(i).source);
    instructions=getOne(obj.vectorizedOperations,'instructions',obj.saves(i).source);
    instructions=obj.memoryLocations(instructions);
    if isempty(instructions)
        fprintf('\nwriteCfunctionpergroup: saving ZERO matrix "%s"\n',obj.saves(i).functionName);
    else
        if all(instructions(2:end)-instructions(1:end-1)==1) && strcmp(obj.scratchbookType,'double')
            % consecutive destination
            fprintf(fid,'  write(fd,m+%d,(size_t)%d*sizeof(double));\n',...
                    instructions(1)-1,length(instructions));
        else
            fprintf(fid,'  write(fd,m+%d,sizeof(double));\n',instructions-1);
        end
    end
    fprintf(fid,'  close(fd);\n}\n;');
    % save subscripts
    osize=getOne(obj.vectorizedOperations,'osize',obj.saves(i).source);
    subscripts=getOne(obj.vectorizedOperations,'subscripts',obj.saves(i).source);
    fis=fopen(fullfile(folder,obj.saves(i).filename),'w');
    % magic=int64(intmax('int64')*rand(1)); 
    fprintf('Saving "%s" (fd=%d,magic=%d)\n',obj.saves(i).filename,fis,obj.saves(i).magic);
    fwrite(fis,obj.saves(i).magic,'int64');
    fwrite(fis,int32(length(osize)),'int32');
    fwrite(fis,int64(osize),'int64');
    fwrite(fis,int64(size(subscripts,2)),'int64');
    fwrite(fis,int64(subscripts),'int64');
    fclose(fis);
end

%% Create functions for externalFunctions
if verboseLevel>0
    fprintf('writeC: creating external functions\n');
end
for i=1:length(obj.externalFunctions)
    % create function header
    if verboseLevel>0
        fprintf('   void %s(); // "%s"\n',obj.externalFunctions(i).functionName,...
                obj.externalFunctions(i).fileName);
    end
    if ~isempty(Hfunction)
        fprintf(fih,'extern void %s(\n',obj.externalFunctions(i).functionName);
        if length(obj.externalFunctions(i).inputs)>0
            fprintf(fih,'  /* inputs */');
        end
        sep='';
        for j=1:length(obj.externalFunctions(i).inputs)
            fprintf(fih,'%c\n  %s *%s',sep,matlab2Ctype(obj.externalFunctions(i).inputs(j).type),...
                    obj.externalFunctions(i).inputs(j).name);
            sep=',';
        end
        if length(obj.externalFunctions(i).outputs)>0
            fprintf(fih,'%c\n  /* outputs */',sep);
            sep='';
        end
        for j=1:length(obj.externalFunctions(i).outputs)
            fprintf(fih,'%c\n  %s *%s',sep,matlab2Ctype(obj.externalFunctions(i).outputs(j).type),...
                    obj.externalFunctions(i).outputs(j).name);
            sep=',';
        end
        fprintf(fih,');\n');
    end
    fprintf(fid,'\n/********************************/\n');
    fprintf(fid,'/* defines for function void %s() */\n',obj.externalFunctions(i).functionName);
    fprintf(fid,'/**********************************/\n\n');

    % defines
    if ~isempty(obj.externalFunctions(i).defines)
        names=fields(obj.externalFunctions(i).defines);
        for j=1:length(names)
            value=getfield(obj.externalFunctions(i).defines,names{j});
            if ischar(value)
                fprintf(fid,'#define %s %s\n',names{j},value);
            elseif isnumeric(value) && length(value)==1
                fprintf(fid,'#define %s %g\n',names{j},value);
            else
                value
                error('define %s has invalid value\n',names{j});
            end
        end
    end

    fprintf(fid,'\n/**********************************************************/\n');
    fprintf(fid,'/* source code for function void %s(), copied from \"%s\" */\n',...
            obj.externalFunctions(i).functionName,obj.externalFunctions(i).fileName);
    fprintf(fid,'/**********************************************************/\n\n');

    % append function body
    fic=fopen(obj.externalFunctions(i).fileName,'r');
    if fid<0
        error('writeCfunctionpergroup: unable to open file ''%s'' with external function ''%s'' (%d)\n',...
              obj.externalFunctions(i).fileName,obj.externalFunctions(i).functionName,i);
end
    while 1
        tline=fgetl(fic);
        if ~ischar(tline)
            break
        end
        fprintf(fid,'%s\n',tline);
    end
    fclose(fic);
    fprintf(fid,'\n');
end

%% Close

fprintf(fid,'\n/******************************************/\n');
fprintf(fid,'/*  End of code generated by compile2C.m  */\n');
fprintf(fid,'/******************************************/\n');
fclose(fid);
if ~isempty(Hfunction)
    fclose(fih);
end
if ~isempty(logFile)
    fclose(fil);
end

if verboseLevel>0
    fprintf('done writeCfunctionpergroup (%.3f sec)\n',etime(clock(),t0));
end

end

%% Convert a matlab type to a c type
function str=matlab2Ctype(str)

    switch (str)
      case {'uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
        str=sprintf('%s_t',str);
      case {'double','single'}
      otherwise
        error('createGateway: unkown type %s\n',str);
    end
end

function includeFile(fid,filename)
% includeFile(fid,filename)
%   includes the file named 'filename' into the file with the given
%   descriptor
    
    fii=fopen(filename,'r');
    if fii<0
        error('unable to open file ''%s''\n',filename);
    end
    fprintf(fid,'/* START OF #included "%s" */\n',filename);
    a=fread(fii,inf);
    fwrite(fid,a);
    fprintf(fid,'/* END OF #included "%s" */\n',filename);
    fclose(fii);
    
end


compileInstructionsTable


%% unload library if it is loaded
if libisloaded('instructionsTable')
    unloadlibrary('instructionsTable')
end
[notfound,warnings]=loadlibrary('instructionsTable','instructionsTable.h');
m=libfunctions('instructionsTable');


s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

N=10000;
nP=30;
maxP=1;
nO=10;
maxO=1;

profile on

initInstructionsTable()
fprintf('Calling appendInstruction %d times...');t0=clock;
Itype=int32(23);
for i=1:N
    k=ceil(nP*rand());
    parameters=double(ceil(maxP*rand(1,k)));
    k=ceil(nO*rand());
    operands=int64(maxO*rand(1,k));
    index=appendInstruction(Itype,parameters,operands);
    [Rtype,Rparameters,Roperands]=getInstruction(index);
    if (Itype~=Rtype || ~isequal(parameters,Rparameters) || ~isequal(operands,Roperands))
        Itype
        Rtype
        parameters
        Rparameters
        operands
        Roperands
        error('mismatch\n');
    end
end
fprintf('table size = %d, done (%.3f sec, %.3f us/call)\n',index,etime(clock,t0),1e6*etime(clock,t0)/N);

initInstructionsTable()
fprintf('Calling appendUniqueInstruction %d times...');t0=clock;
for i=1:N
    k=ceil(nP*rand());
    parameters=double(ceil(maxP*rand(1,k)));
    k=ceil(nO*rand());
    operands=int64(maxO*rand(1,k));
    index=appendUniqueInstruction(Itype,parameters,operands);
    [Rtype,Rparameters,Roperands]=getInstruction(index);
    if (Itype~=Rtype || ~isequal(parameters,Rparameters) || ~isequal(operands,Roperands))
        Itype
        Rtype
        parameters
        Rparameters
        operands
        Roperands
        error('mismatch\n');
    end
end
fprintf('table size = %d, done (%.3f sec, %.3f us/call)\n',index,etime(clock,t0),1e6*etime(clock,t0)/N);

profile off
profile viewer


unloadlibrary('instructionsTable')

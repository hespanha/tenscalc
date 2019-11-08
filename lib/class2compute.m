function varargout=class2compute(varargin)
    
    %% Function global help
    declareParameter(...
        'Help', {
            'Creates a set of cmex C functions for performing a csparse computation.'
            ' '
            'The computation is performed through a matlab class with methods'
            'for the set, get, and copy operations.'
                });
    
    localVariables_=parameters4compute(localVariables_);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Retrieve parameters and inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [stopNow,params]=setParameters(nargout,varargin);
    if stopNow
        return 
    end
    

    %% transfer any folder in classname into folder
    [folder,classname]=fileparts(fsfullfile(folder,classname));
    
    %% create folder if it does not exist
    if ~strcmp(folder,'.') && ~exist(folder,'dir')
        fprintf('class2compute: outputs folder ''%s'' does not exist, creating it\n',folder);
        if mkdir(folder)==0
            error('Unable to create folder ''%s''\n',folder)
        end
    end
    
    rmpath(folder);
    addpath(folder);
    
    %% Fix class when gotten from pedigree
    classname=regexprep(classname,'+TS=','_TS_');
    classname=regexprep(classname,'-','_');
    classname=regexprep(classname,'+','_');

    classHelp=helpFromTemplate(classname,csparseObject.template);

    fprintf(' creating matlab code... ');
    t_compile2matlab=clock();
    compile2matlab(csparseObject,...
                   fsfullfile(folder,sprintf('%s.m',classname)),...
                   fsfullfile(folder,sprintf('%s.log',classname)),...
                   classHelp,profiling);
    statistics.time.compile2matlab=etime(clock,t_compile2matlab);
    
    fprintf(' done creating matlab code (%.3f sec)\n',etime(clock,t_compile2matlab));
    
    rehash path;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    varargout=setOutputs(nargout,params);
end
        

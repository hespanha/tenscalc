fprintf('Seeting up path...');
home=[fileparts(which('install_tenscalc')),'/lib'];
folders={home;[home,'/csparse']};

s=path;
old=regexp(s,'[^:]*tenscalc[^:]*','match');
if ~isempty(old)
    fprintf('removing from path:\n');
    disp(old')
    rmpath(old{:})
end

fprintf('adding to path:\n');
addpath(folders{:});
disp(folders)

fprintf('Compiling...\n');
compileInstructionsTable;

if true %strcmp(lower(computer(,'glnxa64')))
    fprintf('ATTENTION:\n');
    fprintf('If using the bash shell, add the following line to your .bashrc file:\n');
    fprintf('  export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n',fileparts(which('instructionsTable_load')));
    fprintf('If using the csh/tcsh shell, add the following line to your .cshrc file:\n');
    fprintf('  setenv LD_LIBRARY_PATH %s:$LD_LIBRARY_PATH\n',fileparts(which('instructionsTable_load')));
end


fprintf('saving path...');
try
    savepath;
catch me
    fprintf('ATTENTION: unable to save path, add following strings to the matlab path:');
    disp(folders)
    rethrow
end

fprintf('done!\n');
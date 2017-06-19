
fprintf('Seeting up path...');
home=[fileparts(which('install_tenscalc')),'/lib'];

fprintf('removing old...');
s = warning('OFF', 'MATLAB:rmpath:DirNotFound');
rmpath(home);
rmpath([home,'/csparse']);
warning(s);

fprintf('adding new...');
addpath(home);
addpath([home,'/csparse']);

fprintf('saving path...');
try
    savepath;
catch me
    fprintf('ATTENTION: unable to save path, add following string to the matlab path:\n%s\n',home);
    rethrow
end

fprintf('compiling cmex functions...');
compileInstructionsTable;

fprintf('done!\n');
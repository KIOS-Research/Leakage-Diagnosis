function load_paths()
%Load neccesary paths

%------------- BEGIN CODE --------------
% addpath(genpath(pwd));
if strcmp(computer('arch'),'win64')
    rmpath(fileparts(which('gurobi.mexw32')))
else
    rmpath(fileparts(which('gurobi.mexw64')))
end
run 'gurobi_setup.m'
disp('Toolkits Loaded.');    
%------------- END OF CODE --------------



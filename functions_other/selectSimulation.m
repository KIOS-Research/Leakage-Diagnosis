function [simname] = selectSimulation(data_num)
%% choose a network to load from networks folder
clc
dirName = [pwd,'\results\*.mat'];
Allinpnames = dir(dirName);

if isempty(data_num)
    disp(sprintf('\nChoose scenario file:'))
    for i=1:length(Allinpnames)
        disp([num2str(i),'. ', Allinpnames(i).name])
    end
    x = input(sprintf('\nEnter Data Number: '));
else
    x = data_num;
end
simname=['results\',Allinpnames(x).name];
end


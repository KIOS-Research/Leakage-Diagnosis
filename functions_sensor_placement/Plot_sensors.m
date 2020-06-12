try 
d.unload
catch ERR
end 
fclose all;clear class;clear;clc;close all;
addpath(genpath(pwd));

%% Choose Network
[inpname,dispname] = enterNetwork([]);
d=epanet(inpname);

%% Choose sensors:
sens_indstr{1} = d.getNodeIndex({'78',  '204',  '305',  '322',  '501', '601'});
sens_indstr{2} = d.getNodeIndex({'78',  '219',  '306',  '322',  '420', '601'});
sens_indstr{3} = d.getNodeIndex({'78',  '302',  '306',  '420',  '601', '645'});
% sens_indstr{4} = d.getNodeIndex({'78',  '151',  '172',  '177',  '313',  '322',  '527',  '590'});
% sens_indstr{5} = d.getNodeIndex({'78',  '151',  '172',  '177',  '314',  '322',  '527',  '590'});
% sens_indstr{6} = d.getNodeIndex({'78',  '165',  '306',  '314',  '322',  '420',  '550',  '645'});
% sens_indstr{7} = d.getNodeIndex({'78',  '172',  '177',  '314',  '322',  '550',  '586',  '613'});
% sens_indstr{8} = d.getNodeIndex({'78',  '174',  '188',  '306',  '322',  '420',  '550',  '645'});
sens_indstr{4} =d.getNodeIndex({'78', '302', '314', '601'});



%% Plot sensors:
for i =1:length(sens_indstr)
sens_ind=sens_indstr{i};
d.plot
legend('off')
coor=d.getNodeCoordinates;
x=coor{1}(sens_ind);y=coor{2}(sens_ind);
plot(x,y,'o','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',14)
fontweight='bold';
fontsize=11;
for i=1:length(sens_ind)
text(x(i)-5,y(i),d.getNodeNameID{sens_ind(i)},'Color','black','FontWeight',fontweight,'Fontsize',fontsize)
end
title(['Sensors:',num2str(length(sens_ind))])
tightfig()
end
%%
d.unload
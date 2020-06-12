try 
d.unload
catch ERR
end 
fclose all;clear class;clear;clc;close all;
addpath(genpath(pwd));
load_paths();

%% Load network data
[inpname,dispname] = enterNetwork([]);
d=epanet(inpname);
node_node_link_index=double([d.NodesConnectingLinksIndex d.LinkIndex']);
allParameters=d.getComputedTimeSeries;
Qm = mean(allParameters.Flow)';
wadj=node_node_link_index;
wadj(:,3)=abs(Qm); % [node node weight]

%% Define sensors:
% sens_IDs = {};
load(['functions_sensor_placement\SensorIDs_',dispname])
SensorNodes=double(d.getNodeIndex(sens_IDs));

%% Run decompositin algorithm
MaxCardinality=max(max(wadj(:,1:2))) % maximum number of nodes in a group
MinCardinality=2 % mninimum nodes in a group
RelativeDifference=2 %maximum difference of node number between groups (-1 is not active)
Results=WaterNetPartitioning(wadj,SensorNodes,RelativeDifference,MaxCardinality,MinCardinality) % decomposition algorith
wadj(:,4)=node_node_link_index(:,3);

%% Plot results:
gr = Results.NodeSet;

d.plot
legend('off')
coor=d.getNodeCoordinates;
for i = 1:length(gr)
x=coor{1}(gr{i});y=coor{2}(gr{i});
clrdv = (i/length(gr));
clr = [(1/clrdv)/max(1/clrdv) clrdv 1-clrdv];
plot(x,y,'o','LineWidth',2,'MarkerEdgeColor',1-clr,'MarkerFaceColor',clr,'MarkerSize',14)
end
title(['Node groups: ',num2str(length(gr))])
tightfig()

%% Export results
reservoir = double(d.getNodeReservoirIndex); % remove reservoir node
for i = 1:length(gr)
    gr{i}(find(gr{i}==reservoir))=[];
    grIDs{i} = d.getNodeNameID(gr{i}');
end
save(['functions_group_demands\Node_Groups_',dispname],'grIDs')
%%
d.unload

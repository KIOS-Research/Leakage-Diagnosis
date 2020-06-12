function [ nodeTimeSeries, Qepa, Hepa, LinkStatus ] = DataGenerator( d )
%DEMANDGENERATOR 

%% Epanet Simulation

%calculate with bin
% allParameters=d.getBinComputedAllParameters;
% Hepa = allParameters.BinnodeHead';
% Qepa = allParameters.BinlinkFlow';
% LinkStatus=allParameters.BinlinkStatus';
% LinkStatus(find(LinkStatus==1))=0;
% LinkStatus(find(LinkStatus==3))=1;
% nodeTimeSeries=allParameters.BinnodeDemand';

%calculate
EPAall=d.getComputedHydraulicTimeSeries;
Hepa=EPAall.Head;
Qepa=EPAall.Flow;
LinkStatus=EPAall.Status;
nodeTimeSeries=EPAall.Demand;

%process
resIndex=d.NodeReservoirIndex;
tankIndex=d.NodeTankIndex;
nodeTimeSeries(:,[resIndex, tankIndex])=Hepa(:,[resIndex, tankIndex]);
nodeTimeSeries=nodeTimeSeries';
Hepa=Hepa';
Qepa=Qepa';
LinkStatus=LinkStatus';

% Ktemp = (EPAall.HeadLoss./(abs(EPAall.Flow).^2));
% Kepa = Ktemp(5,:);
% Kepa = Kepa';
end


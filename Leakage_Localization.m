%{
 Copyright (c) 2018 KIOS Research and Innovation Centre of Excellence
 (KIOS CoE), University of Cyprus (www.kios.org.cy)
 
 Licensed under the EUPL, Version 1.1 or – as soon they will be approved
 by the European Commission - subsequent versions of the EUPL (the "Licence");
 You may not use this work except in compliance with theLicence.
 
 You may obtain a copy of the Licence at: https://joinup.ec.europa.eu/collection/eupl/eupl-text-11-12
 
 Unless required by applicable law or agreed to in writing, software distributed
 under the Licence is distributed on an "AS IS" basis,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the Licence for the specific language governing permissions and limitations under the Licence.
 
 Author(s)     : Stelios Vrachimis
 
 Work address  : KIOS Research Center, University of Cyprus
 email         : vrachimis.stelios@ucy.ac.cy (Stelios Vrachimis)
 Website       : http://www.kios.ucy.ac.cy
 
 Last revision : June 2020
%}
%%
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all;
addpath(genpath(pwd));
load_paths();
intv=intvclass2();

%% Load Network
[inpname,dispname] = enterNetwork([]);

%% Define leakage parameters:
if contains(inpname,'Hanoi')    
    % Randomize leakage parameters:
    leakStart= randi([1 24]);
    leak_emit = 2+18*rand; % leak emiter coefficient
    leak_node = randi([1 31]);
end

%% Leakage diagnosis simulation
%%% Initialize:
output=[]; outputSeries=[]; cemitSeries=[]; HleakSeries=[]; QleakSeries=[];displeakbnd=[];DemLeakSeries=[];
flag=0;ind=1;

%%% Begin time-steps:
for demandTime=1:10
    clearvars -except leakStart leak_emit leak_node demandTime inpname dispname qleak0 leakStart HleakSeries QleakSeries DemLeakSeries output outputSeries flag LocalIndex ind scenario displeakbnd
    d=epanet(inpname);
    simDuration = 0; % in hours
    d.setTimeSimulationDuration(simDuration*60*60)
    d.setTimePatternStart(demandTime*60*60) %in seconds
    
%% Define simulation parameters:
if contains(inpname,'Hanoi')    
    % Insert model and demand uncertainty
    Kunc=0.05; %pipe parameter uncertainty
    qunc=0.1; %demand uncertainty

    % Insert demand calibration groups:
    quncG = 0.02; % group uncertainty
    load(['functions_group_demands\Node_Groups_',dispname]) % load node groups grIDs
    for i=1:length(grIDs)
        gr{i} = double(d.getNodeIndex(grIDs{i}));
    end

    % Flow sensors:
    res_ind = d.getNodeReservoirIndex; 
    res_link=find(d.NodesConnectingLinksIndex(:,1)==d.getNodeReservoirIndex);
    if isempty(res_link)
    res_link=find(d.NodesConnectingLinksIndex(:,2)==d.getNodeReservoirIndex);
    end
    flow_meas = res_link; % flow sensor at reservoir links
    eq=0.01; % flow sensor uncertainty
    
    % Pressure sensors:
    load(['functions_sensor_placement\SensorIDs_',dispname]) % load sensor node IDs: sens_IDs
    pres_meas = double(d.getNodeIndex(sens_IDs));
%      pres_meas = double(d.getNodeIndex({'13','15','17','22','27','29'})); % Hanoi
%      pres_meas = double(d.getNodeIndex({'5','13','16','22','30'})); % Hanoi
    ep=0.001; % pressure sensor uncertainty
end

%% Create uncertain model:
%%% Generate node time series and simulate with epanet:
[ nodeValues, Qepa, Hepa, linkStatus ] = DataGenerator( d );
Hepa([d.getNodeReservoirIndex d.getNodeTankIndex],:)=[];
nj=double(d.getNodeJunctionCount);
%%% Check for Epanet negative pressures
elevations=repmat(double(d.getNodeElevations'),[1,size(Hepa,2)]);
elevations([d.getNodeReservoirIndex d.getNodeTankIndex],:)=[];
if any(any((Hepa-elevations)<0))
    error('Epanet returned negative pressures')
end

%%% Pipe resistance coefficients (K, n) and uncertainty:
[ K, n ] = Pipe_Coefficients( d );
%check pipe resistance coefficients
Ktest=K;
Ktest(d.LinkPumpIndex)=[];
if any(Ktest<=0)
    error('Invalid pipe resistance coefficient')
end
clear Ktest

%%% Equal bounds within Kunc:
Kl=K-Kunc*K;
Ku=K+Kunc*K;
Kbnd=[Kl Ku];

%%% Pump coefficients:
[ Pcoef ] = Pump_Coefficients( d );

%%% Initialize bounds:
restankIndex=[d.getNodeTankIndex d.getNodeReservoirIndex];
restankIndex(find(restankIndex==0))=[];
Qlower=zeros(d.LinkCount);
Qupper=zeros(d.LinkCount);
hlower=zeros(d.NodeCount);
hlower(restankIndex,:)=[];
hupper=zeros(d.NodeCount);
hupper(restankIndex,:)=[];

%%% External demands and uncertainty:
qext=nodeValues;
qext(restankIndex)=[];
%Equal bounds within qunc:
ql=qext-qunc*qext;
qu=qext+qunc*qext;
qbnd=[ql, qu];

%%% Reservoir and tank levels:
%restankLevels is a vector size n_n of zeros with values only at reservoir and tank index
restankLevels = nodeValues;
temp=ones(size(restankLevels));
temp(restankIndex)=0;
restankLevels(temp==1)=0;

%%% Model uncertainty:
Lengths = d.getLinkLength;
Lengths = Lengths +  Kunc*Lengths.*(1-2*rand(1,length(Lengths)));
d.setLinkLength(Lengths)

%%% Demand uncertainty:
baseDem = d.getNodeBaseDemands{1};
baseDem = baseDem + qunc*baseDem.*(1-2*rand(1,length(baseDem)));
d.setNodeBaseDemands({baseDem})

%%% Leak start time:
if demandTime>=leakStart
    %%% set leak node:
    emit=d.getNodeEmitterCoeff;
    emit=zeros(size(emit));
    emit(leak_node)=leak_emit; % set emitter coefficient (leakage) value
    d.setNodeEmitterCoeff(emit);
    
    emit=d.getNodeEmitterCoeff;
    if abs(emit(leak_node)-leak_emit)>0.01
        keyboard
    end    
else
    %%% no leak:
    emit=d.getNodeEmitterCoeff;
    emit=zeros(size(emit));
    d.setNodeEmitterCoeff(emit);
end
    
%% Simulate real network:
allParameters=d.getComputedTimeSeries;
Hleak = allParameters.Head';
Hleak(d.getNodeReservoirIndex,:)=[];
Qleak = allParameters.Flow';
DemLeak = emit(leak_node)*allParameters.Pressure(leak_node).^0.5;
DemLeakSeries(ind)=DemLeak;
qleak0 = [DemLeak-40 DemLeak+40]; %leak initial bounds
if qleak0(1)<0; qleak0(1)=0; end

%%% Group uncertainty:
actdem = allParameters.Demand';
if demandTime>=leakStart
actdem(leak_node)=nodeValues(leak_node);
end
actdem(double(d.getNodeReservoirIndex))=[];
Mg = zeros(length(gr),length(actdem));
for i=1:length(gr)
    Mg(i,gr{i}) = 1;
end   
g = Mg*actdem;
gl = g-quncG*g;
gu = g+quncG*g;
gbnd = [gl, gu];

%% Iinitial bounds:
%%% Iinitial bounds when fault-free:
maxFlowUnc = 5*max(qbnd(:,2)-qbnd(:,1));
maxHeadUnc = max(Ku)*maxFlowUnc^n;
Qbnd = [Qepa-maxFlowUnc Qepa+maxFlowUnc];
hbnd = [max(elevations,Hepa-maxHeadUnc) Hepa+maxHeadUnc];

%%% Initial bounds in the case of leak:
maxFlowUnc = 5*max(qbnd(:,2)-qbnd(:,1)) + qleak0(2);
maxHeadUnc = max(Ku)*maxFlowUnc^n;
QbndLeak = [Qepa-maxFlowUnc Qepa+maxFlowUnc];
hbndLeak = [max(elevations,Hepa-maxHeadUnc) Hepa+maxHeadUnc];
if any((QbndLeak(:,1)>Qleak)|((QbndLeak(:,2)<Qleak)))
    keyboard
end
if any((hbndLeak(:,1)>Hleak)|((hbndLeak(:,2)<Hleak)))
    keyboard
end

%% Measurements:
%%% Flow measurements:
Qbnd(flow_meas,:)=[ Qleak(flow_meas)-eq Qleak(flow_meas)+eq ];
QbndLeak(flow_meas,:)=[ Qleak(flow_meas)-eq Qleak(flow_meas)+eq ];
QleakSeries = [QleakSeries Qleak(flow_meas)];

%%% Pressure measurements:
hbnd(pres_meas,:) = [Hleak(pres_meas)-ep Hleak(pres_meas)+ep];
hbndLeak(pres_meas,:) = [Hleak(pres_meas)-ep Hleak(pres_meas)+ep];
HleakSeries = [HleakSeries Hleak(pres_meas)];

%% LDL emitter identification

%%% Detection procedure:
if flag==0 % no leak detected
    qleakbnd=[0 0];
    cemitbnd=[];
    [ Qintvl, Hintvl, qleak, cemit, flag] = Iterations_Emitter( d, n, demandTime, restankLevels, linkStatus, Kbnd, qbnd, Mg, gbnd, Qbnd, hbnd, Pcoef, qleak0, qleakbnd, cemitbnd,leak_node,leak_emit,flag);
    if (flag==1)  && (leakStart>=demandTime)
        disp('LEAK DETECTED CORRECTLY!');
        disp(['Leak start: ',num2str(leakStart)])
        disp(['Leak detection: ',num2str(demandTime)])
        pause(2); 
    end
    output = []; % LDL output
    outputSeries{ind} = zeros(size(hbnd));  
end

%%% Localization procedure:
if flag==1 % leak detected in previous step, continue localization
    qleakbnd=qleak0;
    cemitbnd=output;
   [ Qintvl, Hintvl, qleak, cemit, ~] = Iterations_Emitter( d, n, demandTime, restankLevels, linkStatus, Kbnd, qbnd, Mg, gbnd, QbndLeak, hbndLeak, Pcoef, qleak0, qleakbnd, cemitbnd,leak_node,leak_emit,flag);
    if cemit(leak_node,2)<leak_emit
        keyboard
    end
    output = cemit; % LDL output
    outputSeries{ind} = output;
end

%% Prepare next iteration:
ind=ind+1;
d.unload %unloads libraries and deletes temp files

end 

%% Save and Display results
FileName=['results\LeakScenario_',dispname,datestr(now, '_yyyy-mm-dd_HH-MM')];
save(FileName)
run Display_Results.m

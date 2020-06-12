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
fclose all;clear class;clear;clc;close all;
addpath(genpath(pwd));

%% Choose Network
[inpname,dispname] = enterNetwork([]);
d=epanet(inpname);

Sen_filename = [pwd,'\functions_sensor_placement\SensitivityMat_',dispname,'.mat'];
if isfile(Sen_filename)
% File exists:
loadSen = input(sprintf('\nLoad saved Sensitivity Matrix? (1=yes / 0=no): '));
else
% File does not exist:
loadSen = 0;
end
nj = double(d.getNodeJunctionCount);
nn = double(d.getNodeCount);
reservInd = d.getNodeReservoirIndex;
switch loadSen
    
case 0 % Calculate Sensitivity matrix
%% Calculate healthy states in extended time simulation:
d.setTimeSimulationDuration(24*60*60) % greater weight to low demand hours
d.setTimePatternStart(0) %in seconds
allParameters=d.getComputedTimeSeries;
P0 = allParameters.Pressure';
P0(d.getNodeReservoirIndex,:)=[];
Dem0 = allParameters.Demand';
simSteps = size(P0,2);

%% Create Augmented-time Sensitivity Matrix
%%% Simulate all leakages and get all scenarios pressures
leak_mag_desir=mean(mean(Dem0(Dem0>0)));
mean_pressure = mean(mean(P0(P0>0)));
leak_emit = leak_mag_desir/sqrt(mean_pressure);
emit0=d.getNodeEmitterCoeff;
S = zeros(nj,nj);
for leak=1:nj
    clc
    disp('Calculating Sensitivity Matrix...')
    disp(['Simulating leakage ',num2str(leak),' out of ',num2str(nj)])
    emit=zeros(size(emit0));
    emit(leak)=leak_emit; % set emitter coefficient (leakage) value
    d.setNodeEmitterCoeff(emit);
    allParameters=d.getComputedTimeSeries;
    P = allParameters.Pressure';
    P(d.getNodeReservoirIndex,:)=[];
    Dem = allParameters.Demand';
    leak_mag = Dem(leak,:)-Dem0(leak,:);
    Stmp=(P-P0)./(leak_mag);
    rmax = max(abs(Stmp),[],1);
    Stmp = abs(Stmp)./rmax;
    S(:,leak) = max(Stmp')'; % --> main difference in max-min Sensitivity!!!
end
save(['functions_sensor_placement\SensitivityMat_',dispname],'S')

case 1 % Load residual matrix
    load(['functions_sensor_placement\SensitivityMat_',dispname])
otherwise
    disp('Wrong option selected')
    return
end

%% Place sensors
sensors=9; %number of sensors
exist_sens_ind=[]; % existing sensor indices
disp('Solving using GA...')
[res,fvalSenTr] = SolveEnum2GA(S,nj,sensors,exist_sens_ind);

%% Display solution
u=res;
Sm = S(u>0,:);
Smax = max(Sm);
min_cost = min(Smax);        
mean_cost = mean(Smax);

%%% Sensor indices and IDs
sens_ind = find(u>0.1);
sens_IDs = d.getNodeNameID(sens_ind);

%%% Plot sensors:
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
title(['Sensors:',num2str(sensors)])
tightfig()

%% Save and unload
save(['functions_sensor_placement\SensorIDs_',dispname],'sens_IDs')
d.unload

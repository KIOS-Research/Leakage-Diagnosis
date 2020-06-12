%% Display Results
%  Copyright (c) 2018 KIOS Research and Innovation Centre of Excellence
%  (KIOS CoE), University of Cyprus (www.kios.org.cy)
%  
%  Licensed under the EUPL, Version 1.1 or – as soon they will be approved
%  by the European Commission - subsequent versions of the EUPL (the "Licence");
%  You may not use this work except in compliance with theLicence.
%  
%  You may obtain a copy of the Licence at: https://joinup.ec.europa.eu/collection/eupl/eupl-text-11-12
%  
%  Unless required by applicable law or agreed to in writing, software distributed
%  under the Licence is distributed on an "AS IS" basis,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the Licence for the specific language governing permissions and limitations under the Licence.
%  
%  Author(s)     : Stelios Vrachimis
%  
%  Work address  : KIOS Research Center, University of Cyprus
%  email         : vrachimis.stelios@ucy.ac.cy (Stelios Vrachimis)
%  Website       : http://www.kios.ucy.ac.cy
%  
%  Last revision : June 2020

%% Load results and network:
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all;
addpath(genpath(pwd));
simname = selectSimulation([]);
load(simname)

%% Load network
d=epanet(inpname);

%% Print emmiter bounds change over iterations
if isempty(output)
    AllIsolabilityPercent=0;
else
    AllIsolabilityPercent = 100*((output(:,2))./(sum(output(:,2))));
end
outputSeriesU = [];
for i =1:length(outputSeries)
    outputSeriesU = [outputSeriesU outputSeries{i}(:,2)];
end
outputSeriesU

%% Localization index:
for i =1:length(outputSeries)
    LocalIndex(i) = 100*(length(find(outputSeriesU(:,i)<outputSeriesU(leak_node,i)))/(double(d.getNodeJunctionCount)-1));
end

%% Plot network and flow sensors:
close all
d.plot('nodes','no','highlightlink',d.getLinkNameID(flow_meas),'colorlink',{'g'})
title('(a)')
legend('off')
title('Leak localization')
orange =[0.9100    0.4100    0.1700];
% orange= [(output/max(output)) 1-(output/max(output)) zeros(length(output),1)];
%%% plot sensor nodes:
coor=double(d.getNodeCoordinates(pres_meas));
x=coor(:,1);y=coor(:,2);
plot(x,y,'p','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',16)
%%% plot possible leak nodes (extended time):
coor=double(d.getNodeCoordinates(1:length(output)));
x=coor(:,1);y=coor(:,2);
% plot(x,y,'o','LineWidth',2,'MarkerEdgeColor',orange,'MarkerFaceColor',orange,'MarkerSize',5)
markerSize = 150*(output/max(output));
markerSize(markerSize==0)=0.001;
scatter(x,y,markerSize,orange,'filled') % marker proportional to value
%%% plot leak node:
coor=double(d.getNodeCoordinates(leak_node));
x=coor(:,1);y=coor(:,2);
plot(x,y,'x','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',20)
%%% Print isolability percentage:
fontsize=12;
if isempty(output); output=zeros(size(hbnd)); end
label = round(output(:,2),2);
label = num2str(label);
for node=1:length(output)
coor=d.getNodeCoordinates(node);
x=coor(1);y=coor(2);
xchar = length(label(node,:));
text(x-xchar*42,y-200,label(node,:),'Fontsize',fontsize)
end

%%% Save figure:
% figure2=gcf;
% saveas(figure2,['network.png'])

%% Plot localization index, leak magnitude, pressure and flow measurements:
figure
subplot(4,1,4)
LocalIndex(isnan(LocalIndex))=0;
linewidth=1.5;
plot(LocalIndex,'linewidth',linewidth)
title('(e)')
legend('Localization Index','location','NorthWest')
grid on
axis tight
ylim([0 100])
xlabel('Time (hours)')
ylabel('LI (%)')
% create timestamp for x-axis:
t1 = datetime(2020,1,1,0,0,0);
t2 = datetime(2020,1,2,24,0,0);
data = t1:60*minutes(1):t2;
angle = 60;
set(gca,'xtick',[1:length(LocalIndex)],'xticklabel',datestr(data,'HH:MM'))
xtickangle(angle)
% remove every other x label:
ax = gca;labels = string(ax.XAxis.TickLabels); labels(2:2:length(LocalIndex)) = nan; ax.XAxis.TickLabels = labels;
vline(leakStart-1, '--r', 'Leak Start') % plot leak start vertical line

subplot(4,1,1)
plot(HleakSeries','linewidth',linewidth)
grid on
title('(b)')
% title('Pressure Measurements')
% xlabel('Time (hours)')
ylabel('Head (m)')
legend('Head 13','Head 16','Head 22','Head 30','location','SouthWest')
set(gca,'xtick',[1:length(LocalIndex)],'xticklabel',datestr(data,'HH:MM'))
xtickangle(angle)
% remove every other x label:
ax = gca;labels = string(ax.XAxis.TickLabels); labels(2:2:length(LocalIndex)) = nan; ax.XAxis.TickLabels = labels;
vline(leakStart-1, '--r', 'Leak Start') % plot leak start vertical line
axis tight

subplot(4,1,2)
plot(QleakSeries','linewidth',linewidth)
grid on
title('(c)')
legend('Inflow Measurements','location','NorthWest')
% xlabel('Time (hours)')
ylabel('Flow (CMH)')
set(gca,'xtick',[1:length(LocalIndex)],'xticklabel',datestr(data,'HH:MM'))
xtickangle(angle)
% remove every other x label:
ax = gca;labels = string(ax.XAxis.TickLabels); labels(2:2:length(LocalIndex)) = nan; ax.XAxis.TickLabels = labels;
vline(leakStart-1, '--r', 'Leak Start') % plot leak start vertical line
axis tight

subplot(4,1,3)
%DemLeakSeries(1:leakStart-1)=0;
%DemLeakSeries(DemLeakSeries<0)=0;
plot(DemLeakSeries,'linewidth',linewidth)
grid on
title('(d)')
legend('Leak Magnitude','location','NorthWest')
% xlabel('Time (hours)')
ylabel('Flow (CMH)')
set(gca,'xtick',[1:length(LocalIndex)],'xticklabel',datestr(data,'HH:MM'))
xtickangle(angle)
% remove every other x label:
ax = gca;labels = string(ax.XAxis.TickLabels); labels(2:2:length(LocalIndex)) = nan; ax.XAxis.TickLabels = labels;
vline(leakStart-1, '--r', 'Leak Start') % plot leak start vertical line
axis tight

%%
d.unload
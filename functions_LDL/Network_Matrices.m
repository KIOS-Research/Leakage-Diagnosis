function [ A12, A21, hext, RTL ] = Network_Matrices( d, restankLevels, LinkStatus )
%NETWORK_MATRICES 

%% Create incidence matrix (A12, A21)
node_node_link_index=double([d.NodesConnectingLinksIndex d.LinkIndex']);
A21=zeros(d.NodeJunctionCount,d.LinkCount); %size nodes x pipes
for i=1:length(d.NodesConnectingLinksIndex)
    if  ~ismember(node_node_link_index(i,1),d.NodeReservoirIndex) && ...
        ~ismember(node_node_link_index(i,1),d.NodeTankIndex)&& ...
        LinkStatus(node_node_link_index(i,3))==1
    A21(node_node_link_index(i,1),node_node_link_index(i,3))=-1;
    end
    if ~ismember(node_node_link_index(i,2),d.NodeReservoirIndex) && ...
       ~ismember(node_node_link_index(i,2),d.NodeTankIndex)&& ...
       LinkStatus(node_node_link_index(i,3))==1
    A21(node_node_link_index(i,2),node_node_link_index(i,3))=1;
    end
end;
A12=A21';

%% Idetintify nodes with known head - Reservoirs and tanks (hext = A10*h0)
A10= zeros(d.LinkCount,d.NodeCount);
index=[];
for i=1:length(d.NodeReservoirIndex)
    index=find(node_node_link_index(:,1)==d.NodeReservoirIndex(i));
    if LinkStatus(node_node_link_index(index,3))==1
    A10(index,d.NodeReservoirIndex(i))=-1;
    end
    index=find(node_node_link_index(:,2)==d.NodeReservoirIndex(i));
    if LinkStatus(node_node_link_index(index,3))==1
    A10(index,d.NodeReservoirIndex(i))=1;
    end
end
for i=1:length(d.NodeTankIndex)
    index=find(node_node_link_index(:,1)==d.NodeTankIndex(i));
    if LinkStatus(node_node_link_index(index,3))==1
    A10(index,d.NodeTankIndex(i))=-1;
    end
    index=find(node_node_link_index(:,2)==d.NodeTankIndex(i));
    if LinkStatus(node_node_link_index(index,3))==1
    A10(index,d.NodeTankIndex(i))=1;
    end
end
%Known heads of reservoirs and tanks
h0=restankLevels;% tank head
hext = A10*h0; 

%% Identify reservoir and tank connected links
RTL=zeros(1,d.LinkCount);
for i=1:length(d.NodeReservoirIndex)
    index=find(node_node_link_index(:,1)==d.NodeReservoirIndex(i));
    RTL(node_node_link_index(index,3))=1;
    index=find(node_node_link_index(:,2)==d.NodeReservoirIndex(i));
    RTL(node_node_link_index(index,3))=-1;
end
for i=1:length(d.NodeTankIndex)
    index=find(node_node_link_index(:,1)==d.NodeTankIndex(i));
    RTL(node_node_link_index(index,3))=1;
    index=find(node_node_link_index(:,2)==d.NodeTankIndex(i));
    RTL(node_node_link_index(index,3))=-1;
end

end


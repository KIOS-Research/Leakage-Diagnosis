function [ Qbnd, hbnd, y ] = Initial_Bounds(d,intv,n,restankLevels, LinkStatus,Kbnd,qbnd,Pcoef,Qest,hest)
%INITIAL_BOUNDS 
method = 1; % bounds using interval analysis
% method = 2; % bounds using epanet estimates

if method ==1
%% Define
n_org=n;
Kl=Kbnd(:,1);
Ku=Kbnd(:,2);
qu=qbnd(:,2);
Acoef=Pcoef.Acoef;
Bcoef=Pcoef.Bcoef;
Ccoef=Pcoef.Ccoef;
restankIndex=[d.getNodeTankIndex d.getNodeReservoirIndex];
restankIndex(find(restankIndex==0))=[];

%% Construct Network Matrices
[ A12, A21, hext ] = Network_Matrices( d, restankLevels,LinkStatus );

%% Find initial bounds on heads (h)
Qinit = 1*ones(d.LinkCount,1);
hinit = double(d.NodeElevations');
hinit(restankIndex)=[]; %remove nodes with known head
Qhut = Qinit; hhut = hinit;

% lower bound on heads is node elevations
hlow=hinit;
% upper bound on heads is the sum of:
%reservoir heads, tank heads, pump cutoff heads 
% hup=sum(restankLevels)*ones(size(hlow))+ sum(Ccoef);
hup=sum(restankLevels)*ones(size(hlow))+ sum(Acoef)+30;

%% Find initial bounds on (Q)
hbnd=intv.def(hlow,hup);
Kbnd=intv.def(Kl,Ku);
y=intv.mult(-A12,hbnd);
y=intv.sub(y,hext);
y=intv.div(y,Kbnd);
y=intv.int2mat(y);
hbnd=intv.int2mat(hbnd);
pumpNum=1;

for i=1:size(y,1)
    %link is a pump
    if any(i==d.LinkPumpIndex); 
       Q(i,1)=0;
%        Q(i,2)=20*sum(qu);
       Q(i,2)= (Acoef(pumpNum)/Bcoef(pumpNum))^(1/Ccoef(pumpNum));
       pumpNum=pumpNum+1;
       continue
    end
    
    %link is a pipe
    Q(i,1)=(abs(y(i,1)))^(1/n);
    if y(i,1)<0; Q(i,1)=-Q(i,1);end

    Q(i,2)=(abs(y(i,2)))^(1/n);
    if y(i,2)<0; Q(i,2)=-Q(i,2);end

    lowtemp=min(Q(i,1),Q(i,2)); uptemp=max(Q(i,1),Q(i,2));
    Q(i,1)=lowtemp; Q(i,2)=uptemp;
    n=n_org;
end
Qbnd=Q;
[ ~, ~, ~, RTL ] = Network_Matrices( d, restankLevels,LinkStatus );
Qbnd(:,2)=Qbnd(:,2)+abs(Qbnd(:,2)).*abs(RTL'); %double bounds at reservoir links
end

%% Second method
if method==2

Qbnd(:,1)=Qest-10*abs(Qest);
Qbnd(:,2)=Qest+10*abs(Qest);

restankIndex=[d.getNodeTankIndex d.getNodeReservoirIndex];
restankIndex(find(restankIndex==0))=[];
hinit = double(d.NodeElevations');
hinit(restankIndex)=[];
hbnd(:,1)=hinit;
hbnd(:,2)=hest+abs(hest);

end


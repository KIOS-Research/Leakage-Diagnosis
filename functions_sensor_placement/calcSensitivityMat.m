function S = calcSensitivityMat(inpname)

%% Choose Network
d=epanet(inpname);
nj = double(d.getNodeJunctionCount);

%% Calculate healthy states in extended time simulation:
d.setTimeSimulationDuration(24*60*60) % greater weight to low demand hours
d.setTimePatternStart(0) %in seconds
allParameters=d.getComputedTimeSeries;
P0 = allParameters.Pressure';
P0(d.getNodeReservoirIndex,:)=[];
Dem0 = allParameters.Demand';

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

d.unload
end
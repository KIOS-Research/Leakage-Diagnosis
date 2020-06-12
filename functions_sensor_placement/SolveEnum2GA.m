function [res,fval] = SolveEnum2GA(S,nj,sensors,exist_sens_ind)
%% Initialize parameters
u = ones(1,nj);
lb = zeros(size(u));
ub = ones(size(u));

%%% existing sensors:
if ~isempty(exist_sens_ind)
    Aes = zeros(1,nj);
    Aes(exist_sens_ind) = 1;
    bes = length(exist_sens_ind);
else
    Aes=[];
    bes=[];
end

%% GA options
options = gaoptimset(@ga);
% options.MutationFcn = @mutationgaussian; % select mutation function
options.PlotFcn = @gaplotbestf; % plot the fitness function
options.PopulationSize = 20000; % population in each generation
% options.Generations=1; % set maximum generations
options.StallGenLimit=30; % generation stall limit
options.UseParallel=0; % parallel simulation
% options.InitialPopulation=u;

%% Enumeration loops for threshold and sensors
sen_num_l = sensors(1);
if length(sensors)==1
    sen_num_u = sensors(1);
else
    sen_num_u = sensors(2);
end

%% GA Problem struct
problem = struct(...
'fitnessfcn',@(u)GACostFunction(u,S),... %Fitness function
'nvars',nj,... %Number of design variables
'Aineq',[ones(1,nj); -ones(1,nj); Aes; -Aes],... % contraint on number of sensors
'bineq',[sen_num_u; -sen_num_l; bes; -bes ],...
'Aeq',[],...
'beq',[],...
'lb',lb,...
'ub',ub,...
'nonlcon',[],... %Nonlinear constraint function
'intcon',[1:nj],... %Index vector for integer variables
'options',options,...%Options created with gaoptimset
'rngstate',[]); %state of the random number generator

%%%% Solve GA
% GACostFunction(u,S)
[res,fval] = ga(problem);


end

%%% min sensitivity of all leaks cost function:
function cost = GACostFunction(u,S)
Sm = S(u>0,:);
Smax = max(Sm);
% cost = 1-min(Smax);        
% cost = 1-mean(Smax);         
% cost = 1-0.6*min(Smax)-0.4*mean(Smax);
% cost = 1-0.6*min(Smax)-0.4*median(Smax);
if min(Smax)<0.5
    cost =1;
else
    cost=1-median(Smax);
end
end



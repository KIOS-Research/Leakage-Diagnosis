function [res,fval] = SolveMultiObjGA(S,nj,sensors,exist_sens_ind)
%% Initialize parameters
u0 = ones(1,nj);
lb = zeros(size(u0));
ub = ones(size(u0));

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
options = gaoptimset(@gamultiobj);
% options.MutationFcn = @mutationgaussian; % select mutation function
options.PlotFcn = 'gaplotpareto'; % plot the fitness function
options.PopulationSize = 10000; % population in each generation
options.Display = 'iter';
options.PlotInterval=1;
% options.Generations=1; % set maximum generations
options.StallGenLimit=20; % generation stall limit
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
'fitnessfcn',@(u)GAMultiCostFunction(u,S,sen_num_u),... %Fitness function
'nvars',nj,... %Number of design variables
'Aineq',[ones(1,nj);  Aes; -Aes],... % contraint on number of sensors
'bineq',[sen_num_u;  bes; -bes ],...
'Aeq',[],...
'beq',[],...
'lb',lb,...
'ub',ub,...
'nonlcon',[],... %Nonlinear constraint function
'rngstate',[],... %Optional field to reset the state of the random number generator
'solver','gamultiobj',...
'options',options);%Options created with gaoptimset
% 'intcon',[1:nj],... %Index vector for integer variables
 

%%%% Solve GA
% GACostFunction(u,S)
[res,fval] = gamultiobj(problem);


end

%%% min sensitivity of all leaks cost function:
function cost = GAMultiCostFunction(u,S,sen_num_u)
u=u==1;
% u=round(u);
% if sum(u)>5
%    cost(1) = inf;
%    cost(2) = inf;  
% else
Sm = S(u>0,:);
if isempty(Sm)
   cost(1) = inf;
   cost(2) = inf;
elseif sum(u)<sen_num_u
    cost(1) = 2*sen_num_u-sum(u);
    cost(2) = 2*sen_num_u-sum(u);
else
    Smax = max(Sm); 
    cost(1) = 1-min(Smax);
    cost(2) = 1-mean(Smax);
end
end



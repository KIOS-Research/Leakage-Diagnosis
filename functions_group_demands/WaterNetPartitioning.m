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
 
 Author(s)     : Alexis Kyriacou
 
 Work address  : KIOS Research Center, University of Cyprus
 email         : akyria09@ucy.ac.cy (Alexis Kyriacou)
 Website       : http://www.kios.ucy.ac.cy
 
 Last revision : June 2018

 More info can be found at:
 A. Kyriacou, S. Timotheou, M. P. Michaelides, C. Panayiotou, and M. Polycarpou, 
 “Partitioning of Intelligent Buildings for Distributed Contaminant Detection and Isolation,” 
 IEEE Trans. Emerg. Top. Comput. Intell., vol. 1, no. 2, pp. 72–86, Feb. 2017.
%}

function [Results]=WaterNetPartitioning(InitialMatrix,SensorNodes,RelativeDifference,MaxCardinality,MinCardinality)

%%%%%% Graph appearance settings
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 13);
set(0,'defaultlinelinewidth',2)
set(0,'defaultAxesLineWidth',0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model problem in yalmip for gurobi

%Incident matrix definition

arxikos=InitialMatrix;
D=InitialMatrix;

% Group=CoherencyGroups;

%%%%%%%%%%%%%%%%%Matrix formatting in case it is not ready%%%%%%%%%%%%5
size_D=max(max(arxikos(:,1)),max(arxikos(:,2)));
D=zeros(size_D,size_D);
for i=1:size(arxikos)
    D(arxikos(i,1),arxikos(i,2))=abs(arxikos(i,3));
end

D_arxikos=D;
N=size(D,1);


counter=0;
for i=1:N
    for j=1:i
        if (D_arxikos(i,j)>0)
            counter=counter+1;
        end
    end
end

Index=zeros(1,N,N);
objective_final=zeros(N,1);

i=1;
while i<=max(size(D(:,1)))
    
    if nnz(D(i,:))<1 && nnz(D(:,i))<1
        D(i,:)=[];
        D(:,i)=[];
        i=i-1;%%%%%%%%% Auto gia na min peida grammes
    end
    i=i+1;
end

N=size(D,1);
D=D(1:N,1:N);

for i=1:N
    D(i,i)=0;
end


% if ~issymmetric(D)
%      msg='Re gare en ine simmetrikos. Mallon enna exei ke tuto diples grammes!!!!'
% 
% for i=1:length(D)
%     for j=1:length(D)
%         if (D(i,j)~=D(j,i))
%             i
%             j
%         end
%     end
% end
%             err(msg)
% end
  D=D+D.';

% MaxSize=size(D,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Existing edges %%%%%%%%%%%
new_matrix=sparse(N,N);
index=1;
for i=1:N
    for j=1:N
        if D(i,j)~=0
            new_matrix(i,j)=index;
            new_D(index)=D(i,j);
            index=index+1;
        end
    end
end
index=index-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=length(SensorNodes);
% MaxCardinality=MaxCardinality;

clear('yalmip');
%% Variable Definition
w=sdpvar(index,K,'full');
c=sdpvar(index,1,'full');


    x=binvar(N,K,'full');

%%%%%%%%%
y=sdpvar(N,K,'full');
u=sdpvar(N,K,'full');
B=sdpvar(1,K,'full');
g=sdpvar(N,K,'full');
f=sdpvar(index,K,'full');

%% Main Partitioning Constraints
constraints = [];
constraints = [constraints,0<=w<=1];
constraints = [constraints, 0<=c<=1];

vectmp = [];
index2= 1;
for i=1:N
    for j=1:i
        if new_matrix(i,j)~=0
            vectmp(index2,1) = new_matrix(i,j);
            vectmp(index2,2) = new_matrix(j,i);
            index2 = index2 + 1;
        end
    end
end

constraints = [constraints, c(vectmp(:,1))==c(vectmp(:,2))];
index=1;
for i=1:N
    for j=1:N
        if new_matrix(i,j)~=0
            new_array(index,1)=new_matrix(i,j);
            new_array(index,2)=i;
            new_array(index,3)=j;
            index=index+1;
        end
    end
end

for h=1:K
    constraints = [constraints, w(new_array(:,1),h)<=x(new_array(:,2),h)];
    constraints = [constraints, w(new_array(:,1),h)<=x(new_array(:,3),h)];
end

% Xi,k constrains
constraints = [constraints,c==sum(w,2)];
%%% Changed section for the Shared nodes
constraints = [constraints,sum(x.')==1]; %Shows that a vertex can only be included in one cluster
%         constraints = [constraints,sum(x(4,:))<==2];
% constraints = [constraints,sum(x)<=MaxCardinality];
%% Connectivity Constraints
%         Yj,h constrains
WW = zeros(N,N);
for i=1:N
    j=1:i;
    WW(i,j)=1;
end

for h=1:K
    constraints = [constraints,y(:,h)<=WW*x(:,h)];
    constraints = [constraints,y(:,h)>=(WW./N)*x(:,h)];
end
constraints = [constraints,y>=x];

constraints = [constraints,0<=y<=1];





%         Find the node with the smallest index inside a group\

for j=2:N
    for h=1:K
        constraints = [constraints,u(j,h)==y(j,h)-y(j-1,h)];
        
    end
end

constraints = [constraints,u(1,:)==y(1,:)];

constraints = [constraints,sum(u,1)==1];
constraints = [constraints,u<=x];
constraints = [constraints,0<=u<=1];

% Connectivity Linearization

constraints = [constraints,B==sum(x,1)];
constraints = [constraints,MinCardinality<=B<=MaxCardinality];
% constraints = [constraints,MaxCardinality>=sum(x)>=MinCardinality]; % groups cardinality


%%% Equalization of the cluster size with a difference between clusters
if RelativeDifference~=-1
    for h=1:K
        constraints = [constraints,(B(h)-B([(1:(h-1)),((h+1):K)]))<=RelativeDifference];
        constraints = [constraints,(B(h)-B([(1:(h-1)),((h+1):K)]))>=-RelativeDifference];
    end
end


constraints = [constraints,g<=(N*u)];
constraints = [constraints,g<=ones(N,1)*B-(1-u)*MinCardinality];
constraints = [constraints,g>=ones(N,1)*B-(1-u)*N];
constraints = [constraints,g>=(MinCardinality*u)];

% Linearized equation

for j=1:N
    indtmp = find(new_matrix(:,j)>0);
    vecInd_vj = new_matrix(indtmp ,j);
    indtmp = find(new_matrix(j,:)>0);
    vecInd_jv = new_matrix(j,indtmp);
    constraints = [constraints,(1/N)*(g(j,:)-x(j,:))+sum(f(vecInd_vj,:),1)==sum(f(vecInd_jv,:),1)];
end

for h=1:K
    constraints = [constraints,0<=f(:,h)<=c];
end

% Sensor Nodes assignment
%% Generator Assignmnet
for h=1:K    
    constraints = [constraints,x(SensorNodes(h),h)==1];
end

    objective = 0.5*sum((1-c(1:size(new_D,2),1)).*new_D');


%% Solver Options
options = sdpsettings('verbose',1,'solver','gurobi');
%%% Time limit in case it is needed
options.gurobi.method=3;
options.gurobi.MIPGap=0.00;

%    options.gurobi.DisplayInterval=30;
% problem solution
initial_constraints=constraints;
diagnostics = solvesdp(constraints,objective,options);


%% Solution Results and Statistics Gathering
MilpTime=diagnostics.solvertime;
ResultingGroups=double(x);
temp=double(x);

Results.ResultingGroups=ResultingGroups;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Power Islands Calculation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Find the nodes of each island
ResultingIslands=[];
for i=1:size(temp,2)
    counter=1;
    for j=1:size(temp,1)
        if (temp(j,i)>0)
            ResultingIslands(i,counter)=j;
            counter=counter+1;
        end
    end
end

%% CutSet Extraction

% INITIALIZATION
CutSet={};
for i=1:size(ResultingIslands,1)
    
    CutSet{i}=[];
end

full_Initial=D;

Cost=zeros(size(ResultingIslands,1),1);
full_Initial=full(InitialMatrix);
for k=1:size(ResultingIslands,1)-1
    for n1=ResultingIslands(k,:)
        if n1==0
            break;
        else
            for c=k+1:size(ResultingIslands,1)
                for n2=ResultingIslands(c,:)
                    if n2==0
                        break;
                    else
                        if D(n1,n2)~=0
                            CutSet{c-1}(end+1,1)=n1;
                            CutSet{c-1}(end,2)=n2;
                            Cost(c-1)=Cost(c-1)+D(n1,n2);
                        end
                    end
                end
            end
        end
    end
end
Cost_total=sum(Cost);



% Results.CutSet=CutSet;
% Results.Cost=Cost;
% Results.ResultingIslands=ResultingIslands;
% Results.MILP_Time=overall_time;
% Results.Recursive=Recursive;
% Results.overall_recursion_overhead=overall_recursion_overhead;
% Results.added_load_nodes=added_load_nodes;
Results.Cost_total=Cost_total;
Results.CutSet=CutSet;
for h=1:K
    Results.NodeSet{h}=find(temp(:,h));
end
end


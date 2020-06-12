function [Q, h, qleak, cemitbnd, skipCriterion] = Iterations_Emitter(d,n,timeStep,restankLevels,LinkStatus,Kbnd,qbnd,Mg,gbnd,Qbnd,hbnd,Pcoef,...
                                                                            qleak0,qleakbnd,cemitbnd,leak_node,leak_emit,loc)
%{
 LDL Iterations

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
%% Initialize
Kl=Kbnd(:,1);
Ku=Kbnd(:,2);
qextl=qbnd(:,1);
qextu=qbnd(:,2);
gl=gbnd(:,1);
gu=gbnd(:,2);
Acoef=Pcoef.Acoef;
Bcoef=Pcoef.Bcoef;
Ccoef=Pcoef.Ccoef;
pumpInd=d.getLinkPumpIndex;
np=length(Qbnd);
nu=length(hbnd);
ngc=length(gl);
elev = d.getNodeElevations';
elev(d.getNodeReservoirIndex)=[];
alpha=0.5; % leak emitter exponent
iter=2;

%% Leak matrices:
if loc==0
    qleak0=[0 0];
    qleakbnd=[0 0];
end
nleak=nu;
Bleak = -eye(nu,nleak);
qleakLB = qleak0(1)*ones(size(hbnd,1),1);

if length(qleakbnd)<nu
    qleakLB = qleakbnd(1)*ones(size(hbnd,1),1);
    qleakbnd(1) = 0;
    qleakbnd=repmat(qleakbnd,nu,1);
end

if isempty(cemitbnd)
    cemit = zeros(size(qleakbnd));
    cemit(:,1) = qleakbnd(:,1)./ ((hbnd(:,1)-elev).^alpha);
    cemit(isnan(cemit(:,1)),1)=0;
    cemit(:,2) = qleakbnd(:,2)./ ((hbnd(:,2)-elev).^alpha);
    if (cemit(leak_node,2)<leak_emit)&&loc
%         keyboard
    end
else
    cemit=cemitbnd;
    qleakbnd(:,1)=cemitbnd(:,1).*((hbnd(:,1)-elev).^alpha);
    qleakbnd(:,2)=cemitbnd(:,2).*((hbnd(:,2)-elev).^alpha);
end

%% Construct Network Matrices
[ A12, A21, hext, RTL ] = Network_Matrices( d, restankLevels,LinkStatus  );
    
%% Find initial approximate interval lines
linkInd=1:length(Qbnd);
[ slope, beta ] = Interval_Lines(n, Kl, Ku, Qbnd, linkInd, Pcoef, pumpInd, iter);

%% Define Lambda, A and B matrices 
Ll = diag(slope(:,1));
Lu = diag(slope(:,2));

Aeq = [eye(np)  zeros(np,np+nleak)  A12];
Aeq=sparse(Aeq);
Beq = -hext;


% Vleak=ones(1,nleak);
% bVleak=max(max(qleakbnd));

Vleak = (1./qleakbnd(:,2))';
Vleak(isinf(Vleak))=0;
bVleak=1;

Aineq =[zeros(nu,np)   A21     Bleak    zeros(nu,nu);
        zeros(nu,np)  -A21    -Bleak    zeros(nu,nu);
        zeros(1,np)    zeros(1,np) Vleak  zeros(1,nu);
        zeros(ngc,np)  Mg*A21   Mg*Bleak  zeros(ngc,nu);
        zeros(ngc,np) -Mg*A21  -Mg*Bleak  zeros(ngc,nu);
        -eye(np)       Ll      zeros(np,nu+nleak);
         eye(np)      -Lu      zeros(np,nu+nleak)];

Aineq=sparse(Aineq);
Bineq =[  qextu
         -qextl
          bVleak
          gu
         -gl
         -beta(:,1)
          beta(:,2)];

zbnd(:,1) = (slope(:,1).*Qbnd(:,1))+beta(:,1);
zbnd(:,2) = (slope(:,2).*Qbnd(:,2))+beta(:,2);
LB = [ zbnd(:,1);
       Qbnd(:,1);
       qleakbnd(:,1);
       hbnd(:,1)];
UB = [ zbnd(:,2);
       Qbnd(:,2);
       qleakbnd(:,2);
       hbnd(:,2)];    

model.A = [Aineq;Aeq];
model.rhs = [Bineq;Beq];
model.lb = LB;
model.ub = UB;
sizeIneq = size(Aineq,1);
sizeEq   = size(Aeq,1);
model.sense(1:sizeIneq) = '<';
model.sense((sizeIneq+1):(sizeIneq + sizeEq)) = '=';
model.modelsense = 'min';
params.method = 2;
params.OutputFlag = 0;
    
%% Iterations
tolerance=1;
maxIter=20;
tic
iter=2;
criterion(1)=-inf;
criterion(2)=inf;
while  (iter < maxIter) && (abs(criterion(iter)-criterion(iter-1)) > tolerance)
%% Initialize:
    iter=iter+1;
    skipCriterion = 0;
    Q=Qbnd;
    h=hbnd;    
    qleak=qleakbnd;
    
%% Run linear program
    % Flows:
    for i=(np+1):(2*np)
        if mod(i,10)==0
        clc
        if loc disp('Localization phase'); else disp('Detection phase'); end
        disp(['Time step: ',num2str(timeStep)])
        disp(['Iteration: ',num2str(iter)])
        disp(['State (Pipe) ',num2str(i-np),' of ',num2str(np+nu+nleak)])
        end
        
        if any(LB>UB); error('Lower bound is larger than upper bound');end
        % Minimization Pipe (i-np):
        W=zeros(1,2*np+nu+nleak);
        W(i)=1;
        model.obj = W;
        result = gurobi(model,params);        
        if strcmp(result.status,'OPTIMAL')
            xl=result.x;
            Q(i-np,1)=xl(i);
            LB(i)=xl(i);
            LB(i-np)=slope(i-np,1)*xl(i)+beta(i-np,1);
            model.lb = LB;            
        else
            skipCriterion = 1;
            disp(result.status) 
        end
        
        % Maximization Pipe (i-np):
        model.obj = -W;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xu=result.x;
            Q(i-np,2)=xu(i);
            UB(i)=xu(i);
            UB(i-np)=slope(i-np,2)*xu(i)+beta(i-np,2);
            model.ub = UB;
        else
            skipCriterion = 1;
            disp(result.status)
        end
        
        % Recalculate pipe approximate lines:
        linkInd=(i-np);
        [ slope((linkInd),:), beta((linkInd),:)] = Interval_Lines(n, Kl(linkInd), Ku(linkInd), Q(linkInd,:), linkInd, Pcoef, pumpInd, iter);
        
        % Redifine inequality matrices:
        Ll = diag(slope(:,1));
        Lu = diag(slope(:,2));
        Aineq =[zeros(nu,np)  A21       Bleak      zeros(nu,nu);
                zeros(nu,np) -A21      -Bleak      zeros(nu,nu);
                zeros(1,np)    zeros(1,np)  Vleak zeros(1,nu);
                zeros(ngc,np)  Mg*A21   Mg*Bleak   zeros(ngc,nu);
                zeros(ngc,np) -Mg*A21  -Mg*Bleak   zeros(ngc,nu);
                -eye(np)      Ll        zeros(np,nu+nleak);
                 eye(np)     -Lu        zeros(np,nu+nleak)];
        Aineq=sparse(Aineq);
        Bineq =[  qextu
                 -qextl
                  bVleak
                  gu
                 -gl
                 -beta(:,1)
                  beta(:,2)];
        model.A = [Aineq;Aeq];
        model.rhs = [Bineq;Beq];
    end

    % Heads:
    for i=(2*np+nleak+1):(2*np+nu+nleak)
        if mod(i,10)==0
        clc
        if loc disp('Localization phase'); else disp('Detection phase'); end
        disp(['Time step: ',num2str(timeStep)])
        disp(['Iteration: ',num2str(iter)])
        disp(['State (Node) ',num2str(np+ i-(2*np+nleak)),' of ',num2str(num2str(np+nu+nleak))])
        end
        
        if any(LB>UB); error('Lower bound is larger than upper bound');end
        % Minimization node (i-(2*np+nleak)):
        W=zeros(1,2*np+nu+nleak);
        W(i)=1;
        model.obj = W;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xl=result.x;
            h(i-(2*np+nleak),1)=xl(i);
            LB(i)=xl(i);
            model.lb = LB;           
        else
            skipCriterion = 1;
            disp(result.status)
        end
        
        % Maximization node (i-(2*np+nleak)):
        model.obj = -W;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xu=result.x;            
            h(i-(2*np+nleak),2)=xu(i);
            UB(i)=xu(i);
            model.ub = UB;             
        else
            skipCriterion = 1;
            disp(result.status)
        end
        
    end
        
    % Update leak bounds:
    qleak(:,1)=cemit(:,1).*((h(:,1)-elev).^alpha);
    qleak(:,2)=cemit(:,2).*((h(:,2)-elev).^alpha);
    LB((2*np+1):(2*np+nleak))=qleak(:,1);
    UB((2*np+1):(2*np+nleak))=qleak(:,2);
    model.ub=UB;
    model.lb=LB;    
    
    % Leakages:
    UBtemp = UB((2*np+1):(2*np+nleak));
    LBtemp = LB((2*np+1):(2*np+nleak));
    for i=(2*np+1):(2*np+nleak)
        if mod(i,10)==0
        clc
        if loc disp('Localization phase'); else disp('Detection phase'); end
        disp(['Time step: ',num2str(timeStep)])
        disp(['Iteration: ',num2str(iter)])
        disp(['State (Leak) ',num2str(np+nu+ i-(2*np)),' of ',num2str(np+nu+nleak)])
        end
        
        if any(LB>UB); error('Lower bound is larger than upper bound');end
        
        % Maximization leak (i-2*np) when all other leaks are zero:
        if UBtemp(i-(2*np))==0; continue; end % skip if upper bound is zero
        W=zeros(1,2*np+nu+nleak);
        W(i)=-1;
        model.obj = W;
        UB((2*np+1):(2*np+nleak))=0;
        UB(i)=UBtemp(i-(2*np));
        LB(i)=qleak0(1);
        model.ub=UB;
        model.lb=LB;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xu=result.x;
            UBtemp(i-2*np)=xu(i);
            qleak(i-2*np,2)=xu(i);
        else
            UBtemp(i-2*np)=0;
            LBtemp(i-2*np)=0;
            qleak(i-2*np,:)=0;
            if (qleak(leak_node,2)==0)&&loc
                keyboard
            end 
        end
        UB((2*np+1):(2*np+nleak))=UBtemp;
        LB((2*np+1):(2*np+nleak))=LBtemp;
        model.lb=LB;
        model.ub=UB;
        
        % Minimize leak (i-2*np) when all other leaks are zero:
        if UBtemp(i-(2*np))==0; continue; end % skip if upper bound is zero
        W=zeros(1,2*np+nu+nleak);
        W(i)=1;
        model.obj = W;
        UB((2*np+1):(2*np+nleak))=0;
        UB(i)=UBtemp(i-(2*np));
        LB(i)=qleak0(1);
        model.ub=UB;
        model.lb=LB;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xl=result.x;
            qleakLB(i-2*np)=xl(i); % used here in next iteration
        else
            UBtemp(i-2*np)=0;
            LBtemp(i-2*np)=0;
            qleak(i-2*np,:)=0; 
            if (qleak(leak_node,2)==0)&&loc
                keyboard
            end 
        end
        UB((2*np+1):(2*np+nleak))=UBtemp;
        LB((2*np+1):(2*np+nleak))=LBtemp;
        model.lb=LB;
        model.ub=UB;
    end
    
    % Update emitter bounds:
    cemit(:,1)=qleak(:,1)./((h(:,1)-elev).^alpha);
    cemit(:,2)=qleak(:,2)./((h(:,2)-elev).^alpha); 
    if (cemit(leak_node,2)<leak_emit) && loc
%         keyboard
    end    
    %% Redefine variable bounds
    Qbnd(:,1)=Q(:,1);
    Qbnd(:,2)=Q(:,2);
    qleakbnd(:,1)=qleak(:,1);
    qleakbnd(:,2)=qleak(:,2);
    Vleak = (1./qleakbnd(:,2))';
    Vleak(isinf(Vleak))=0;
    hbnd(:,1)=h(:,1);
    hbnd(:,2)=h(:,2);
    if ~isempty(cemitbnd)
    cemitbnd(:,1)=max(cemit(:,1),cemitbnd(:,1));
    cemitbnd(:,2)=min(cemit(:,2),cemitbnd(:,2));
    else
    cemitbnd(:,1)=cemit(:,1);
    cemitbnd(:,2)=cemit(:,2);    
    end
    if skipCriterion==0
        criterion(iter) = sum([Q(:,2); h(:,2)] - [Q(:,1); h(:,1)]); %1-norm of Q and h
    else
        criterion(iter) = criterion(iter-1)+2*tolerance;
    end

end

end


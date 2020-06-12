function [ slope, beta ] = Interval_Lines(n, Kl, Ku, Q, linkInd, Pcoef, pumpInd, iter)
%INTERVAL_LINES
%%Pump coefficients
Acoef=Pcoef.Acoef;
Bcoef=Pcoef.Bcoef;
Ccoef=Pcoef.Ccoef;
    
%% Find approximate lines
np=size(Q,1);
slope=zeros(np,2);
beta=zeros(np,2);
% pumpNum=1;
for i=1:np
    y1=[];y2=[];
    %%%link is a PUMP
    if any(linkInd(i)==pumpInd); 
        pumpNum=find(pumpInd==linkInd(i)); %find which pump
        x1=Q(i,1);
        y(i,1)=Acoef(pumpNum)-Bcoef(pumpNum)*x1^Ccoef(pumpNum);
%         y(i,1)=Acoef(pumpNum)*x1^2+Bcoef(pumpNum)*x1+Ccoef(pumpNum);
%         y(i,1)=Acoef(pumpNum)*x1^2+Ccoef(pumpNum);
        x2=Q(i,2);
        y(i,2)=Acoef(pumpNum)-Bcoef(pumpNum)*x2^Ccoef(pumpNum);
%         y(i,2)=Acoef(pumpNum)*x2^2+Bcoef(pumpNum)*x2+Ccoef(pumpNum);
%         y(i,2)=Acoef(pumpNum)*x2^2+Ccoef(pumpNum);

        slope(i,1)=(y(i,2)-y(i,1))/(Q(i,2)-Q(i,1)); %strictly decreasing function
        slope(i,2)=slope(i,1);
        if (y(i,2)-y(i,1))==0; slope(i,1)=0; end
        
        x1=Q(i,1);
        x2=Q(i,2);
        x3= + abs(slope(i,1)/(-Bcoef(pumpNum)*Ccoef(pumpNum)))^(1/(Ccoef(pumpNum)-1));
        x4= - abs(slope(i,1)/(-Bcoef(pumpNum)*Ccoef(pumpNum)))^(1/(Ccoef(pumpNum)-1));
%         x3=(slope(i,1)-Bcoef(pumpNum))/(2*Acoef(pumpNum));
%         x3=(slope(i,1))/(2*Acoef(pumpNum));

        fx1=Acoef(pumpNum)-Bcoef(pumpNum)*x1^Ccoef(pumpNum)-slope(i,1)*x1;
%         fx1=Acoef(pumpNum)*x1^2+Bcoef(pumpNum)*x1+Ccoef(pumpNum)-slope(i,1)*x1;
%         fx1=Acoef(pumpNum)*x1^2+Ccoef(pumpNum)-slope(i,1)*x1;

        fx2=Acoef(pumpNum)-Bcoef(pumpNum)*x2^Ccoef(pumpNum)-slope(i,1)*x2;
%         fx2=Acoef(pumpNum)*x2^2+Bcoef(pumpNum)*x2+Ccoef(pumpNum)-slope(i,1)*x2;
%         fx2=Acoef(pumpNum)*x2^2+Ccoef(pumpNum)-slope(i,1)*x2;

        if (x3<Q(i,1) || x3>Q(i,2)); fx3=fx1;else
            fx3=Acoef(pumpNum)-Bcoef(pumpNum)*x3^Ccoef(pumpNum)-slope(i,1)*x3;            
%             fx3=Acoef(pumpNum)*x3^2+Bcoef(pumpNum)*x3+Ccoef(pumpNum)-slope(i,1)*x3;
%             fx3=Acoef(pumpNum)*x3^2+Ccoef(pumpNum)-slope(i,1)*x3;
        end
        
        if (x4<Q(i,1) || x4>Q(i,2)); fx4=fx2;else
            fx4=Acoef(pumpNum)-Bcoef(pumpNum)*x4^Ccoef(pumpNum)-slope(i,1)*x4;
        end        
        
        beta(i,1)=min([fx1;fx2;fx3;fx4]);
        beta(i,2)=max([fx1;fx2;fx3;fx4]); 
        
        %Reverse bounds because pump gives head
        temp1=slope(i,1);
        temp2=beta(i,1);
        slope(i,1)=-slope(i,2);
        beta(i,1)=-beta(i,2);
        slope(i,2)=-temp1;
        beta(i,2)=-temp2;
        
        continue
    end

    %link is a PIPE
    x1=Q(i,1); %q_l
    x2=Q(i,2); %q_u
    y1(1)=Kl(i)*abs(x1)^(n-1)*x1;
    y1(2)=Kl(i)*abs(x2)^(n-1)*x2; 
    y2(1)=Ku(i)*abs(x1)^(n-1)*x1;
    y2(2)=Ku(i)*abs(x2)^(n-1)*x2;
    
    % check if Q interval is thin
    if (x1==x2)
        slope(i,1)=0;
        slope(i,2)=0;
        beta(i,1)=min(y1(1),y2(1));
        beta(i,2)=max(y1(1),y2(1));
        continue
    end
    
    %%% Lower line
    slope(i,1)=(min(y1(2),y2(2))-min(y1(1),y2(1)))/(x2-x1);
    % First function
    fx1=Kl(i)*abs(x1)^(n-1)*x1-slope(i,1)*x1;
    fx2=Kl(i)*abs(x2)^(n-1)*x2-slope(i,1)*x2;
    x3=-(slope(i,1)/(Kl(i)*n))^(1/(n-1)); 
    if (x3<x1 || x3>x2); fx3=fx1;else
    fx3=Kl(i)*abs(x3)^(n-1)*x3-slope(i,1)*x3;
    end
    x4= (slope(i,1)/(Kl(i)*n))^(1/(n-1));
    if (x4<x1 || x4>x2); fx4=fx2;else
    fx4=Kl(i)*abs(x4)^(n-1)*x4-slope(i,1)*x4;
    end
    % Second function
    fx5=Ku(i)*abs(x1)^(n-1)*x1-slope(i,1)*x1;
    fx6=Ku(i)*abs(x2)^(n-1)*x2-slope(i,1)*x2;
    x7=-(slope(i,1)/(Ku(i)*n))^(1/(n-1));
    if (x7<x1 || x7>x2); fx7=fx5;else
    fx7=Ku(i)*abs(x7)^(n-1)*x7-slope(i,1)*x7;
    end
    x8= (slope(i,1)/(Ku(i)*n))^(1/(n-1));
    if (x8<x1 || x8>x2); fx8=fx6;else
    fx8=Ku(i)*abs(x8)^(n-1)*x8-slope(i,1)*x8;
    end
    % Lower line offset
    beta(i,1)=min([fx1;fx2;fx3;fx4;fx5;fx6;fx7;fx8]);

    %%% Upper line
    slope(i,2)=(max(y1(2),y2(2))-max(y1(1),y2(1)))/(x2-x1);        
    % First function
    fx1=Kl(i)*abs(x1)^(n-1)*x1-slope(i,2)*x1;
    fx2=Kl(i)*abs(x2)^(n-1)*x2-slope(i,2)*x2;
    x3=-(slope(i,2)/(Kl(i)*n))^(1/(n-1)); 
    if (x3<x1 || x3>x2); fx3=fx1;else
    fx3=Kl(i)*abs(x3)^(n-1)*x3-slope(i,2)*x3;
    end
    x4= (slope(i,2)/(Kl(i)*n))^(1/(n-1));
    if (x4<x1 || x4>x2); fx4=fx2;else
    fx4=Kl(i)*abs(x4)^(n-1)*x4-slope(i,2)*x4;
    end
    % Second function
    fx5=Ku(i)*abs(x1)^(n-1)*x1-slope(i,2)*x1;
    fx6=Ku(i)*abs(x2)^(n-1)*x2-slope(i,2)*x2;
    x7=-(slope(i,2)/(Ku(i)*n))^(1/(n-1));
    if (x7<x1 || x7>x2); fx7=fx5;else
    fx7=Ku(i)*abs(x7)^(n-1)*x7-slope(i,2)*x7;
    end
    x8= (slope(i,2)/(Ku(i)*n))^(1/(n-1));
    if (x8<x1 || x8>x2); fx8=fx6;else
    fx8=Ku(i)*abs(x8)^(n-1)*x8-slope(i,2)*x8;
    end
    % Upper line offset
    beta(i,2)=max([fx1;fx2;fx3;fx4;fx5;fx6;fx7;fx8]);

end
    
%% Line check - upper line is above lower line within feasible range of Q
if any((slope(:,1).*Q(:,1)+beta(:,1))>(slope(:,2).*Q(:,1)+beta(:,2)))||...
   any((slope(:,1).*Q(:,2)+beta(:,1))>(slope(:,2).*Q(:,2)+beta(:,2))) 
   error('Wrong line approximation.')
end


%% Plot approximate lines    
flag=0;
if flag==1
p=1; %pump number (1,2)
if linkInd== 6%pumpInd(p)

%     
% if np>9
%     plotperfig =9;
% else
%     plotperfig=np;
% end
% pumpNum=1;
% for j=1:ceil(np/plotperfig)
% figure('units','normalized','outerposition',[0 0 0.95 1])
% 
% for i=1:plotperfig
%     subplot(3,3,i)
%     i=i+(j-1)*plotperfig;
%     if i>np; break; end
if ismember(iter,[3,4,5,6,7,8])
    spnum=(iter-2);
    subplot(2,3,spnum)
    xt=linspace(Q(i,1),Q(i,2),1000);
    if linkInd==pumpInd(p); 
        pumpNum=find(pumpInd==linkInd(i));
        yt1=-(ones(1,length(xt))*Acoef(pumpNum)-Bcoef(pumpNum)*xt.^Ccoef(pumpNum));
        yt2=-(ones(1,length(xt))*Acoef(pumpNum)-Bcoef(pumpNum)*xt.^Ccoef(pumpNum));
%         yt1=-(Acoef(pumpNum)*xt.^2+ Bcoef(pumpNum)*xt + ones(1,length(xt))*Ccoef(pumpNum));
%         yt2=-(Acoef(pumpNum)*xt.^2+ Bcoef(pumpNum)*xt + ones(1,length(xt))*Ccoef(pumpNum));
%         yt1=-(Acoef(pumpNum)*xt.^2+ones(1,length(xt))*Ccoef(pumpNum));
%         yt2=-(Acoef(pumpNum)*xt.^2+ones(1,length(xt))*Ccoef(pumpNum));
    else
        yt1=Kl(i)*abs(xt).^(n-1).*xt;
        yt2=Ku(i)*abs(xt).^(n-1).*xt;
    end  

    plot(xt,yt1,'r','LineWidth',1.2)
    hold all
    plot(xt,yt2,'Color',[1 0.7 0.1],'LineWidth',1.2)
    yt=slope(i,1).*xt+beta(i,1);
    plot(xt,yt,'--','Color',[0 0.5 0.1],'LineWidth',1.2)
    yt=slope(i,2).*xt+beta(i,2);
    plot(xt,yt,'--b','LineWidth',1.2)    
    xlim([Q(i,1) Q(i,2)])
    xlabel(['q_',num2str(linkInd)],'FontSize',12)
    ylabel(['{\bff_',num2str(linkInd),'}( q_',num2str(linkInd),' )'],'FontSize',12)
    axis tight
    %Area between lines:
    ql=Q(1);
    qu=Q(2);
    fl_ql= slope(1)*ql+beta(1);
    fl_qu= slope(1)*qu+beta(1);
    fu_ql= slope(2)*ql+beta(2);
    fu_qu= slope(2)*qu+beta(2);
    Area= 0.5* (qu-ql)* (fu_qu+fu_ql-fl_qu-fl_ql);
    title(['Iteration ',num2str(iter-2),' , Unc.Area = ', num2str(Area)])
    grid on
%     legend('K_l*f(q)','K_u*f(q)', 'y_l=\lambda_l q + b_l', 'y_u=\lambda_u q + b_u','Location','Northwest')
end
end
end

end

function [ Pcoef ] = Pump_Coefficients( d )
%PUMP_COEFFICIENTS
%Can handle single point curves and 3-point curves

Acoef=[];
Bcoef=[];
Ccoef=[];

%Calculate coefficients with lsqcurvefit
pumpIndex=d.getCurveIndex;
if ~isempty(pumpIndex)
    for i=1:d.LinkPumpCount
        points=d.getCurveValue(pumpIndex(i));
        xdata=points(:,1);
        ydata=points(:,2);
        % add points when a single point curve is given:
        if length(xdata)==1
            pointx = xdata;
            pointy = ydata;
            xdata = zeros(3,1);
            ydata = zeros(3,1);
            ydata(1)=1.33*pointy;%shutoff head at zero flow equal to 133% of the design head
            xdata(2)=pointx;
            ydata(2)=pointy;
            xdata(3)=2*pointx;%maximum flow at zero head equal to twice the design flow
            x0=[50 0.1];
            CurPar=lsqcurvefit(@(x,xdata) x(1)-x(2)*xdata.^2,x0,xdata,ydata);
            Acoef(i)=CurPar(1);
            Bcoef(i)=CurPar(2);
            Ccoef(i)=2;
            continue
        end
        x0=[50 0.1 2];
        CurPar=lsqcurvefit(@(x,xdata) x(1)-x(2)*xdata.^x(3),x0,xdata,ydata);
        Acoef(i)=CurPar(1);
        Bcoef(i)=CurPar(2);
        Ccoef(i)=CurPar(3);
        
        %Plot curves
%         figure
%         xt=floor(min(xdata)):ceil(max(xdata));
%         yt=ones(1,length(xt))*Acoef(i)-Bcoef(i)*xt.^Ccoef(i);
%         plot(xdata,ydata,'r*')
%         hold all
%         plot(xt,yt)
    end
end

% Create struct
Pcoef = struct('Acoef',Acoef,'Bcoef',Bcoef,'Ccoef',Ccoef);
end


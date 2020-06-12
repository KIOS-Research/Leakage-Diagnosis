function [ K, n ] = Pipe_Coefficients( d, Kepa )
%PIPE_COEFFICIENTS 
formula = d.getBinOptionsInfo.BinOptionsHeadloss;
switch formula
case 'H-W'
    %Hazen-Williams formula
    n=1.852;
    C=double(d.LinkRoughnessCoeff'); %unitless
    D=double(d.LinkDiameter'); %Diameter in mm
    L=double(d.LinkLength'); % Length in meters
    % K = 4.727*(C.^(-n)).*(D.^-4.871).*(L); %Hazen-Williams in imperial units
    K = 1.1088*10^11*(C.^(-n)).*(D.^-4.867).*(L/100); %Hazen-Williams in SI
    K(d.LinkPumpIndex)=-1; %Identify pumps

case 'D-W'
    %Darcy-Weisbach formula
    n=2;
    C=double(d.LinkRoughnessCoeff'); %unitless
    D=double(d.LinkDiameter')/1000; %Diameter in meters
    L=double(d.LinkLength'); % Length in meters
    g=9.8*(12960000); %meters per hour squared
    K = 1.3086*(8*L.*C*(1/60)) ./ (g*pi^2*D.^5); %Darcy-Weisbach in SI
    K =        (8*L.*C)        ./ (g*pi^2*D.^5); %Darcy-Weisbach in SI

%%Cheat   
%     hl=double(d.getLinkHeadloss);
%     flow=double(d.getLinkFlows);
%     Kepa= (hl'./(abs(flow).^n));
    H=diag(K);
    z=Kepa;
    x=pinv(H)*z;
%     K = x.*((8*L.*C) ./ (g*pi^2*D.^5)); 
    K_KEPA=[K  Kepa]
    
otherwise
   error('Unknown headloss formula type given.') 
end

end


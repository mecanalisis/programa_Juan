function [mu2,density1,Cv,v2]=visCal(v,density,t,ti)
% This code is based on Seeton Fit approximation 
% it is valid for a large temperature range 
% t and ti are in C, t is the reference temperature, ti is the decire
% temperature 
% v and v2 are in cSt, v is the oil viscosity at t(1) and t(2) in C, v2 is oil
% viscosity at ti
% mu2 is the dynamic viscosity in Pascal-second Pa-s or kg/m/s

Cv=2000;        %J/Kg/K
beta=.633e-3;   %density/temp coeficient for most mineral oils
ti=ti+273;
to=15.5+273;    %15.5C=60F standard API temp. for SG

if isempty(t)
    t=[40 100]+273; %standard ASTM temps.
else
    t=t+273;
end

V=[log(log(v(1)+.7)); 
   log(log(v(2)+.7))];
M=[1 log(t(1));
   1 log(t(2))];
X=M\V;

density1=density/(1+beta*(ti-to));  %oil volumetric termal expantion

v2=-.7+exp(exp(X'*[1;log(ti)]));    %cSt
mu2=v2*density1*1e-6;               %dynamic viscosity in Pa-s [Kg/m/s]

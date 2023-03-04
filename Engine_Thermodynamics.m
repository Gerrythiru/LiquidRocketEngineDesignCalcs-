function [MachArray,Axarray,PgArray,pRatioArray,Tx,Taw] = Engine_Thermodynamics(gamma,TCcurve,Lc,At,Pc_ns,Tc_ns)
%clear; clc;
%Thermodynamics for the BLT Regen Engine 
%This code plots the Mach number, temperature, and pressure across our 
%thrust chamber. Reference each individual function for more details.  
%
% Source for equations: Nasa Isentropic Flow Equations

% Mach number calculation module
term1 = (gamma+1)/2;
term2 = (gamma+1)/(2*(gamma-1));
term3 = (gamma-1)/2;

%TCcurve(:,1) --> this plots all the x coordinates across the thrust chamber
%TCcurve(:,2); --> this plots all the y coordinates across the thrust chamber

l=length(TCcurve(:,2));
nt = find(TCcurve(:,1) == Lc);           % index of throat
% disp(nt)
Axarray = []; % array of areas for each x coordinate
for k = 1:nt-1

    Axarray(k) = pi * TCcurve(k,2)^2;
    f = @(Mx) ((1/Mx) * (((2/(gamma + 1))*( 1 + (((gamma-1)/2)*(Mx^2)) ))^((gamma+1)/(2*(gamma - 1))))) - (Axarray(k)/At);
    MachArray(k) = fzero(f, [0.001,1.001]);

end

Axarray(nt) = At;
MachArray(nt) = 1;

for k = nt+1:l
    Axarray(k) = pi * TCcurve(k,2)^2;
    f = @(Mx) ((1/Mx) * (((2/(gamma + 1))*( 1 + (((gamma-1)/2)*(Mx^2)) ))^((gamma+1)/(2*(gamma - 1))))) - (Axarray(k)/At);
    MachArray(k) = fzero(f, [1,15]);
    
end

% Pressure calculation 
PgArray = Pc_ns * (1 + term3*MachArray.^2).^(-gamma/(gamma-1));
pRatioArray = Pc_ns ./ PgArray;

%Temperature calculation 
Pr = 4*gamma/(9*gamma - 5);
r = Pr^0.33; 
Tx = Tc_ns ./ (1 + term3*MachArray.^2);
Taw = Tx.*((1+r*((gamma-1)/2)*MachArray.^2)/(1+((gamma-1)/2)*MachArray.^2));


end
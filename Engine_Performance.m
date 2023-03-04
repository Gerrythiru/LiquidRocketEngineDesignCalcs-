function [c_star,cf,E,Ve,mdot,mdotO,mdotF,At,Isp,Engine_efficiency,R] = Engine_Performance(Thrust,Pc_ns,mix_ratio,Pe,Pa,Tc_ns,gamma,MW,Reaction_efficiency,nozzle_efficiency)
% Engine_Performance calculates the expected performance of the engine 
% 
% OUTPUTS:
%       c_star = characteristic velocity (ft/s)
%       cf = thrust coefficient
%       E = expansion ratio
%       Ve = exit velocity
%       mdot = total mass flow rate
%       mdotO = oxidizer mass flow rate
%       mdotF = fuel mass flow rate
%       At = throat area
%       Isp = specific impulse (s)
%       Engine_efficiency = overall engine efficiency
%       R = gas constant
%
% INPUTS:
%       Thrust = engine thrust (lbf)
%       Pc_ns = stagnation pressure (chamber pressure)
%       mix_ratio = O/F ratio
%       Pe = design exit pressure (psi)
%       Pa = ambient pressure (psi)
%       Tc_ns = stagnation temperature (Rankine)
%       gamma = specific heat ratio
%       MW = molecular weight
%       Reaction_efficiency = c_star efficiency
%       nozzle_efficiency = cf efficiency
%       thetaE = nozzle exit angle (deg)



g = 32.17405; %acceleration due to gravity in english units

%terms 1-5 are configurations of the way gamma shows up in the equations
term1 = 2/(gamma+1);
term2 = (gamma+1)/(gamma-1);
term3 = (gamma-1)/gamma;
term4 = 1/gamma;
term5 = (2*gamma^2)/(gamma-1);
term6 = 1/(gamma-1);

R = 1545.349/MW;         %R is a conversion factor to get molecular weight into english units

% characteristic velocity, bluebook Eq. 1-32a 
c_star = sqrt(g*gamma*R*Tc_ns)/(gamma*sqrt((term1)^(term2)))*Reaction_efficiency/100; 

% expansion ratio epsilon, bluebook Eq. 1-20
E = ((term1)^(term6)*((Pc_ns/Pe)^(term4)))/sqrt((term2)*(1-(Pe/Pc_ns)^(term3)));  %expansion ratio

% Thrust coefficient, bluebook Eq. 1-33a
cf = (sqrt((term5)*(term1)^(term2)*(1-(Pe/Pc_ns)^(term3)))+E*((Pe-Pa)/Pc_ns))*nozzle_efficiency/100;  % thrust coefficient

Engine_efficiency = Reaction_efficiency*nozzle_efficiency/100;

% bluebook Eq. 1-34
At = Thrust/(cf*Pc_ns);  %throat area
Ae = At * E;

Ve = c_star * cf;
% bluebook Eq. 1-5
mdot = g*(Thrust - (Pe - Pa)*Ae)/Ve; % mass flow rates
mdotF = mdot/(1+mix_ratio);
mdotO = mix_ratio*mdotF;
Isp = Thrust/mdot;

end
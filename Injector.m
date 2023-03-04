function [Vi, Ai] = Injector(mdot,rho,Cd,dp)
%INJECTOR
%
% OUTPUTS: 
%       Vi = propellant injection velocity (ft/s)
%       Ai = propellant orifice area (in^2)
%
% INPUTS:
%       mdot = propellant mass flow rate
%       rho = prop density
%       Cd = orifice discharge coefficient
%       dp = delta-P (pressure drop)


v0 = 0; % initial velocity
% get injection velocity using bernoulli's eqn
Vi = sqrt((2*dp*1728/rho) + v0^2);

% RPA eq. 8-2 and ModEng eq 4-43
Ai = mdot/(Cd*sqrt(rho*dp/2.238));

end
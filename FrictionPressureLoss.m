function [Velocity, ReK, ReD, FF, FrictionLoss] = FrictionPressureLoss(Density, muK, muD, mdot, ID,...
    RelRough, tubelength)
%This function calculates the pressure drop in the system given its
%dimensions and flow conditions using the Darcy-Weisbach Equation
%Inputs
%density = density of fluid
%muK = kinematic Viscosity
%muD = dynamic Viscosity
%mdot = mass flow rate
%ID = inner diameter of tube
%RelRough = relative roughness of tube
%tubelength = length of tubing

g = 32.17405; %ft/s^2

%Velocity
Velocity = mdot/(Density*(pi/4*ID^2)); %ft/s

%Reynolds Numbers
ReK = Velocity*ID/muK;
ReD = Density/g*Velocity*ID/muD;

%Friction Factor
LaminarFF = 64/ReK; %simple, assuming laminar flow
HaalandFF = (1/(-1.8*log((RelRough/3.7)^1.11+(6.9/ReK))))^2; %complex, assuming trubulent flow, using haaland equation
ColebrookFF = ColebrookEquationBisection(ReK, RelRough);

FF =  ColebrookFF; %set this to the FF you want to use

%Darcy-Weisbach Equation
FrictionLoss = (FF*Density/g*tubelength/ID/2*Velocity^2)/144;
end
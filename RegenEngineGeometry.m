%% Engine Contour
clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%Engine Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<<<<<<< HEAD
Thrust = 1450;  % lbf
Pc_ns=300;      % chamber pressure (psia)
mix_ratio=1.5;  % O/F ratio (a.k.a. mixture ratio)
=======
Thrust = 2000;  % lbf
Pc_ns=350;      % chamber pressure (psia)
mix_ratio=1.6;  % O/F ratio (a.k.a. mixture ratio)
>>>>>>> e591050a14fe50acb0b8dc6065db6c527c832427
Pe = 12;        % exit pressure (psia)
Pa = 13.6;      % atmospheric pressure (psia)
g = 32.17405;   % acceleration due to gravity in english units 

% Numbers specified by designer
Reaction_efficiency = 85;    % predefined reaction efficiency
nozzle_efficiency = 97;      % Cf efficiency
Ec = 4.35;                      % contraction ratio (Ac/At)
L_star = 40;                 % L* (in.)
Lf = 80;                     % Length of bell nozzle compared to conical nozzle (Le/c15)
exportcheck = 1;             % if set to 1, excel file is created
RCcheck = 1;                 % if set to 1, Regen_Channels is run

%Materials -  whichever material is set to 1 will be used.  Only label one
%of the values below as 1.  The rest should remain as a 0.
<<<<<<< HEAD
Cu = 0;                     
In718 = 0;
In625 = 0;
SS = 1;

% Cooling channel parameters
ts = 0.9/25.4;
N = 90;
Rw = 1.3/25.4;
Ch = 1.3/25.4;
=======
matselect = 'Inconel 625';

% Cooling channel parameters
ts = 0.5/25.4;
N = 100;
Rw = 1.2/25.4;
Ch = 0.9/25.4;
>>>>>>> e591050a14fe50acb0b8dc6065db6c527c832427

% Number of divisions per inch for thrust chamber.  The more divisions you
% have the longer the code takes but the more accurate results you get.
% Less divisions runs the code faster, but may not be as accurate.
dx = 100;

% Solid deposit thermal resistance, soot forming on chamber walls lowers thermal conductivity, changes with area ratio
% Used to find overall gas side thermal conductance
% highest thermal conductivity with zero solid deposit
Rd = 0;

%%%%%%%%%%%%%%%%%%%End of Engine Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numbers from RPA given by the engine parameters (make sure to change
% these when changing the O/F ratio)
gamma = 1.2287;     % gamma = isentropic exponent = specific heat rato Cp/Cv
Tc_ns = 4852.9405;  % flame temperature in Rankine
MW = 18.3307;       % molecular weight at the chamber


% Run Engine Peformance module
[c_star,cf,E,Ve,mdot,mdotO,mdotF,At,Isp,Engine_efficiency,R]=Engine_Performance(Thrust,Pc_ns,mix_ratio,Pe,Pa,Tc_ns,gamma,MW,Reaction_efficiency,nozzle_efficiency);

% Run Injector module
%[VO,AO] = Injector(mdotO,rhoO,CdO,dpO);
%[VF,AF] = Injector(mdotF,rhoF,CdF,dpF);

% Run Engine Geometry module
[Rt,Dt,Ae,Me,Re,Ac,Rc,L_cyl_barrel,Lc,Vc,Rcc,Rcd,Rct,thetaN,Ltotal,TCcurve,x_1,x_2,x_3,a,b,c,alpha,thetaE] = Engine_Geometry(dx,At,E,Ec,L_star,Lf,Ve,gamma,exportcheck,Tc_ns,Pc_ns,Pe,R,Engine_efficiency); 

% Run Thermodynamics Module 
[MachArray,Axarray,PgArray,pRatioArray,Tx,Taw] = Engine_Thermodynamics(gamma,TCcurve,Lc,At,Pc_ns,Tc_ns);
% Combustion chamber instability
[L1,T1,R1] = Instability(Lc,Rc,g,gamma,R,Tc_ns);

% Run Regen Channels
if RCcheck == 1
[Twg,Twc,Tco,Vco,deltaPtotal,deltaPinj,PPipe,qg,Dh,Cw,hgc,kw,muF,kF,Re_Dh,ff,hc,Pr,maxTwg,maxTwc,maxTco,minCw]=Regen_Channels(dx,TCcurve,PgArray,MW,Rcd,Rct,gamma,Tc_ns,MachArray,Taw,Dt,Pc_ns,g,c_star,At,Rd,mdotF,Lc,ts,N,Rw,Ch,matselect);
% Run Stress calculator
[Smech,Sth,Stotal,k,cte,ModE,v,maxThermal_Stress,maxTotal_Stress] = StressCalculator(TCcurve,Lc,qg,Dh,ts,PPipe,PgArray,Twg,matselect);
end

%tables to store the data (for reference)
Initial_Data = table(Thrust,Pc_ns,mix_ratio,Pe,Pa,Reaction_efficiency,nozzle_efficiency)
RPA_Data = table(gamma,Tc_ns,MW)
Calculated_Data = table(Isp,c_star,E,cf,Ve,mdot,mdotO,mdotF)
EngineGeometry = table(At,Rt,Ae,Re,Ac,Rc,L_cyl_barrel,Lc,Vc,Rcc,Rcd,Rct,Ltotal,thetaN,thetaE)
EngineGeometry2= table(alpha,x_1,x_2,x_3,a,b,c)
InstTable = table(L1,T1,R1)
Data_Outputs = table(maxTwg,maxTwc,maxTco,deltaPtotal,minCw,maxThermal_Stress,maxTotal_Stress)
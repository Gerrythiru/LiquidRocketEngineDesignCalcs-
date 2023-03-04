%Propulsion Pressure Drop and Tank Sizing Code
%--------------------------------------------------------------------------
%Change hardware sizes, Material, propellant properties and system settings
%Gives the speed of each propellant in the tubes, Re numbers, pressure
%losses, and tank sizes
%--------------------------------------------------------------------------

format long;
close all;
clear;
clc;

%Input Values ====
%Tube Dimensions
%Change these to test different sizes
ODinch = 8/16; %outer diameter or dash size of the tube
Thickness = .035; 

IDft = (ODinch-2*Thickness)/12;

%Tank Dimensions
TankOD = 8/12; %ft
TankThickness = 1/8/12; %ft
TankID = TankOD-2*TankThickness;
TankDomeH = 1.5/12; %ft height of dome modeled as an ellipsoid

TankArea = pi/4*TankID^2;

%System Settings
Tburn = 15; %time of the burn in seconds
LOXUllage = .15; %percent of volume that you want to be ullage
RP1Ullage = .1;

LOXmdot = 2.8; %lbm/s requirement from engine
RP1mdot = 1.27; 

%Equivalent Length of Fittings, more may be added but probably not
TeeRunEL8 = 12; %in
TeeRunEL10 = 14.4; %in
TeeRunEL12 = 16.8; %in

%Tube material roughness: Currently Aluminum
AbsRough = .001; %ft
RelRough = AbsRough/IDft;

%Constants
g = 32.17405; %ft/s^2

%Colebrook Equation Function Values
a = 0;
b = 1;
err = .000001;
iter = 100;

%RP-1 ==============================================================================================
%RP-1 Properties
RP1density = 51.1909; %lbf/ft^3
RP1visKinematic = 21.097*10^-6; %ft^2/s
RP1visDynamic = 3.3567*10^-5;  %lbm*s/ft^2 

%RP-1 System Size
RP1height = 16; %in
RP1numTees = 1;
RP1length = (RP1height + TeeRunEL8*RP1numTees)/12; %ft

%RP1 Bend 1
RP1angleB1 = 15; %degrees
RP1radiusB1 = 1.5/12; %ft
RP1numberB1 = 4;
RP1K1 = .5; %see chart

%RP1 Bend 2
RP1angleB2 = 20; %degrees
RP1radiusB2 = 3/12; %ft
RP1numberB2 = 0;
RP1K2 = .5; %see chart

%RP1 Bend 3
RP1angleB3 = 30; %degrees
RP1radiusB3 = 3/12; %ft
RP1numberB3 = 0;
RP1K3 = .8; %see chart

%RP-1 Calculations ---------------------------------------------------------------------------------
%Velocity
RP1velocity = RP1mdot/(RP1density*(pi/4*IDft^2)); %ft/s

%Reynolds Numbers
RP1reK = RP1velocity*IDft/RP1visKinematic;
RP1reD = RP1density/g*RP1velocity*IDft/RP1visDynamic;

%Friction Factor
RP1laminarFF = 64/RP1reK; %simple assuming laminar flow
RP1haalandFF = (1/(-1.8*log((RelRough/3.7)^1.11+(6.9/RP1reK))))^2; %complex assuming trubulent flow, using haaland equation
RP1colebrookFF = ColebrookEquatiuonBisection(a, b,err, iter, RP1reK, RelRough);

RP1FF =  RP1colebrookFF;%set this to the FF you want to use

%Darcy-Weisbach Equation
RP1frictionLosses = (RP1FF*RP1density/g*RP1length/IDft/2*RP1velocity^2)/144;

%Tube Bend Losses
RP1bendLossB1 = RP1numberB1*(RP1FF/2*RP1density/g*RP1velocity^2*pi*RP1radiusB1/IDft*RP1angleB1/180+RP1K1/2*RP1density/g*RP1velocity^2)/144;
RP1bendLossB2 = RP1numberB2*(RP1FF/2*RP1density/g*RP1velocity^2*pi*RP1radiusB2/IDft*RP1angleB2/180+RP1K2/2*RP1density/g*RP1velocity^2)/144;
RP1bendLossB3 = RP1numberB3*(RP1FF/2*RP1density/g*RP1velocity^2*pi*RP1radiusB3/IDft*RP1angleB3/180+RP1K3/2*RP1density/g*RP1velocity^2)/144;

RP1bendLosses = RP1bendLossB1+RP1bendLossB2+RP1bendLossB3;

%Dynamic Losses
RP1tankVelocity = RP1mdot/(RP1density*TankArea);
RP1dynamicLosses = (RP1density/g/2*(RP1velocity^2-RP1tankVelocity^2)-RP1density*RP1height/12)/144;

%Total
RP1totalLosses = RP1frictionLosses+RP1bendLosses+RP1dynamicLosses;

%LOX ===============================================================================================
%LOX Properties
LOXdensity = 71.2303; %lbf/ft^3
LOXvisKinematic = 1.844*10^-6; %ft^2/s
LOXvisDynamic = 4.087*10^-6; %lbm*s/ft^2

%resets length for LOX
ODinch = 10/16; %outer diameter or dash size of the tube
Thickness = .049; 
IDft = (ODinch-2*Thickness)/12;
RelRough = AbsRough/IDft;

%LOX System Size
LOXtubeHeight = 25; %length in in
LOXnumTees = 1;
LOXlength = (LOXtubeHeight + TeeRunEL10*LOXnumTees)/12; %ft

%LOX Bend 1
LOXAngleB1 = 20; %degrees
LOXradiusB1 = 3/12; %ft
LOXnumberB1 = 0;
LOXKB1 = .5; %see chart

%LOX Calculations ----------------------------------------------------------------------------------
%Gets Velocity > Re > Friction Factor > Pressure Drop, Bend Losses, Dynamic changes
%Velocity
LOXvelocity = LOXmdot/(LOXdensity*(pi/4*IDft^2)); %ft/s

%Reynolds Numbers
LOXreK = LOXvelocity*IDft/LOXvisKinematic;
LOXreD = LOXdensity/g*LOXvelocity*IDft/LOXvisDynamic;

%Friction Factor
LOXlaminarFF = 64/LOXreK; %simple assuming laminar flow
LOXhaalandFF = (1/(-1.8*log((RelRough/3.7)^1.11+(6.9/LOXreK))))^2; %complex assuming trubulent flow, using haaland equation
LOXcolebrookFF = ColebrookEquatiuonBisection(a, b, err, iter, LOXreK, RelRough); %actual value using Colebrook-White equation

LOXFF = LOXcolebrookFF; %set this to the FF you want to use
    
%Darcy-Weisbach Equation
LOXfrictionLosses = (LOXFF*LOXdensity/g*LOXlength/2/IDft*(LOXvelocity)^2)/144;

%Tube Bend Losses
LOXbendLossB1 = LOXnumberB1*(LOXFF/2*LOXdensity/g*LOXvelocity^2*pi*LOXradiusB1/IDft*LOXAngleB1/180+LOXKB1/2*LOXdensity/g*LOXvelocity^2)/144;

LOXbendLosses = LOXbendLossB1;

%Dynamic Losses
LOXtankVelocity = LOXmdot/(LOXdensity*TankArea);
LOXdynamicLosses = (LOXdensity/g/2*(LOXvelocity^2-LOXtankVelocity^2)-LOXdensity*LOXtubeHeight/12)/144;

%LOX SECONDARY CALCULATION FOR DIAMATER CHANGE
%resets length for LOX
ODinch = 8/16; %outer diameter or dash size of the tube
Thickness = .035; 
IDft = (ODinch-2*Thickness)/12;
RelRough = AbsRough/IDft;
LOXvelocity = LOXmdot/(LOXdensity*(pi/4*IDft^2)); %ft/s
LOXreK = LOXvelocity*IDft/LOXvisKinematic;
LOXtubeHeight2 = 8;
LOXlength2 = (LOXtubeHeight2)/12; %ft
LOXFF2 = ColebrookEquatiuonBisection(a, b, err, iter, LOXreK, RelRough); %actual value using Colebrook-White equation
LOXfrictionLosses2 = (LOXFF2*LOXdensity/g*LOXlength2/2/IDft*(LOXvelocity)^2)/144;

%Total
LOXtotalLosses = LOXfrictionLosses+LOXbendLosses+LOXdynamicLosses+LOXfrictionLosses2;

%Tank Sizing =======================================================================================
LOXDCod = 1; %in tube inside lox tank %tube + 1/4 + Tank wall thickness
LOXDowncomerArea = pi/4*(LOXDCod/12)^2;       %ft^2

RP1DCod = 1.125; %in tube inside RP1 tank
RP1DowncomerArea = pi/4*(RP1DCod/12)^2;       %ft^2

DomeVolume = pi/6*TankID^2*TankDomeH; %ft^3

%LOX calculations
LOXmass = LOXmdot*Tburn; %lbm
LOXvolume = LOXmass/(LOXdensity); %ft^3
LOXheight = (LOXvolume-(DomeVolume-LOXDowncomerArea*TankDomeH))/(TankArea-LOXDowncomerArea); %ft
LOXullageVolume = LOXvolume*LOXUllage; %ft^3
LOXullageHeight = (LOXullageVolume-(DomeVolume-LOXDowncomerArea*TankDomeH))/(TankArea-LOXDowncomerArea);%ft
LOXcylinderHeight = (LOXheight+LOXullageHeight)*12; %in
LOXtankHeight = LOXcylinderHeight+2*TankDomeH*12; %in total height including domes

%RP-1 calculations
RP1mass = RP1mdot*Tburn; %lbm
RP1volume = RP1mass/(RP1density);%ft^3
RP1height = (RP1volume-(DomeVolume-RP1DowncomerArea*TankDomeH))/(TankArea-RP1DowncomerArea);%ft
RP1ullageVolume = RP1volume*RP1Ullage;%ft^3
RP1ullageHeight = (RP1ullageVolume-(DomeVolume-RP1DowncomerArea*TankDomeH))/(TankArea-RP1DowncomerArea);%ft
RP1cylinderHeight = (RP1height + RP1ullageHeight)*12;%in
RP1tankHeight = RP1cylinderHeight + 2*TankDomeH*12;%in total height including domes

%COPV Sizing

%tank pressures based on chamber pressure + engine pressure drop + prop pressure drop
ChambPress = 450; %PSI
EngineDelP = 90;
LOXpressure = 600;%ChambPress + EngineDelP + LOXtotalLosses; %PSI
RP1pressure = 600;%ChambPress + EngineDelP + RP1totalLosses;

CCF = 2.8; %Cryo collapse factor, no clue what this should be

%COPVpressure = 5000; %PSI 
%COPVvolume = (LOXpressure*LOXvolume*CCF+RP1pressure*RP1volume)/COPVpressure*28.317; %Liters

COPVvolume = 6.8; %Liters
COPVpressure = (LOXpressure*LOXvolume*CCF+RP1pressure*RP1volume)/(COPVvolume/28.317); %PSI

%TO DO ---------------------------------------------------------------------------------------------
%Helium Numbers
%Cryo Collapse

%Prints out info for each fluid ====================================================================
fprintf('    Properties      | LOX Numbers   | RP-1 Numbers\n');
fprintf('--------------------|---------------|----------------\n');
fprintf('Velocity\t\t\t| %4.2f ft/s\t| %.2f ft/s\n', LOXvelocity, RP1velocity);
fprintf('Re \t\t\t\t\t|%7.0f\t\t| %.0f \n', LOXreK, RP1reK);
%fprintf('Re \t\t\t| %7.0f\t\t| %.0f \n', LOXreD, RP1reD);
fprintf('Friction Factor\t\t| %0.6f\t\t| %.6f \n', LOXFF, RP1FF);
fprintf('Friction Losses\t\t| %4.2f PSI  \t| %.2f PSI\n', LOXfrictionLosses, RP1frictionLosses);
fprintf('Tube Bend Losses\t| %4.2f PSI\t\t| %.2f PSI\n', LOXbendLosses, RP1bendLosses);
fprintf('Dynamic Losses\t\t| %4.2f PSI\t\t| %.2f PSI\n', LOXdynamicLosses, RP1dynamicLosses);
fprintf('Total Pressure Drop | %4.2f PSI  \t| %.2f PSI\n', LOXtotalLosses, RP1totalLosses);
fprintf('--------------------|---------------|----------------\n');
fprintf('Propellant Mass\t\t| %3.3f lb   \t| %.3f lb\n', LOXmass, RP1mass);
fprintf('Propellant Volume\t| %1.3f ft^3\t| %.3f ft^3\n', LOXvolume, RP1volume);
%fprintf('Propellant Height\t| %1.3f ft\t\t| %.3f ft\n', LOXheight, RP1height);
fprintf('Diptube Length\t\t| %1.3f in\t\t| %.3f in\n', (LOXullageHeight+TankDomeH)*12, (RP1ullageHeight+TankDomeH)*12);
fprintf('Cylinder Height\t\t| %1.3f in\t\t| %.3f in\n', LOXcylinderHeight, RP1cylinderHeight);
fprintf('Tank Height\t\t\t| %1.3f in\t\t| %.3f in\n', LOXtankHeight, RP1tankHeight);
fprintf('--------------------|---------------|----------------\n');
fprintf('COPV Pressure is %.3f PSI for a %.3f liter tank\n', COPVpressure, COPVvolume);
fprintf('Required LOX Pressure is %.2f\n', LOXpressure);
fprintf('Required RP1 Pressure is %.2f\n', RP1pressure);

%SV Cv test stuff
sg = RP1density/62.4;
q = RP1mdot/RP1density*7.481;
delp = sg/(.06^2)*q^2;


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

%Input Values ======================================================================================
%Tube Dimensions
%Change these to test different sizes
ODinch = 8/16; %outer diameter or dash size of the tube
Thickness = .035; 

IDft = (ODinch-2*Thickness)/12;

%Tank Dimensions
%TankOD = 8/12; %ft
%TankThickness = 1/8/12; %ft
TankID = 7.981/12; %TankOD-2*TankThickness; %ft
TankDomeH = 1.5/12; %ft height of dome modeled as an ellipsoid

TankArea = pi/4*TankID^2; %ft^2

%System Settings
Tburn = 20; %time of the burn in seconds
LOXUllage = .2; %percent of volume that you want to be ullage
CH4Ullage = .2;

LOXmdot = 2.81453; %lbm/s requirement from engine
CH4mdot = .9039; 

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

%CH4 ==============================================================================================
%CH4 Properties
CH4density = 26.38; %lbf/ft^3
CH4visKinematic = 2.985*10^-6; %ft^2/s
CH4visDynamic = 2.4478*10^-5;  %lbm*s/ft^2 

%CH4 System Size
CH4height = 20; %in
CH4numTees = 1;
CH4length = (CH4height + TeeRunEL8*CH4numTees)/12; %ft

%CH4 Bend 1
CH4angleB1 = 15; %degrees
CH4radiusB1 = 1.5/12; %ft
CH4numberB1 = 4;
CH4K1 = .5; %see chart

%RP-1 Calculations ---------------------------------------------------------------------------------
%Velocity
CH4velocity = CH4mdot/(CH4density*(pi/4*IDft^2)); %ft/s

%Reynolds Numbers
CH4reK = CH4velocity*IDft/CH4visKinematic;
CH4reD = CH4density/g*CH4velocity*IDft/CH4visDynamic;

%Friction Factor
CH4laminarFF = 64/CH4reK; %simple assuming laminar flow
CH4haalandFF = (1/(-1.8*log((RelRough/3.7)^1.11+(6.9/CH4reK))))^2; %complex assuming trubulent flow, using haaland equation
CH4colebrookFF = ColebrookEquatiuonBisection(a, b,err, iter, CH4reK, RelRough);

CH4FF =  CH4colebrookFF;%set this to the FF you want to use

%Darcy-Weisbach Equation
CH4frictionLosses = (CH4FF*CH4density/g*CH4length/IDft/2*CH4velocity^2)/144;

%Tube Bend Losses
CH4bendLosses = CH4numberB1*(CH4FF/2*CH4density/g*CH4velocity^2*pi*CH4radiusB1/IDft*CH4angleB1/180+CH4K1/2*CH4density/g*CH4velocity^2)/144;

%Dynamic Losses
CH4tankVelocity = CH4mdot/(CH4density*TankArea);
CH4dynamicLosses = (CH4density/g/2*(CH4velocity^2-CH4tankVelocity^2)-CH4density*CH4height/12)/144;

%Total
CH4totalLosses = CH4frictionLosses+CH4bendLosses+CH4dynamicLosses;

%LOX ===============================================================================================
%LOX Properties
LOXdensity = 71.2303; %lbf/ft^3
LOXvisKinematic = 1.844*10^-6; %ft^2/s
LOXvisDynamic = 4.087*10^-6; %lbm*s/ft^2

%LOX System Size
LOXtubeHeight = 11; %length in in
LOXnumTees = 1;
LOXlength = LOXtubeHeight/12; %ft

%CH4 Bend 1
LOXangleB1 = 0; %degrees
LOXradiusB1 = 1.5/12; %ft
LOXnumberB1 = 0;
LOXK1 = .5; %see chart

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
LOXbendLosses = LOXnumberB1*(LOXFF/2*LOXdensity/g*LOXvelocity^2*pi*LOXradiusB1/IDft*LOXangleB1/180+LOXK1/2*LOXdensity/g*LOXvelocity^2)/144;

%Dynamic Losses
LOXtankVelocity = LOXmdot/(LOXdensity*TankArea);
LOXdynamicLosses = (LOXdensity/g/2*(LOXvelocity^2-LOXtankVelocity^2)-LOXdensity*LOXtubeHeight/12)/144;

%LOX SECONDARY CALCULATION FOR DIAMATER CHANGE
%resets length for LOX
ODinch = 10/16; %outer diameter or dash size of the tube
Thickness = .049; 
IDft = (ODinch-2*Thickness)/12;
RelRough = AbsRough/IDft;
LOXvelocity = LOXmdot/(LOXdensity*(pi/4*IDft^2)); %ft/s
LOXreK = LOXvelocity*IDft/LOXvisKinematic;
LOXtubeHeight2 = 30;
LOXlength2 = (LOXtubeHeight2 + TeeRunEL10*LOXnumTees)/12; %ft
LOXFF2 = ColebrookEquatiuonBisection(a, b, err, iter, LOXreK, RelRough); %actual value using Colebrook-White equation
LOXfrictionLosses2 = (LOXFF2*LOXdensity/g*LOXlength2/2/IDft*(LOXvelocity)^2)/144;

%Total
LOXtotalLosses = LOXfrictionLosses+LOXbendLosses+LOXdynamicLosses+LOXfrictionLosses2;

%Tank Sizing =======================================================================================
LOXDCod = .5; %in, tube inside lox tank 
LOXDowncomerArea = pi/4*(LOXDCod/12)^2;       %ft^2

CH4DCod = 10/16; %in, tube inside CH4 tank
CH4DowncomerArea = pi/4*(CH4DCod/12)^2;       %ft^2

DomeVolume = pi/6*TankID^2*TankDomeH; %ft^3
BoilOffFraction = 1; %multiplyer to account for boil off of propellant loaded first

%only works if ullage volume is greater than dome volume to avoid complex modeling of dome volumes
%LOX calculations
LOXmass = LOXmdot*Tburn*BoilOffFraction; %lbm
LOXvolume = LOXmass/(LOXdensity); %ft^3
LOXullageVolume = LOXvolume*LOXUllage; %ft^3
LOXheight = (LOXvolume)/(TankArea-LOXDowncomerArea); %ft
LOXullageHeight = (LOXullageVolume-(DomeVolume-LOXDowncomerArea*TankDomeH))/(TankArea-LOXDowncomerArea); %ft
LOXcylinderHeight = (LOXheight+LOXullageHeight)*12; %in
LOXtankHeight = LOXcylinderHeight+TankDomeH*12; %in, total height including domes

%CH4 calculations
CH4mass = CH4mdot*Tburn; %lbm
CH4volume = CH4mass/(CH4density); %ft^3
CH4ullageVolume = CH4volume*CH4Ullage; %ft^3
CH4height = (CH4volume-(DomeVolume-CH4DowncomerArea*TankDomeH))/(TankArea-CH4DowncomerArea); %ft
CH4ullageHeight = (CH4ullageVolume)/(TankArea-CH4DowncomerArea); %ft
CH4cylinderHeight = (CH4height + CH4ullageHeight)*12; %in
CH4tankHeight = CH4cylinderHeight + TankDomeH*12; %in, total height including domes

%COPV Sizing

%tank pressures based on chamber pressure + engine pressure drop + prop pressure drop
ChambPress = 450; %PSI
EngineDelP = 90;
LOXpressure = 600;%ChambPress + EngineDelP + LOXtotalLosses; %PSI
CH4pressure = 600;%ChambPress + EngineDelP + RP1totalLosses;

CCF = 1.7; %Cryo collapse factor, no clue what this should be

%COPVpressure = 5000; %PSI 
%COPVvolume = (LOXpressure*LOXvolume*CCF+RP1pressure*RP1volume)/COPVpressure*28.317; %Liters

COPVvolume = 6.8; %Liters
COPVpressure = (LOXpressure*LOXvolume*CCF+CH4pressure*CH4volume)/(COPVvolume/28.317); %PSI

%TO DO ---------------------------------------------------------------------------------------------
%Helium Numbers
%Cryo Collapse

%Prints out info for each fluid ====================================================================
fprintf('    Properties      | LOX Numbers   | CH4 Numbers\n');
fprintf('--------------------|---------------|----------------\n');
fprintf('Velocity\t\t\t| %4.2f ft/s\t| %.2f ft/s\n', LOXvelocity, CH4velocity);
fprintf('Re \t\t\t\t\t|%7.0f\t\t| %.0f \n', LOXreK, CH4reK);
%fprintf('Re \t\t\t| %7.0f\t\t| %.0f \n', LOXreD, RP1reD);
fprintf('Friction Factor\t\t| %0.6f\t\t| %.6f \n', LOXFF, CH4FF);
fprintf('Friction Losses\t\t| %4.2f PSI  \t| %.2f PSI\n', LOXfrictionLosses+LOXfrictionLosses2, CH4frictionLosses);
fprintf('Tube Bend Losses\t| %4.2f PSI\t\t| %.2f PSI\n', LOXbendLosses, CH4bendLosses);
fprintf('Dynamic Losses\t\t| %4.2f PSI\t\t| %.2f PSI\n', LOXdynamicLosses, CH4dynamicLosses);
fprintf('Total Pressure Drop | %4.2f PSI  \t| %.2f PSI\n', LOXtotalLosses, CH4totalLosses);
fprintf('--------------------|---------------|----------------\n');
fprintf('    Properties      | LOX Numbers   | CH4 Numbers\n');
fprintf('--------------------|---------------|----------------\n');
fprintf('Propellant Mass\t\t| %3.3f lb   \t| %.3f lb\n', LOXmass, CH4mass);
fprintf('Propellant Volume\t| %1.3f ft^3\t| %.3f ft^3\n', LOXvolume, CH4volume);
%fprintf('Propellant Height\t| %1.3f ft\t\t| %.3f ft\n', LOXheight, RP1height);
fprintf('Diptube Length\t\t| %1.3f in\t\t| %.3f in\n', (LOXullageHeight+TankDomeH)*12, (CH4ullageHeight)*12);
fprintf('Cylinder Height\t\t| %1.3f in\t\t| %.3f in\n', LOXcylinderHeight, CH4cylinderHeight);
fprintf('Tank Height\t\t\t| %1.3f in\t\t| %.3f in\n', LOXtankHeight, CH4tankHeight);
fprintf('--------------------|---------------|----------------\n');
fprintf('COPV Pressure: %.3f PSI \nCOPV Size %.3f liters \nBurn Time: %.1f s\n', COPVpressure, COPVvolume, Tburn);
%fprintf('Required LOX Pressure is %.2f\n', LOXpressure);
%fprintf('Required CH4 Pressure is %.2f\n', CH4pressure);


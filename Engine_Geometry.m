function [Rt,Dt,Ae,Me,Re,Ac,Rc,L_cyl_barrel,Lc,Vc,Rcc,Rcd,Rct,thetaN,Ltotal,TCcurve,x_1,x_2,x_3,a,b,c,alpha,thetaE] = Engine_Geometry(dx,At,E,Ec,L_star,Lf,Ve,gamma,exportcheck,Tc_ns,Pc_ns,Pe,R,Engine_efficiency)
% Engine_Geometry returns thrust chamber geometry
%   
%       OUTPUTS:
%       Rt = radius at throat (in)
%       Ae = exit area (in^2)
%       Me = exit mach number
%       Re = radius at exit (in)
%       Ac = contraction area (cylindrical chamber area) (in^2)
%       Rc = chamber radius (in)
%       L_cyl_barrel = length of cylindrical barrel (in)
%       Lc = chamber length or axial d
%       Vc = chamber volume
%       Rcc = first converging curve radius
%       Rcd = converging throat radiusistance to throat (length of cylindrical barrel + converging section)
%       Rct = diverging throat radius
%       thetaN = diverging nozzle angle, calculated using exit mach number (deg)
%       Lfcalc = length fraction of the generated curve (%)
%       Ltotal = total length of thrust chamber + nozzle
%       TCcurve = full coordinates of the thrust chamber curve
%       x_1 = The x distance of Section II
%       x_2 = The x distance for Section III
%       x_3 = The x distance for Section IV
%       a = The leading coefficient of the parabolic bell nozzle (the
%           nozzle follows the form of x = ay^2+by+c)
%       b = The coefficient of the second term of the parabolic bell nozzle 
%           (see above equation)
%       c = The coefficient of the third term of the parabolic nozzle (see
%           above equation)
%       alpha = converging half angle (values range from 20 - 45 deg)
%
%       INPUTS:
%       At = area at throat (in^2)
%       E = expansion ratio (Ae/At)
%       Ec = contraction ratio (Ac/At)
%       Lf = length fraction of nozzle (only used for cubic nozzle)
%       Ve = exit velocity (ft/s)
%       gamma = specific heat ratio
%       thetaE = nozzle exit angle (degrees)
%       exportcheck - if set to 1, the excel file is created
%       curveSelect - if set to 1, parabolic curve is created, if set to 0, cubic curve is created
%       Tc_ns = stagnation temperature (initial temp)
%       Pc_ns = stagnation pressure (chamber pressure)
%       Pe = exit pressure
%       R = gas constant
%       Engine_efficiency = overall engine efficiency


g = 32.17405;
Dt = sqrt(4*At/pi);
Rt = Dt/2;

Ae = At*E;           % nozzle exit area
De = sqrt(4*Ae/pi);  % diameter at exit
Re = De/2;           % radius at exit
            
Ac = Ec*At;          % contraction area
Dc = sqrt(4*Ac/pi);  % chamber diameter
Rc = Dc/2;           % radius at chamber


%%%%Engine Geometry Parameter 

term1 = gamma+1;
term2 = gamma-1;

% Mach no at exit
Te = Tc_ns/((Pc_ns/Pe)^(term2/gamma));      % exit temperature Eq. 1-13
ae = sqrt(g*gamma*R*Te);                    % speed of sound at exit (bluebook page 12)
Me = Ve/(ae*Engine_efficiency/100);         % mach no at exit (bluebook page 12)

% Calculate diverging nozzle half angle (thetaN) using Prandtl-Meyer/2
% Prandtl-Meyer function gives angle of expansion fan starting at Mach 1
PMang = sqrt(term1/term2)*atan(sqrt((term2/term1)*(Me^2 - 1))) - atan(sqrt(Me^2 - 1));
thetaN = (PMang*180/pi)/2;


alpha = 30;                %converging half angle (values range from 20 - 45 deg)
Vc = L_star*At ;           %chamber volume
Rcd = 1.5*Rt;
x_3 = Rcd*sind(alpha);
%Assume Rcc/Dc = 1
Rcc = Dc;
x_1 = Rcc*sind(alpha);
x_2 = (Dc/2-(Rt+Rcc*(1-cosd(alpha))+Rcd*(1-cosd(alpha))))/tand(alpha);
L_conv = x_1+x_2+x_3;
Vconvergent_cone = (pi/3)*L_conv*(Rc^2+Rt^2+Rc*Rt);
V_cyl_barrel = Vc-Vconvergent_cone;

L_cyl_barrel = V_cyl_barrel/(pi*Rc^2);
Lc = L_conv+L_cyl_barrel;

%calculations for the diverging nozzle section
Rct = 0.382*Rt;
Lct = Rct*sind(thetaN);
L15 = (Rt*(sqrt(E)-1)+Rct*(secd(15)-1))/tand(15);    % length of 15 degree conical nozzle, bluebook Eq. 4-7
L_bell = (Lf/100)*L15;                               % bell length 
Ltotal = Lc+L_bell;                                  % total length


% Engine Contour for LOX/RP-1 Regen Engine

%The engine contour is split into 6 sections

%Section I
%Straight line that makes the combustion chamber
x1 = linspace(0,L_cyl_barrel,L_cyl_barrel*dx);
y1 = linspace(Rc,Rc,L_cyl_barrel*dx);
plot(x1,y1,'k')
hold on

T1 = [x1', y1']; % array for this section
TCcurve = T1; % add section to thrust chamber curve array

%Section II
%circular arc making the first part of the converging section
x2 = [L_cyl_barrel,L_cyl_barrel+x_1];
xc1 = linspace(L_cyl_barrel,L_cyl_barrel+x_1,x_1*dx);
xtheta1 = linspace(0,x_1,x_1*dx);
theta1 = acos(xtheta1/Rcc);
%theta1 = pi/3:0.01:pi/2; %based on an alpha of 30 degrees
yCent1 = -Rc;
xCoord1 = xc1;
%xCoord1 = L_cyl_barrel + Rcc*cos(theta1);
yCoord1 = yCent1 + Rcc*sin(theta1);
xCoord1(1) = []; xCoord1(end) = []; % delete repeated points
yCoord1(1) = []; yCoord1(end) = [];

plot(xCoord1,yCoord1,'k')
hold on
T2 = [xCoord1', yCoord1']; % store x and y in single array
TCcurve = [TCcurve; T2];   % store this segment in TCcurve

%Section III
%straight line of the converging section
x3 = linspace(L_cyl_barrel+x_1,L_cyl_barrel+x_1+x_2,x_2*dx);
y3 = linspace(Rc-(Rcc-Rcc*sind(90-alpha)),Rc-(Rcc-Rcc*sind(90-alpha))-(x_2*tand(alpha)), x_2*dx);
plot(x3,y3,'k')
hold on

T3 = [x3', y3'];
TCcurve = [TCcurve; T3];

%Section IV
%the throat section part 1
x4 = [L_cyl_barrel+x_1+x_2,Lc];
xthroat1 = linspace(L_cyl_barrel+x_1+x_2,Lc,x_3*dx);
xttheta2 = linspace(0,x_3,x_3*dx);
theta2 = acos(xttheta2/Rcd);
yCent2 = Rt+Rcd;
xCoord2 = xthroat1;
yCoord2 = yCent2 - Rcd*sin(theta2);
yCoord2 = fliplr(yCoord2);
xCoord2(1) = [];
yCoord2(1) = [];
plot(xCoord2,yCoord2,'k')
hold on

T4 = [xCoord2', yCoord2'];
TCcurve = [TCcurve; T4]; 

%Section V
%the throat section part 2
x5 = [Lc,Lc+Lct];
xthroat2 = linspace(Lc,Lc+Lct,Lct*dx);
xtheta3 = linspace(0,Lct,Lct*dx);
theta3 = acos(xtheta3/Rct);
yCent3 = Rt+Rct;
xCoord3 = xthroat2;
yCoord3 = yCent3 - Rct*sin(theta3);
xCoord3(1) = [];
yCoord3(1) = [];
xCoord3(end) = [];
yCoord3(end) = [];
plot(xCoord3,yCoord3,'k')
hold on

T5 = [xCoord3', yCoord3'];
TCcurve = [TCcurve; T5];

%Section VI
% Diverging nozzle curve
xN = Lc+Lct;                  % x coordinate where diverging throat curve ends/ where parabolic curve begins
x6 = [xN,Ltotal];
Rn = Rt+Rct-Rct*cosd(thetaN); % radius where diverging throat curve ends

% Quadratic curve
% using angles and starting point as initial conditions calculate parabolic nozzle
A = [2*Rn 1 0; Re^2 Re 1; Rn^2 Rn 1];  
B = [1/tand(thetaN); Ltotal; xN];
abc = A\B;
a = abc(1);
b = abc(2);
c = abc(3);

x_exit = linspace(xN,Ltotal,L_bell*dx);
y_exit = (-b+sqrt(b^2 - 4*a*(c-x_exit)))/(2*a);
x_div = x_exit;
y_div = y_exit;

% calculate exit angle
thetaE = atand((y_exit(end) - y_exit(end -1))/(x_exit(end) - x_exit(end-1)));

plot(x_div,y_div,'k')
hold off
axis([0, 16, 0, 9])
xlabel('Thrust Chamber Length (in)');
ylabel('Thrust Chamber Height (in)');
title('Thrust Chamber 1-D Cross Section')

T6 = [x_div', y_div'];
TCcurve = [TCcurve; T6];
tccsize = size(TCcurve);
lentcc = tccsize(1);
zeroarray = zeros(lentcc, 1);
TCcurve = [TCcurve zeroarray];    % zero for every z vector, necessary for SolidWorks to be happy

if exportcheck == 1
xlswrite('TCcurve.xlsx', TCcurve);
end


end
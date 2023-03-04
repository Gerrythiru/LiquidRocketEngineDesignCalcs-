function [Smech,Sth,Stotal,k,a,ModE,v,maxThermal_Stress,maxTotal_Stress] = StressCalculator(TCcurve,Lc,q,Dh,t,Pco,PgArray,Tw,matselect)
% Heat stress calculation function
%Sc = Max combined compressive stress, lb/in^2
%q = Heat flux, Btu/in.2-s
%Dh = channel hydraulic diameter
%R = Radius of the inner shell, in
%t = Thickness of the inner wall, in
%Pco = Coolant pressure, lb/in.2
%Pg = Combustion-gas pressure, lb/in.2
%E = Modulus of elasticity of inner-shell material, lb/in.2
%a = Thermal expansion coefficient of innershell material, in./in.-°F
%k = Thermal conductivity of inner-shell material, Btu/in.2-s-°F/in
%v = Poisson's ratio of inner-shell material

% Flip regen channel arrays to do stress calculation starting at inlet
Dh = flipud(Dh);
R = Dh/2;
Tw = fliplr(Tw);
q = fliplr(q);
for i = 1:length(q)
    [kp,ap,Ep,vp] = Materials(matselect);
    k(i) = polyval(kp,Tw(i));
    a(i) = polyval(ap,Tw(i));
    ModE(i) = polyval(Ep,Tw(i));
    v(i) = polyval(vp,Tw(i));
    Smech(i) = (Pco(i)-PgArray(i))*R(i)/t;
    Sth(i) = ModE(i)*a(i)*q(i)*t/(2*(1-v(i))*k(i));
end

Stotal = Sth+Smech;
% material properties 
%figure(7)
%plot(TCcurve(:,1),k)
%figure(8)
%plot(TCcurve(:,1),a)
%figure(9)
%plot(TCcurve(:,1),E)
%figure(10)
%plot(TCcurve(:,1),v)

figure(7)
plot(TCcurve(:,1),Sth,'DisplayName','Thermal Stresses')
xline(Lc,'k','DisplayName','Throat Location')
xlabel('Thrust Chamber Length (in)');
ylabel('Thermal Stress (lb/in^2)');
title('Thermal Stress Across The Thrust Chamber')
lgd = legend('Location','northwest');

figure(8)
plot(TCcurve(:,1),Smech,'DisplayName', 'Mechanical Stresses')
xlabel('Thrust Chamber Length (in)');
ylabel('Mechanical Stress (lb/in^2)');
title('Mechanical Stress Across The Thrust Chamber')
lgd = legend('Location','northwest');

figure(9)
plot(TCcurve(:,1),Stotal,'DisplayName','Mechanical & Thermal Stress')
xline(Lc,'k','DisplayName','Throat Location')
xlabel('Thrust Chamber Length (in)');
ylabel('Total Stress (lb/in^2)');
title('Total Stress Across The Thrust Chamber')
lgd2 = legend('Location','northwest');

maxThermal_Stress = max(Sth);
maxTotal_Stress = max(Stotal);
% minimumk = min(k)
% minimuma = min(a)
% minimumE = min(E)
% minimumv = min(v)
end 
function [L1,T1,R1] = Instability(Lc,Rc,g,gamma,R,Tc)

Dc = Rc * 2; 
% Speed of sound in combustion chamber
ac = sqrt(g*gamma*R*Tc);

% Frequencies of acoustic instabilities
% Source: Bluebook page 129

% First Longitudinal Mode
L1 = 0.5 * ac / Lc;

% First Tangential Mode
T1 = 0.59 * ac / Dc;

% First Radial Mode
R1 = 1.22 * ac / Dc;
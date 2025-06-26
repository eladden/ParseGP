function A = RCS2size(RCS)
%% A = RCS2size(RCS)
%returns the cross section of a debris from a known Radar Cross Section 
%(RCS) calculation is acc. to Gautam D. Badhawar "Determination of the Area 
%and Mass Distribution of Orbital Debris Fragments" 1989.
%this finction is good for debris which size is at least 8cm.
%input:  Radar CrosRCS in m^2
%output: cross section in m^2


A = 0.5712*(RCS^0.7666);
if A < RCS
    A = RCS;
end
    
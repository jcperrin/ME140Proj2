%% Project 2 Analysis
% This script is the main controller for the SR30 analysis in Project 2 of
% ME140. It should import data, and call the component functions for each
% of the individual parts of the project.
%
% Authors: Jean-Christophe Perrin, Beck Goodloe, Richard Randall, Jason
% Trinidad
%
% Created: 2018-04-11
% Edited: 2018-04-16
clear all;
clc;

%% Useful Constants
% These are mostly conversion factors so that we can convert from imperial
% units collected into metric.

const.insqToMsq = 0.00064516; % [m^2/in^2]
const.lbfToN = 4.44822; % [N/lbf]
const.R = 287; % [kJ/kg/k] Gas constant of air
const.KCdiff = 273; % [deg]
const.kPaToPa = 1e3; % [Pa/kPa]
const.kJkmol2Jkg = 1e3/28.97; %[J/kg * kmol/kJ]
const.airCoefs = const.kJkmol2Jkg .* ...
    fliplr([28.11, 0.1967e-2, 0.4802e-5,-1.966e-9]); % [J/kg/k] cp air
const.LHVJetA = 42800; % [kJ/kg]
const.kg2g = 1e3; % [kg/g]

%% Given Data
% The following values were supplied to us in the original project
% specifications and should not be changed.

given.pitotEffectiveArea = 6.4*const.insqToMsq; % [m^2]
given.jetAHeating = 42.8e6; % [J/kg/K]
given.A = [27.3, 6.4, 9.0, 7.2, 4.7, 3.87]'.*const.insqToMsq; % [m^2]

%% Import Collected Data
% We import the collected data into an table for use throughout our
% function. This table should never be changed, compromising the validity
% of our collected data. Instead, all deried values should be copied out
% the table into a seperate vector or table.

fname = 'data.txt';
collectedData = readtable(fname);
collectedData.Properties.VariableUnits = {'', 'C', 'C', 'C', 'C', 'C', ...
    'C', 'kPa', 'kPa', 'kPa', 'kPa', 'kPa', 'kg/s', 'lbs'};

nObservations = height(collectedData); % Number of observations
% disp(collectedData);

%%
% We also measured the following values seperately.
T1 = 23; % [C] room temp
P1 = 101.4e3; % [Pa] atmospheric pressue

%% Solve for Each RPM Observation
% Loop over every observation taken. For each set of Tm or Pm (the measured
% values of Temp and Pressure) call Richard's Analysis to calculated the
% values of interest (T0, P0, M, V). These are stored in a 2D array. Each
% row of the array corresponds to an observation (i.e. one RPM). 
% Each column is a location of interest.

nLocations = 6; % locations of interest 1-5, 8
T0 = NaN(nObservations, nLocations);
P0 = NaN(nObservations, nLocations);
Ma = NaN(nObservations, nLocations);
V = NaN(nObservations, nLocations);

mdotAir = NaN(nObservations, 1); % [kg/s]

for iObservation = 1:nObservations
    Tm = collectedData{iObservation, 2:6};
    Pm = collectedData{iObservation, 8:12};
    mdotFuel = collectedData{iObservation, 13};
    thrust = collectedData{iObservation, 14};
    [calculatedValues, thisMdotAir] = solveEachLocation(Tm, Pm, mdotFuel, thrust);
    calculatedValues([6, 7], :) = [];
    mdotAir(iObservation) = thisMdotAir;
    T0(iObservation, :) = cell2mat(calculatedValues.T0)';
    P0(iObservation, :) = cell2mat(calculatedValues.P0)';
    Ma(iObservation, :) = cell2mat(calculatedValues.M)';
    V(iObservation, :) = cell2mat(calculatedValues.V)';
end

%% Deliverable 2
% Use your performance data and area measurements (listed in Appendix) to
% construct a series of plots showing how the following quantities vary
% with spool speed: 
%
% * station stagnation temperature (K), 
% * station stagnation pressure (kPa, absolute), 
% * station Mach number
% * station velocity (m/s). 
%
% Make sure to include Station 1 in all of your plots.
%
% First will calculate the deisred values, then plot them later on.

krpm = collectedData.RPM/1000; % rescale RPM for plotting

mdotFuel = collectedData.FuelFlow; % [kg/s]
airFuelRatio = mdotAir./mdotFuel; % [1]

%% Stagnation Temp vs Spool Speed
legendString = {'Station 1', 'Station 2', 'Station 3', 'Station 4', ...
    'Station 5', 'Station 8'};

plot(krpm, T0, '-o');
xlabel('Spool speed [kRPM]');
ylabel('Stagnation Temp [k]');
legend(legendString, 'Location', 'bestoutside');
title('Stagnation Temperature v. Spool Speed');
plotFixer();
print('-depsc','-tiff','-r300','plots/stagTVsRpm');

%% Stagnation Pressure vs Spool Speed
plot(krpm, 1e-3.*P0, '-o');
xlabel('Spool speed [kRPM]');
ylabel('Stagnation Pressure [kPa]');
legend(legendString, 'Location', 'bestoutside');
title('Stagnation Pressure v. Spool Speed');
plotFixer();
print('-depsc','-tiff','-r300','plots/stagPVsRpm');

%% Mach Number vs Spool Speed
plot(krpm, Ma, '-o');
xlabel('Spool speed [kRPM]');
ylabel('Mach Number [1]');
legend(legendString, 'Location', 'bestoutside');
title('Mach Number v. Spool Speed');
plotFixer();
print('-depsc','-tiff','-r300','plots/MaVsRpm');

%% Velocity vs Spool Speed
plot(krpm, V, '-o');
xlabel('Spool speed [kRPM]');
ylabel('Velocity [m/s]');
legend(legendString, 'Location', 'bestoutside');
title('Velocity v. Spool Speed');
plotFixer();
print('-depsc','-tiff','-r300','plots/velVsRpm');

%% Air mass flow rate, fuel mass flow rate, air-fuel ratio
% All three of these quantities to be plotted on the same plot.
% TA (Isabel) recommends plotting mdot_air/10 and mdot_fuel*10 so that the
% plot is meaningful, and indicating the change in the legend
%
% TODO: Convert to plotyy

scaledMdotAir = mdotAir.*const.kg2g; % [g/s]
scaledMdotFuel = mdotFuel.*const.kg2g * 100; % 100 [g/s]

close;
figure;
[hAx,hLine1,hLine2] = plotyy(krpm, airFuelRatio, [krpm, krpm],...
                                [scaledMdotAir, scaledMdotFuel]);

xlabel('Spool speed [kRPM]');
ylabel(hAx(1), 'Rate_{Air} [g/s],     Rate_{Fuel} [.01 g/s]') % left y-axis 
ylabel(hAx(2), 'Air:Fuel') % right y-axis

legendString = {'Air-fuel ratio', 'Air mass flow rate', ...
    'Fuel mass flow rate'};
legend(legendString, 'Location', 'southeast');
title('Air Mass Flow, Fuel Mass Flow, and Air-Fuel Ratio');
plotFixer();
print('-depsc','-tiff','-r300','plots/mdot,ratioVsRpm');

%% Deliverable 3
% Construct performance-metric plots showing how specific thrust, thrust
% specific fuel consumption, and thermal efficiency vary with spool speed.
% Use the lower heating value of Jet A (42,800 kJ/kg) for your thermal
% efficiency calculation, and base your performance metrics on the
% calculated thrust.
idealThrust = (mdotAir+mdotFuel).*V(:, end);

sp_thrust = idealThrust ./ mdotAir;
TSFC = idealThrust ./ mdotFuel; % Thrust Specific Fuel Consumption

powerIn = mdotFuel*const.LHVJetA; % [kW]
powerThrust = idealThrust.*V(:,end); 
thermal_eff = powerThrust ./ powerIn;

%% Specific thrust vs spool speed
plot(krpm, sp_thrust, '-o');
xlabel('Spool speed [kRPM]');
ylabel('Specific thrust [Ns/kg]');
title('Specific Thrust v. Spool Speed');
plotFixer();
print('-depsc','-tiff','-r300','plots/spthrustVsRpm');

%% Thrust specific fuel consumption
plot(krpm, 1e-3.*TSFC, '-o');
xlabel('Spool speed [kRPM]');
ylabel('Specific thrust [kN s/kg]');
title('Thrust Specific Fuel Consumption v. Spool Speed');
plotFixer();
print('-depsc','-tiff','-r300','plots/tsfcVsRpm');

%% Thermal Efficiency
% The thermal efficiency of the engine is defined as the ratio between the
% useful power of the engine over the power supplied.

plot(krpm, thermal_eff, '-o');
xlabel('Spool speed [kRPM]');
ylabel('Thermal Efficiency [1]');
title('Thermal Efficiency v. Spool Speed');
plotFixer();
print('-depsc','-tiff','-r300','plots/thermalEffVsRpm');

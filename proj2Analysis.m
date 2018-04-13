%% Project 2 Analysis
% This script is the main controller for the SR30 analysis in Project 2 of
% ME140. It should import data, and call the component functions for each
% of the individual parts of the project.
%
% Authors: Jean-Christophe Perrin, Beck Goodloe, Richard Randall, Jason
% Trinidad
%
% Created: 2018-04-11
% Edited: 2018-04-11
clear all;
clc;

%% Usefull Constants
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

%% Given Data
% The following values were supplied to us in the original project
% specifications and should not be changed.

given.pitotEffectiveArea = 6.4*const.insqToMsq; % [m^2]
given.jetAHeating = 42.8e6; % [J/kg/K]
given.A = [27.3, 6.4, 9.0, 7.2, 4.7, 3.87]'.*const.insqToMsq; % [m^2]

%% Import Collected Data
% Our collected data is in an excel file named `Data_proj2`. Here we import
% that data into a vector for use in our MatLab script. The excel file
% should never be modified by our script.

fname = 'data.txt';
%[NUM, TXT, RAW] = xlsread(fname, sname, range); % basic reading
stateData = readtable(fname);

%% Convert Data
% Before we start to do calculations on this data, let's convert it to the
% basic units that we will use throughout the course of our analysis. That
% means converting all temperature values into Kelvin and all pressure
% values into Pascals.

stateData.Properties.VariableUnits = {'', 'C', 'C', 'C', 'C', 'C', 'C', ...
    'kPa', 'kPa', 'kPa', 'kPa', 'kPa', 'kg/s', 'lbs'};
summary(stateData); % need to update values

%% Call a Helper function for each Deliverable
% Each function should take in the current table of calculated data and
% then add to it the component values that it calculates.

stateData = deliverable2(stateData, given, const);
stateData = deliverable3(stateData, given, const);
stateData = deliverable4(stateData, given, const);

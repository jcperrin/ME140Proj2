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

insqToMsq = 0.00064516; % [m^2/in^2]
lbfToN = 4.44822; % [N/lbf]

%% Given Data
% The following values were supplied to us in the original project
% specifications and should not be changed.

pitotEffectiveArea = 6.4*insqToMsq; % [m^2]
jetAHeating = 42.8e6; % [J/kg/K]
A = [27.3, 6.4, 9.0, 7.2, 4.7, 3.87]'.*insqToMsq; % [m^2]

%% Import Collected Data


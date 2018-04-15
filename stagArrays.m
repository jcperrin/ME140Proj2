function [T0, P0, Ma, U] = stagArrays(Data, const, given)
%stagArrays generates initial calculations

state(1).T = 

%% Stagnation Temperature v. RPM
% Let's collect all of the stagnation temperatures that we collect into a
% table. We can slice the variables of interest out of the table we
% imported earlier. However, we are still missing T1. This is the same for
% all trials -- the ambient air temperature of the lab. We generate a
% vector of this measurement (one for each trial) and collected all
% measurements in an array.

vars = {'T2', 'T3', 'T4', 'T5', 'T8'};
T0 = Data{:,vars} + const.KCdiff; % [K]
T1Vec = const.T1 * ones(nObservations, 1) + const.KCdiff; % [K]
T0 = [T1Vec, T0]; % [K] all T0 observations

end

function stateOut = nozzle(stateIn)
    disp('nozzle');
end

function stateOut = compressor(stateIn)
    disp('compressor');
end

function stateOut = combustor(stateIn)
    disp('combustor');
end

function stateOut = turbine(stateIn)
    disp('turbine');
end

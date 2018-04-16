function [outputTable, mdotAir] = solveEachLocation(Tm, Pm, mdot_fuel, real_thrust_lbs)
%% solveEachLocation
% This function solves for the values of interest at each location along
% the length of the jet engine that we are trying to analyze.
%
%% Syntax
%# table = solveEachLocation(Tm, Pm, mdot_fuel, real_thrust_lbs)
%
%% Description
% This function accepts vectors for the meausred values of temperature and
% pressure at each location of interest in the jet engine for one trial (a
% single value of RPM).
% 
% * Tm - A horizontal vector of the temperatures measured at locations 2-5, 8
% * Pm - A vector of the temperatures measured at locations 2-5, 8
% * mdot_fuel - A double for the rate of fuel consumption
% * real_thrust_lbs - Measured thrust produced by the engine
%
% * outputTable - all of the values of interest
% * mdot_air - mass of air flowing through engine
%
%% Example
% TODO
%
%% See Also
% TODO

% Testing solving for T0, P0, M from experimental values
%% Ensure all input values are formatted correctly
assert(isequal(size(Tm),  [1, 5]), "Tm must be a horizontal vector 5 wide");
assert(isequal(size(Pm),  [1, 5]), "Pm must be a horizontal vector 5 wide");
assert(isfloat(mdot_fuel), "mdot_fuel must be of type float");
assert(isfloat(real_thrust_lbs), "real_thrust_lbs must be of type float");

%% Constants
R = 287; % [kJ/kg/k] Gas constant of air
P_atm = 101325; % [Pa] atmospheric pressure
ZK = 273.15; % [deg] Celsius to metric conversion.
conv = 1/1550; % [m^2/in^2] square meters per square inch

% Collected Data
Tm_2 = Tm(1) + ZK; %21.3838+ZK;
Tm_3 = Tm(2) + ZK; %113.8878+ZK;
Tm_4 = Tm(3) + ZK; %581.9337+ZK;
Tm_5 = Tm(4) + ZK; %493.1745+ZK;
Tm_8 = Tm(5) + ZK; %478.7055+ZK;

Pd_2 = Pm(1)*1e3; %1000*1.5733;
P0_3 = Pm(2)*1e3 + P_atm; %1000*108.9267+P_atm;
P_4 = Pm(3)*1e3 + P_atm; %1000*103.8528+P_atm;
P0_5 = Pm(4)*1e3 + P_atm; %1000*12.4216+P_atm;
P0_8 = Pm(5)*1e3 + P_atm; %1000*9.7488+P_atm;

% Areas
A_1 = 27.3*conv;
A_2 = 6.4*conv;
A_3 = 9.0*conv;
A_4 = 7.2*conv;
A_5 = 4.7*conv;
A_8 = 3.87*conv;

% Recovery Factors
Rf_2 = 0.68;
Rf_3 = 0.68;
Rf_4 = 0.68;
Rf_5 = 0.86;
Rf_8 = 0.68;

%mdot_fuel = 0.0029;
%real_thrust_lbs = 6.7000;

% Code is decomposed so that all the iteration (though
% not all the math) is in the functions at the bottom.

%% Solve for state 2
P0_2 = P_atm;%Assumed same as atmospheric for intake
[P_2, M_2, T0_2, T_2, P0P_2] = solve2(P0_2, Pd_2, Tm_2, Rf_2, A_2, R);
%Do MFP & V calculations
MFP_2 = (M_2*sqrt(k(T_2))/P0P_2)*sqrt(1+k(T_2)*R*0.5*M_2*M_2/cp(T_2));
mdot_2 = A_2*MFP_2*P0_2/sqrt(R*T0_2);
V_2 = M_2*sqrt(k(T_2)*R*T_2);

%Store 2 in struct
states(2).M = M_2;
states(2).P0 = P0_2;
states(2).T0 = T0_2;
states(2).T = T_2;
states(2).V = V_2;
states(2).mdot = mdot_2;
mdotAir = mdot_2; % returned value

%% Solve state 1
[M_1,T_1] = solve1(mdot_2,P0_2,T0_2, A_1, R);%Ideal intake
V_1 = M_1*sqrt(k(T_1)*R*T_1);

%Store 1 in struct
states(1).M = M_1;
states(1).P0 = P0_2;
states(1).T0 = T0_2;
states(1).T = T_1;
states(1).V = V_1;

%% Solve state 3
[M_3, T_3, T0_3] = solve358(mdot_2, P0_3, Tm_3, Rf_3, A_3, R);
V_3 = M_3*sqrt(k(T_3)*R*T_3);
states(3).M = M_3;
states(3).P0 = P0_3;
states(3).T0 = T0_3;
states(3).T = T_3;
states(3).V = V_3;

%% Solve state 4
[M_4, T_4, T0_4, P0_4] = solve4(mdot_2+mdot_fuel, P_4, Tm_4, Rf_4, A_4, R);
V_4 = M_4*sqrt(k(T_4)*R*T_4);
states(4).M = M_4;
states(4).P0 = P0_4;
states(4).T0 = T0_4;
states(4).T = T_4;
states(4).V = V_4;

%% Solve state 5
[M_5, T_5, T0_5] = solve358(mdot_2+mdot_fuel, P0_5, Tm_5, Rf_5, A_5, R);
V_5 = M_5*sqrt(k(T_5)*R*T_5);
states(5).M = M_5;
states(5).P0 = P0_5;
states(5).T0 = T0_5;
states(5).T = T_5;
states(5).V = V_5;

%% Solve state 8
[M_8, T_8, T0_8] = solve358(mdot_2+mdot_fuel, P0_8, Tm_8, Rf_8, A_8, R);
V_8 = M_8*sqrt(k(T_8)*R*T_8);
states(8).M = M_8;
states(8).P0 = P0_8;
states(8).T0 = T0_8;
states(8).T = T_8;
states(8).V = V_8;

%Thrust as sanity check
thrust = (mdot_2+mdot_fuel)*V_8 - (mdot_2)*V_1;
thrust_lbs = thrust/4.45;
%disp(thrust_lbs);
%disp(real_thrust_lbs);

outputTable = struct2table(states);
end %overall function

%Begin functions for iteration

%Solve state 2 given Pd, P0, Tm, Rf, and A
function [P, M, T0, T, P0P] = solve2(P0, Pd, Tm, Rf, A, R)
    %Find M from P0/P; Find new T from M; repeat with new k(T).
    P = P0 - Pd;
    P0P = P0/P;
    precision = 0.1;
    M = 0;
    T0 = 0;
    %Guess M, iterate to get T, compute & check P0/P
    while(precision>0.00001)
        M = M + precision;
        T = 300;%guess
        n = 10;%In practice converges really fast
        while n>0
            n = n-1;
            k_ = k(T);
            cp_ = cp(T);
            T = Tm/(1+Rf*k_*R*0.5*M*M/cp_);
        end
        T0 = T*(1+k(T)*R*M*M*0.5/cp(T));
        integ = @(x) cp(x)./(R*x);
        P0P_est = exp(integral(integ, T, T0));%stag pressure ratio
        if(P0P_est>P0P)
            M = M - precision;
            precision = precision*0.1;
        end
    end
end


%Solve state 1 given mdot, P0, T0, A
function [M,T] = solve1(mdot, P0, T0, A, R)
    %Find M from MFP; find T from M; repeat with adjusted k(T), cp(T)
    T = T0; k_ = k(T); cp_ = cp(T);
    MFP = (mdot/A)*sqrt(R*T0)/P0;
    M = 0;
    n = 10;%Always enough to converge in practice
    while (n>0)
        T = T0/(1+k_*R*0.5*M*M/cp_);
        n = n - 1;
        k_ = k(T);
        cp_ = cp(T);
        precision = 0.1;
        while(precision>0.00001)
            M = M + precision;
            integ = @(x) cp(x)./(R*x);
            P0P_est = exp(integral(integ, T, T0));
            mfp_est = (M*sqrt(k_)/P0P_est)*sqrt(1+0.5*k_*R*M*M/cp_);
            if (MFP < mfp_est)
                M = M - precision;
                precision = precision/10;
            end
        end
    end
end


%Solve state 3, 5, or 8 given P0, Tm, A, mdot
function [M,T, T0] = solve358(mdot, P0, Tm, Rf, A, R)
    %Find M from MFP. Find T0,T,P0/P from M. Find MFP from P0,T0. Repeat.
    T = Tm; k_ = k(T); cp_ = cp(T);
    M = 0;
    n = 10;%enough to converge in practice
    while (n>0)
        n = n - 1;
        m = 10;
        while m>0%Iterate to find T given M
            m = m-1;
            k_ = k(T);
            cp_ = cp(T);
            T = Tm/(1+Rf*k_*R*0.5*M*M/cp_);
        end
        T0 = T*(1+k(T)*R*M*M*0.5/cp(T));%Find T0
        MFP = (mdot/A)*sqrt(R*T)/P0;
        
        precision = 0.1;
        while(precision>0.00001)
            M = M + precision;
            integ = @(x) cp(x)./(R*x);
            P0P_est = exp(integral(integ, T, T0));
            mfp_est = (M*sqrt(k_)/P0P_est)*sqrt(1+0.5*k_*R*M*M/cp_);
            if (MFP < mfp_est)
                M = M - precision;
                precision = precision/10;
            end
        end
    end
end

%Solve state 4 given P, Tm, A, mdot
function [M,T, T0, P0] = solve4(mdot, P, Tm, Rf, A, R)
    % Find M from MFP. Find T0,T,P0/P from M. Find P0 from P0/P. 
    % Find MFP from P0,T0. Repeat.
    T = Tm; k_ = k(T); cp_ = cp(T);
    P0 = P;%Estimate to begin iteration
    M = 0;
    n = 10;%enough to converge in practice
    while (n>0)
        n = n - 1;
        m = 10;
        while m>0%Iterate to find T given M
            m = m-1;
            k_ = k(T);
            cp_ = cp(T);
            T = Tm/(1+Rf*k_*R*0.5*M*M/cp_);
        end
        T0 = T*(1+k(T)*R*M*M*0.5/cp(T));%Find T0
        MFP = (mdot/A)*sqrt(R*T)/P0;
        P0P_est = 1;
        precision = 0.1;
        while(precision>0.00001)
            M = M + precision;
            integ = @(x) cp(x)./(R*x);
            P0P_est = exp(integral(integ, T, T0));
            mfp_est = (M*sqrt(k_)/P0P_est)*sqrt(1+0.5*k_*R*M*M/cp_);
            if (MFP < mfp_est)
                M = M - precision;
                precision = precision/10;
            end
        end
        P0 = P0P_est*P;
    end
end

%Helper methods for variable cp
function out = cp(T)
    airCoeffs = fliplr([28.11, 0.1967e-2, 0.4802e-5,-1.966e-9]);
    out = 1000*polyval(airCoeffs, T)/28.97;
end

function out = integral_cp(T)
    airCoeffs = fliplr([28.11, 0.1967e-2, 0.4802e-5,-1.966e-9]);
    intCoeffs = polyint(airCoeffs);
    1000*polyval(intCoeffs, T)/28.97;
end

function out = cp_avg(T1, T2)
    out = (integral_cp(T2) - integral_cp(T1)) / (T2-T1);
end

function out = integrate_cp(T1, T2)
    out = (integral_cp(T2) - integral_cp(T1));
end

function out = cv(T)
    R = 287;
    out = cp(T)-R;
end
function out = k(T)
    out = cp(T)/cv(T);
end
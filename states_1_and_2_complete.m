%Testing solving for T0, P0, M from experimental values

%Test values
ZK = 273.15;
Tm_2 = 21.3838+ZK;
Tm_3 = 113.8878+ZK;
Tm_4 = 581.9337+ZK;
Tm_5 = 493.1745+ZK;
Tm_8 = 478.7055+ZK;
Pd_2 = 1.5733;
P0_3 = 108.9267;
P_4 = 103.8528;
P0_5 = 12.4216;
P0_8 = 9.7488;
conv = 1/1550;%square meters per square inch
A_1 = 27.3*conv;
A_2 = 6.4*conv;
A_3 = 9.0*conv;
A_4 = 7.2*conv;
A_5 = 4.7*conv;
A_8 = 3.87*conv;
Rf_2 = 0.68;

%Constants
R = 287;

% Code is decomposed so that all the iteration (though
% not all the math) is in the functions at the bottom.

%Solve for state 2
P0_2 = 101.325;%Assumed same as atmospheric for intake
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

%Solve state 1
[M_1,T_1] = solve1(mdot_2,P0_2,T0_2, A_1, R);%Ideal intake
V_1 = M_1*sqrt(k(T_1)*R*T_1);

%Store 1 in struct
states(1).M = M_1;
states(1).P0 = P0_2;
states(1).T0 = T0_2;
states(1).T = T_1;
states(1).V = V_1;


%Solve state 2 given Pd, P0, Tm, Rf, and A
function [P, M, T0, T, P0P] = solve2(P0, Pd, Tm, Rf, A, R)
    P = P0 - Pd;
    P0P = P0/P;
    precision = 0.1;
    M = 0;
    T0 = 0;
    %Solve for M
    while(precision>0.00001)
        M = M + precision;
        T = 300;%guess
        n = 5;%In practice converges really fast
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
    T = T0;
    MFP = (mdot/A)*sqrt(R*T0)/P0;
    M = 0;
    n = 10;%Enough to converge in practice
    while (n>0)
        n = n - 1;
        k_ = k(T);
        cp_ = cp(T);
        precision = 0.1;
        while(precision>0.00001)
            M = M + precision;
            integ = @(x) cp(x)./(R*x);
            P0P_est = exp(integral(integ, T, T));
            mfp_est = (M*sqrt(k_)/P0P_est)*sqrt(1+0.5*k_*R*M*M/cp_);
            if (MFP < mfp_est)
                M = M - precision;
                precision = precision/10;
            end
        end
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
function [f,O] = dXdT_myo(t,X,Pdata,m,flg)
% Combined and corrected version of dXdT_myo, now strictly adhering to
% ode15s requirements for scalar t and column vector X.

% X is expected to be a column vector (N_states x 1) from ode15s

% If X is a column vector (from ode15s), transpose it to a row for consistency.
if iscolumn(X)
    X = X';
end

% INPUT - Interpolate input data at current time 't'
% These will be scalars as t is a scalar for ode15s calls
P_in = (interp1(Pdata(1,:),Pdata(2,:),t));
P_LV = (interp1(Pdata(1,:),Pdata(3,:),t));
dPdT = (interp1(Pdata(1,:),Pdata(4,:),t));
SIP  = m(1).*(interp1(Pdata(1,:),Pdata(5,:),t)); % m(1) scales P_im
dSIPdt  = m(1).*(interp1(Pdata(1,:),Pdata(6,:),t));
P_RA = (interp1(Pdata(1,:),Pdata(7,:),t));

% Intramyocardial Pressures (P_im) - these will be scalars
P_im1 = 0.167.*P_LV + SIP;
P_im2 = 0.500.*P_LV + SIP;
P_im3 = 0.833.*P_LV + SIP;
dPim1_dt = 0.167.*dPdT + dSIPdt;
dPim2_dt = 0.500.*dPdT + dSIPdt;
dPim3_dt = 0.833.*dPdT + dSIPdt;

% PARAMETERS
% Hardcoded parameters
C_PA = 0.0013*3;
L_PA = 0.2;
R_PA = 2.5;
R_PV = 1.5;
C_PV = 0.0254/3;
V01 = 2.5/3;
V02 = 8.0/3;
C2 = 0.254/9;
gamma = 0.75;

% Parameters from m-vector
cf1 = m(2); % epi/endo compliance scaling
rf1 = m(3); % epi/endo resistance scaling
R0m = m(4)*44;
R01 = m(5)*1.2*R0m;
R02 = 0.5*R0m;
C1 = m(6)*0.013/9;

% STATE VARIABLES (X is a column vector)
P_PA = X(:,1); Q_PA = X(:,2); P11 = X(:,3); P21 = X(:,4);
P12 = X(:,5); P22 = X(:,6); P13 = X(:,7); P23 = X(:,8); P_PV = X(:,9);

% CALCULATIONS (These will be scalars)
% Epicardial Layer (using 1/cf1 and 1/rf1 based on dXdT_myo logic)
V11 = (1./cf1).*((P11 - P_im1).*C1+V01);
V21 = (1./cf1).*((P21 - P_im1).*C2+V02);
R11 = (1./rf1).*R01.*(V01./V11).^2;
R21 = (1./rf1).*R02.*(V02./V21).^2;
Rm1 = R0m.*(gamma.*R11./R01 + (1-gamma).*R21./R02);
Q11 = (P_PA - P11)./R11;
Qm1 = (P11 - P21)./Rm1;
Q21 = (P21 - P_PV)./R21;

% Mid-wall Layer (reference layer)
V12 = ((P12 - P_im2).*C1+V01);
V22 = ((P22 - P_im2).*C2+V02);
R12 = R01.*(V01./V12).^2;
R22 = R02.*(V02./V22).^2;
Rm2 = R0m.*(gamma.*R12./R01 + (1-gamma).*R22./R02);
Q12 = (P_PA - P12)./R12;
Qm2 = (P12 - P22)./Rm2;
Q22 = (P22 - P_PV)./R22;

% Endocardial Layer (using cf1 and rf1 based on dXdT_myo logic)
V13 = cf1.*((P13 - P_im3).*C1+V01);
V23 = cf1.*((P23 - P_im3).*C2+V02);
R13 = rf1.*R01.*(V01./V13).^2;
R23 = rf1.*R02.*(V02./V23).^2;
Rm3 = R0m.*(gamma.*R13./R01 + (1-gamma).*R23./R02);
Q13 = (P_PA - P13)./R13;
Qm3 = (P13 - P23)./Rm3;
Q23 = (P23 - P_PV)./R23;

% Total Flows
Q_ima = Q11 + Q12 + Q13;
Q_imv = Q21 + Q22 + Q23;
Q_out = (P_PV - P_RA)./R_PV;

% Volumes
V_PA = P_PA.*C_PA;
V_PV = P_PV.*C_PV;

% DIFFERENTIAL EQUATIONS (f = dX/dt)
% Initialize f as a column vector
f = zeros(size(X));
f(:,1) = (Q_PA - Q_ima)./C_PA; % P_PA
f(:,2) = (P_in - P_PA - Q_PA.*R_PA)./L_PA; % Q_PA
f(:,3) = (Q11-Qm1)./(C1./cf1) + dPim1_dt; % P11
f(:,4) = (Qm1-Q21)./(C2./cf1) + dPim1_dt; % P21
f(:,5) = (Q12-Qm2)./(C1) + dPim2_dt; % P12
f(:,6) = (Qm2-Q22)./(C2) + dPim2_dt; % P22
f(:,7) = (Q13-Qm3)./(cf1.*C1) + dPim3_dt; % P13
f(:,8) = (Qm3-Q23)./(cf1.*C2) + dPim3_dt; % P23
f(:,9) = (Q_imv - Q_out)./C_PV; % P_PV

% Transpose f back to a column vector if the input was a column (for ode15s).
if iscolumn(t)
    f = f';
end

% Conditional Output 'O'
if flg == 0
    O = [];
else
    % Return intermediate values as a row vector for a single time point
    O= [Q11, Q12, Q13, Qm1, Qm2, Qm3, Q21, Q22, Q23, ...
        V11, V12, V13, V21, V22, V23, P_im1, P_im2, P_im3, ...
        R11, R12, R13, R21, R22, R23, Rm1, Rm2, Rm3, ...
        V_PA, V_PV, Q_out];
end
end
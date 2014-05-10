function [A,Q] = characteristicOutflow()

global nnodes                   % number of spatial nodes
global nsteps                   % number of time-steps
global Tcard numPeriods         % the cardiac period being simulated
global CardiacOut
global CardiacPeak
global X0 X1                    % starting and ending coordinate
global Xvec
global omega0 omega1            % range of frequencies
global delK delH                % time-step size and mesh grid spacing
global rootR isR0function       % properties of the geometry
global Elast wallH wallDelta    % mechanical properties of vessel wall
global P0
global Nu rho Mu                % physical properties of the fluid
global Rtop Rbottom
global errTol 
global solA solQ                % arrays for storing solution

%----------------------------------------------------
% Each of these arrays are of size numtimesteps x 2. 
% Column 1 stores solutions for previous time-period
% Column 2 stores solutions for current time-period
%----------------------------------------------------
global solQ_mp_np_half          % array to store Q_{m+1/2}^{n+1/2}
global solA_mp_np_half          % array to store A_{m+1/2}^{n+1/2}
global solQ_m_n                 % array to store Q_m^n
global solA_m_n                 % array to store A_m^n
global solQ_mm1_n               % array to store Q_{m-1}^n
global solA_mm1_n               % array to store A_{m-1}^n

global isMatlabNonLinSolve 
global isMaterialExpModel
global expModelK1 expModelK2 expModelK3

% Calculate R0 and dR0dX:
%------------------------
LenVessel = X1-X0;
if ( isR0function )
    R0      = getR0function(Xvec, Rtop, Rbottom, LenVessel);
    dR0dX   = diffR0function(Xvec, Rtop, Rbottom, LenVessel);
else
    R0 = getR0data(Xvec);       % NEED TO IMPLEMENT THIS FUNCTION
    dR0dX = 0;
end

% Calculate A0:
%--------------
A0 = pi*R0*R0;

% We are going from a_N-1 to a_N in this entire function:
%--------------------------------------------------------

% Calculate A_A = A(time-1,endnode-1):
%-------------------------------------
A_A = solA(a_N-1, end-1);

% Calculate Q_A = Q(time-1, endnode-1):
%--------------------------------------
Q_A = solQ(a_N-1, end-1);

% Calculate C_A = wavSpeed(time-1, end-1):
%-----------------------------------------
C_A = getWavespeed(Xvec(end-1), A_A);

% Calculate A_T = A(time-1, end):
%--------------------------------
A_T = solA(a_N-1, end);

% Calculate Q_T = Q(time-1, end):
%--------------------------------
Q_T = solQ(a_N-1, end);

% Calculate C_T = waveSpeed(time-1,endNode):
%-------------------------------------------
C_T = getWavespeed(Xvec(end),A_T);

if ( isCopyOlufsen == 1 )
else
    x0 = [A_T, Q_T, C_T];
    g1 = @(x) A_T - ((x(2)/x(1)) + x(3))*(delK/delH)*(A_T - A_A);
    g2 = @(x) Q_T - ((x(2)/x(1)) + x(3))*(delK/delH)*(Q_T - Q_A);
    g3 = @(x) C_T - ((x(2)/x(1)) + x(3))*(delK/delH)*(C_T - C_A);
end
% Calculate Q_R

% Calculate R_end
%----------------
R_end = sqrt(A_R/pi);

% Calculate speedC0

if ( isMaterialExpModel == 1 ) 
    f       = (4.0/3.0)*(expModelK1*exp(expModelK2*R0)+expModelK3);
    dfdR0   = (4.0/3.0)*(expModelK1*expModelK2*exp(expModelK2*R0));
else
    f       = (4.0*Elast*wallH)/(3.0*R0);
    dfdR0   = -(4.0*Elast*wallH)/(3.0*R0*R0);
end

coeff_C1 = dfdR0*dR0dX;

coeff_C2 = -(dfdR0*sqrt(A0) + f*sqrt(pi))*dR0dX;

coeff_C3 = 2.0*pi*Nu*Q_R*R_end/wallDelta;

H_R_plus = ((coeff_C1*A_R + coeff_C2*sqrt(A_R))/rho) + (coeff_C3/A_R);
H_R_plus = H_R_plus/(-((Q_R/A_R) + speedC0));

% Calculate Q_tms

% Calculate Q_q

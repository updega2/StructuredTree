function f3 = getNonLinearEq3(a_SolX, a_N, a_Admittances)
%--------------------------------------------------------------------------
% Function to calculate the non-linear equation corresponding to the
% outflow boundary condition for variable Q_{m}^{n+1} as defined by
% Olufsen's thesis
% Argument definitions:
% a_SolX:           The current set of guessed values for the 
%                   iteration variables
% a_N:              The index of the current time step in 
%                   the cardiac period
% a_Admittances:    The inverse fourier transformed outflow 
%                   admittance values
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------

%Define global variables
global Xvec
global delK delH                % time-step size and mesh grid spacing
global isR0function             % property of the geometry
global solA solQ                % arrays for storing solution
global P0
global Rtop Rbottom 
global X1 X0

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

% Define the length of the vessel from X1 and X0:
%------------------------------------------------
LenVessel = X1 - X0;

% Define theta value from discretization stencil:
%------------------------------------------------
%theta = delK/(2.0*delH);   %%%% OLD CODE
theta  = delK/delH;         %%%% NEW CODE

% A_m_n   = solA(a_N, end);     %%%% OLD CODE
% A_mm1_n = solA(a_N, end-1);   %%%% OLD CODE
% Q_m_n   = solQ(a_N, end);     %%%% OLD CODE
% Q_mm1_n = solQ(a_N, end-1);   %%%% OLD CODE

A_m_n   = solA_m_n(a_N,2);      %%%% NEW CODE
A_mm1_n = solA_mm1_n(a_N,2);    %%%% NEW CODE
Q_m_n   = solQ_m_n(a_N,2);      %%%% NEW CODE
Q_mm1_n = solQ_mm1_n(a_N,2);    %%%% NEW CODE

% calculate R0(X_m):
%-------------------
X_m = Xvec(end);
if ( isR0function )
    R0_m    = getR0function(X_m,Rtop, Rbottom, LenVessel);
    dR0dX_m = diffR0function(X_m,Rtop, Rbottom, LenVessel);
else
    R0_m    = getR0data(X_m);    %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_m = 0.0;               %%%% HAVE NOT ENCODED THIS
end

% calculate R(2)_{m}^{n}:
%--------------------------
R_m_n     = getFluxVector(A_m_n,Q_m_n,R0_m,P0);
R2_m_n    = R_m_n(2);

% calculate S(2)_{m}^{n}:
%--------------------------
S_m_n     = getSourceVector(A_m_n,Q_m_n,R0_m,dR0dX_m,P0);
S2_m_n    = S_m_n(2);

% calculate R0(X_{m-1}):
%-----------------------
X_mm1 = Xvec(end-1);
if ( isR0function )
    R0_mm1    = getR0function(X_mm1,Rtop, Rbottom, LenVessel);
    dR0dX_mm1 = diffR0function(X_mm1,Rtop, Rbottom, LenVessel);
else
    R0_mm1    = getR0data(X_mm1);     %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_mm1 = 0.0;                  %%%% HAVE NOT ENCODED THIS
end

% calculate R(2)_{m-1}^{n}:
%----------------------------
R_mm1_n     = getFluxVector(A_mm1_n,Q_mm1_n,R0_mm1,P0);
R2_mm1_n    = R_mm1_n(2);

% calculate S(2)_{m-1}^{n}:
%----------------------------
S_mm1_n     = getSourceVector(A_mm1_n,Q_mm1_n,R0_mm1,dR0dX_mm1,P0);
S2_mm1_n    = S_mm1_n(2);

% calculate Q_{m-1/2}^{n+1/2} using the discretization scheme:
% Q_{m-1/2}^{n+1/2} = \frac{Q_m_n + Q_{m-1}^n}{2} 
%           - \frac{k}{2}\frac{R(2)_m^n - R(2)_{m-1}^n}{h}
%           + \frac{k}{2}\frac{S(2)_m^n + S(2)_{m-1}^n}{2}
%-------------------------------------------------------------
Q_mm_np_half = 0.5*(Q_m_n + Q_mm1_n) - ...
    (delK/(2*delH))*(R2_m_n - R2_mm1_n) + ...
    (delK/4)*(S2_m_n + S2_mm1_n);

% calculate k_5 = A_m^n + theta*Q_{m-1/2}^{n+1/2}:
%-------------------------------------------------
k5 = A_m_n + theta*Q_mm_np_half

a_SolX(1)
a_SolX(4)

f3 = k5 - theta*a_SolX(1) - a_SolX(4);
function Df = getJacobianX(a_SolX,a_N,a_Admittances)
%--------------------------------------------------------------------------
% Function to calculate the Jacobian term for the non-linear equation 
% corresponding to the outflow boundary condition for variableS
% Q_{m+1/2}^{n+1/2}, A_{m+1/2}{n+1/2}, Q_m^{n+1}, and A_m^{n+1}
% as defined by Olufsen's thesis
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
global delK delH nsteps               % time-step size and mesh grid spacing
global isR0function             % properties of the geometry
global Elast wallH wallDelta    % mechanical properties of vessel wall
global Nu                       % physical property of the fluid
global P0
global solA solQ

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

global Rtop Rbottom
global X1 X0

LenVessel = X1 - X0;

Df = ones(4,4)*10^-30;

% Define theta and gamma from the discretization stencil:
%--------------------------------------------------------
%theta = delK/(2*delH); %%%% OLD CODE
%gamma = delK/4;        %%%% OLD CODE

theta = delK/delH;      %%%% NEW CODE
gamma = delK/2.0;       %%%% NEW CODE

% get X_m, and R0(X_m), and (dR0/dX) evaluated at X_m
%----------------------------------------------------
X_m = Xvec(end);
if ( isR0function )
    R0_m    = getR0function(X_m,Rtop, Rbottom, LenVessel);
    dR0dX_m = diffR0function(X_m,Rtop, Rbottom, LenVessel);
else
    R0_m    = getR0data(X_m);       %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_m = 0.0;                  %%%% HAVE NOT ENCODED THIS
end
A0_m    = pi*R0_m*R0_m;

% X_{m+1/2} = 0.5*(X_m + X_{m-1}), and get R0(X_{m+1/2}), and dR0/dX
% evaluated at X_{m+1/2}: (IS THIS FORMULA CORRECT ???)
%-------------------------------------------------------------------
%X_mp_half   = 0.5*(Xvec(end) + Xvec(end-1));
X_mp_half = Xvec(end) + 0.5*delH;            %%%% CHECK THIS IMPLEMENTATION
if ( isR0function )
    R0_mp_half    = getR0function(X_mp_half,Rtop, Rbottom, LenVessel);
    dR0dX_mp_half = diffR0function(X_mp_half,Rtop, Rbottom, LenVessel);
else
    R0_mp_half    = getR0data(X_mp_half); %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_mp_half = 0.0;                  %%%% HAVE NOT ENCODED THIS
end
A0_mp_half  = pi*R0_mp_half*R0_mp_half;

% A_m_n   = solA(a_N,end);    %%%% OLD CODE
% Q_m_n   = solQ(a_N, end);   %%%% OLD CODE
% A_mm1_n = solA(a_N,end-1);  %%%% OLD CODE
% Q_mm1_n = solQ(a_N,end-1);  %%%% OLD CODE
A_m_n   = solA_m_n(a_N,2);      %%%% NEW CODE
Q_m_n   = solQ_m_n(a_N,2);      %%%% NEW CODE
A_mm1_n = solA_mm1_n(a_N,2);    %%%% NEW CODE
Q_mm1_n = solQ_mm1_n(a_N,2);    %%%% NEW CODE

X_mm1   = Xvec(end-1);
X_m     = Xvec(end);
if ( isR0function )
    R0_m    = getR0function(X_m, Rtop, Rbottom, LenVessel);
    dR0dX_m = diffR0function(X_m, Rtop, Rbottom, LenVessel);
else
    R0_m    = getR0data(X_m);          %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_m = 0;                                 %%%% HAVE NOT ENCODED THIS
end
if ( isR0function )
    R0_mm1      = getR0function(X_mm1, Rtop, Rbottom, LenVessel);
    dR0dX_mm1   = diffR0function(X_mm1, Rtop, Rbottom, LenVessel);
else
    R0_mm1      = getR0data(X_mm1);    %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_mm1   = 0;                   %%%% HAVE NOT ENCODED THIS
end

R_m_n       = getFluxVector(A_m_n, Q_m_n, R0_m, P0);
S_m_n       = getSourceVector(A_m_n, Q_m_n, R0_m, dR0dX_m, P0);
R1_m_n      = R_m_n(1);
S1_m_n      = S_m_n(1);
R_mm1_n     = getFluxVector(A_mm1_n, Q_mm1_n, R0_mm1, P0);
S_mm1_n     = getSourceVector(A_mm1_n, Q_mm1_n, R0_mm1, dR0dX_mm1, P0);
R1_mm1_n    = R_mm1_n(1);
S1_mm1_n    = S_mm1_n(1);

A_mm_np_half = 0.5*(A_m_n + A_mm1_n) - ...
    (delK/(2.0*delH))*(R1_m_n - R1_mm1_n) + ...
    (delK/4.0)*(S1_m_n + S1_mm1_n);

k3 = A_mm_np_half;

% calculate (using analytical derivatives based on the assumed form of
% pressure-area law) (dP/dA) at A = 0.5*(k3+x_2), at X_m:
%----------------------------------------------------------------------
DpDx2_m = 0.5*(2*Elast*wallH/(3*R0_m))*...
    (sqrt(A0_m)/sqrt((0.5*(k3+a_SolX(2)))^3));

% calculate (using analytical derivatives based on the assumed form of
% pressure-area law) (dP/dA) at A = x_4, at X_m:
%----------------------------------------------------------------------
DpDx4_m = (4*Elast*wallH/(3*R0_m))*(sqrt(A0_m)/sqrt((a_SolX(4))^3));

% calculate (using analytical derivatives based on the assumed form of
% pressure-area law and the derived form of B) (dB/dx_2) at X_{m+1/2}:
%----------------------------------------------------------------------
DBDx2_mp_half = (4*Elast*wallH/(3*R0_mp_half))*...
    (0.5*A0_mp_half)/sqrt(a_SolX(2)*A0_mp_half);

% calculate R(X_{m+1/2}):
%------------------------
R_mp_half = sqrt(a_SolX(2)/pi);

% calculate (using the assumed form of the wall shear stress) (dF/dQ) at
% Q=x_1, at X_{m+1/2}:
%-----------------------------------------------------------------------
DFDx1_mp_half = -2.0*pi*Nu*R_mp_half/(wallDelta*a_SolX(2));

% calculate (using the assumed form of the wall shear stress) (dF/dA) at
% A+x_2, at X_{m+1/2}:
%-----------------------------------------------------------------------
DFDx2_mp_half = 2.0*pi*Nu*R_mp_half*a_SolX(1)/(wallDelta*a_SolX(2)*a_SolX(2));

% evaluate f = 4Eh/3R0 at X_{m+1/2}:
%-----------------------------------
f_mp_half = 4.0*Elast*wallH/(3.0*R0_mp_half);

% evaluate df/dR0 for f defined above, at X_{m+1/2}:
%---------------------------------------------------
DfDR0_mp_half = -4.0*Elast*wallH/(3.0*R0_mp_half*R0_mp_half);

% evaluate d(dB/dx)/dx_2 at X_{m+1/2}:
%-------------------------------------
D2BDxDx2_mp_half = ((1.0/sqrt(a_SolX(2)))*(sqrt(pi)*f_mp_half ...
    + sqrt(A0_mp_half)*DfDR0_mp_half) - DfDR0_mp_half)*dR0dX_mp_half;

% now for evaluation of the fourier term (Q_{tms})_m^{n+1}:
%----------------------------------------------------------
Qtms_m_np1 = 0;

for countK = 2:nsteps
    
    if ( a_N - countK + 0.5 >= 0 )
        % calculate Q_{m+1/2}^{n-k+1/2} from current period:
        %---------------------------------------------------
        %Q_mp_nmkp_half  = solQ_mphalf(a_N-countK,2); %%%% CHECK INDEXING 
        Q_mp_nmkp_half  = solQ_mp_np_half(countK,2);      %%%% CHECK INDEXING
        
        % calculate A_m^{n-k} from current period:
        %-----------------------------------------
        %A_m_nmk         = solA_m(a_N-countK,2); %%%% CHECK INDEXING
        A_m_nmk         = solA_m_n(countK,2);      %%%% CHECK INDEXING
    else
        % calculate Q_{m+1/2}^{n-k+1/2} from previous period:
        %----------------------------------------------------
        %Q_mp_nmkp_half  = solQ_mphalf(a_N-countK,1);  %%%% CHECK INDEXING
        Q_mp_nmkp_half  = solQ_mp_np_half(countK,1);       %%%% CHECK INDEXING
        
        % calculate A_m^{n-k} from previous period:
        %------------------------------------------
        %A_m_nmk         = solA_m(a_N-countK,1); %%%% CHECK INDEXING
        A_m_nmk         = solA_m_n(countK,1);      %%%% CHECK INDEXING
    end
    
    if ( a_N - countK >= 0 )
        % calculate Q_m^{n-k} from current period:
        %-----------------------------------------
        %Q_m_nmk = solQ_m(a_N-countK,2);            %%%% CHECK INDEXING
        Q_m_nmk = solQ_m_n(countK,2);               %%%% CHECK INDEXING
        
        % calculate Q_{m-1}^{n-k} from current period:
        %---------------------------------------------
        %Q_mm1_nmk = solQ_mm1(a_N-countK,2);        %%%% CHECK INDEXING
        Q_mm1_nmk = solQ_mm1_n(countK,2);           %%%% CHECK INDEXING
        
        % calculate A_m^{n-k} from current period:
        %-----------------------------------------
        %A_m_nmk = solA_m(a_N-countK,2);            %%%% CHECK INDEXING
        A_m_nmk = solA_m_n(countK,2);               %%%% CHECK INDEXING
        
        % calculate A_{m-1}^{n-k} from current period:
        %---------------------------------------------
        %A_mm1_nmk = solA_mm1(a_N-countK,2);        %%%% CHECK INDEXING
        A_mm1_nmk = solA_mm1_n(countK,2);           %%%% CHECK INDEXING
    else
        % calculate Q_m^{n-k} from previous period:
        %------------------------------------------
        %Q_m_nmk = solQ_m(a_N-countK,1);            %%%% CHECK INDEXING
        Q_m_nmk = solQ_m_n(countK,1);               %%%% CHECK INDEXING
        
        % calculate Q_{m-1}^{n-k} from previous period:
        %----------------------------------------------
        %Q_mm1_nmk = solQ_mm1(a_N-countK,1);        %%%% CHECK INDEXING
        Q_mm1_nmk = solQ_mm1_n(countK,1);           %%%% CHECK INDEXING
        
        % calculate A_m^{n-k} from previous period:
        %------------------------------------------
        %A_m_nmk = solA_m(a_N-countK,1);            %%%% CHECK INDEXING
        A_m_nmk = solA_m_n(countK,1);               %%%% CHECK INDEXING
        
        % calculate A_{m-1}^{n-k} from previous period:
        %----------------------------------------------
        %A_mm1_nmk = solA_mm1(a_N-countK,1);        %%%% CHECK INDEXING
        A_mm1_nmk = solA_mm1_n(countK,1);           %%%% CHECK INDEXING
    end
    
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
    
    % calculate R(2)_{m}^{n-k}:
    %--------------------------
    R_m_nmk     = getFluxVector(A_m_nmk,Q_m_nmk,R0_m,P0);
    R2_m_nmk    = R_m_nmk(2);
    
    % calculate S(2)_{m}^{n-k}:
    %--------------------------
    S_m_nmk     = getSourceVector(A_m_nmk,Q_m_nmk,R0_m,dR0dX_m,P0);
    S2_m_nmk    = S_m_nmk(2);
    
    % calculate R0(X_{m-1}):
    %-----------------------
    X_mm1 = Xvec(end-1);
    if ( isR0function )
        R0_mm1    = getR0function(X_mm1,Rtop, Rbottom, LenVessel);
        dR0dX_mm1 = diffR0function(X_mm1,Rtop, Rbottom, LenVessel);
    else
        R0_mm1    = getR0data(X_mm1); %%%% NEED TO IMPLEMENT THIS FUNCTION
        dR0dX_mm1 = 0.0;              %%%% HAVE NOT ENCODED THIS
    end
    
    % calculate R(2)_{m-1}^{n-k}:
    %----------------------------
    R_mm1_nmk     = getFluxVector(A_mm1_nmk,Q_mm1_nmk,R0_mm1,P0);
    R2_mm1_nmk    = R_mm1_nmk(2);
    
    % calculate S(2)_{m-1}^{n-k}:
    %----------------------------
    S_mm1_nmk     = getSourceVector(A_mm1_nmk,Q_mm1_nmk,...
        R0_mm1,dR0dX_mm1,P0);
    S2_mm1_nmk    = S_mm1_nmk(2);
    
    % use these calculated values to plug into the discretisation for
    % Q_{m-1/2}^{n-k+1/2}:
    %----------------------------------------------------------------
    Q_mm_nmkp_half = 0.5*(Q_m_nmk + Q_mm1_nmk) - ...
        (delK/(2.0*delH))*(R2_m_nmk - R2_mm1_nmk) + ...
        (delK/4.0)*(S2_m_nmk + S2_mm1_nmk);
    
    % use the value of Q_{m-1/2}^{n-k+1/2} to get A_m^{n-k+1}:
    %---------------------------------------------------------
    A_m_nmkp1 = A_m_nmk - (delK/4.0)*(Q_mp_nmkp_half - Q_mm_nmkp_half);
    
    % use this to get P(X_m, A_m^{n-k+1}):
    %-------------------------------------
    P_nmkp1 = getPressure(A_m_nmkp1,R0_m,P0);
    
    % use this to get (Q_{tms})_{m}^{n+1} = sum P(X_m,
    % A_m^{n-k+1})y_m^k\Delta t:
    %-------------------------------------------------
    Qtms_m_np1 = Qtms_m_np1 + P_nmkp1*a_Admittances(countK)*delK;
end

k2 = a_Admittances(1)*delK;
k4 = Qtms_m_np1;

% use the definitions of z1,z2,z3,z4 as given in Olufsen's thesis:
%-----------------------------------------------------------------
z1 = k2*DpDx2_m;
%z2 = k4*DpDx4_m; %%%% OLD CODE
z2 = k2*DpDx4_m;  %%%% NEW CODE
z3 = -theta*2*a_SolX(1)/a_SolX(2) + gamma*DFDx1_mp_half;
%z4 = theta*(((a_SolX(1)^2)/(a_SolX(2)^2)) + DBDx2_mp_half) + ...
%    gamma*(DFDx2_mp_half + D2BDxDx2_mp_half);         %%%% OLD CODE
z4 = theta*(((a_SolX(1)^2)/(a_SolX(2)^2)) - DBDx2_mp_half) + ...
    gamma*(DFDx2_mp_half + D2BDxDx2_mp_half);          %%%% NEW CODE

Df(1,:) = [-0.5, z1, 0.0, 0.0];
Df(2,:) = [0.0, 0.0, -1.0, z2];
Df(3,:) = [-theta, 0.0, 0.0, -1.0];
Df(4,:) = [z3, z4, -1.0, 0.0];
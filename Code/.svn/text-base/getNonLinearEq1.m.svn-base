function f1 = getNonLinearEq1(a_SolX, a_N, a_Admittances)
%--------------------------------------------------------------------------
% Function to calculate the non-linear equation corresponding to the
% outflow boundary condition for variable Q_{m+1/2}^{n+1/2} as defined by
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
global solQ solA
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

global Rtop Rbottom X1 X0
global nsteps
global P0

% Define length of the vessel from X1 and X0:
%--------------------------------------------
LenVessel = X1 - X0;

% Access Q_{m+1/2}^{n+1/2}:
%--------------------------
Q_mp_np_half    = a_SolX(1);

% Access A_{m+1/2}^{n+1/2}:
%--------------------------
A_mp_np_half    = a_SolX(2);

% Older Version:
%if a_N ~= 1
%    A_m_n           = solA_m_n(a_N,2);
%    A_m_nm1         = solA_m_n(a_N-1,2);
%else
%    A_m_n           = solA_m_n(a_N,2);
%    A_m_nm1         = solA_m_n(nsteps,2);
%end

%A_mm_half       = (A_m_n + A_m_nm1)/2.0;

% Begin Newer Version:
A_m_n   = solA_m_n(a_N, 2);     % value from current period
Q_m_n   = solQ_m_n(a_N, 2);     % value from current period
A_mm1_n = solA_mm1_n(a_N, 2);   % value from current period
Q_mm1_n = solQ_mm1_n(a_N, 2);   % value from current period

% get nodal coordinate X_m, and R0_m = R0(X_m):
%----------------------------------------------
Xm = Xvec(end);
if ( isR0function )
    R0_m    = getR0function(Xm, Rtop, Rbottom, LenVessel);
    dR0dX_m = diffR0function(Xm, Rtop, Rbottom, LenVessel);
else
    R0_m    = getR0data(Xm); %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_m = 0;
end

Xmm1 = Xvec(end-1);
if ( isR0function )
    R0_mm1 = getR0function(Xmm1, Rtop, Rbottom, LenVessel);
    dR0dX_mm1 = diffR0function(Xmm1, Rtop, Rbottom, LenVessel);
else
    R0_mm1 = getR0data(Xmm1); %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_mm1 = 0;
end
% calculate R_{j+1/2}^T for j = M-1/2, T = N to get R_M^N:
%---------------------------------------------------------
R_m_n   = getFluxVector(A_m_n, Q_m_n, R0_m, P0);

% calculate S_{j+1/2}^T for j = M-1/2, T = N to get S_M^N:
%---------------------------------------------------------
S_m_n   = getSourceVector(A_m_n, Q_m_n, R0_m, dR0dX_m, P0);

% extract the individual vector components out:
%----------------------------------------------
R1_m_n  = R_m_n(1);
R2_m_n  = R_m_n(2); % can comment this out
S1_m_n  = S_m_n(1);
S2_m_n  = S_m_n(2); % can comment this out

% calculate R_{j-1/2}^T for j = M-1/2, T = N to get R_{M-1}^N:
%-------------------------------------------------------------
R_mm1_n = getFluxVector(A_mm1_n, Q_mm1_n, R0_mm1, P0);

% calculate S_{j-1/2}^T for j = M-1/2, T = N to get S_{M-1}^N:
%-------------------------------------------------------------
S_mm1_n = getSourceVector(A_mm1_n, Q_mm1_n, R0_mm1, dR0dX_mm1, P0);

% extract the individual vector components out:
%----------------------------------------------
R1_mm1_n = R_mm1_n(1);
R2_mm1_n = R_mm1_n(2); % can comment this out
S1_mm1_n = S_mm1_n(1);
S2_mm1_n = S_mm1_n(2); % can comment this out

% calculate A_{M-1/2}^{T+1/2} using the difference scheme
% A_{M-1/2}^{N+1/2} = \frac{A_M^N + A_{m-1}^N}{2} 
%           -\frac{k}{2}\frac{(R(1)_M^N - R(1)_{M-1}^N)}{h}
%           \frac{k}{2}\frac{S(1)_M^N + S(1)_{M-1}^N}{2}
%-----------------------------------------------------------
A_mm_np_half = (A_m_n + A_mm1_n) ...
    -(delK/(2.0*delH))*(R1_m_n - R1_mm1_n)...
    +(delK/4.0)*(S1_m_n + S1_mm1_n);
% End Newer Version

% Define k_3 = A_{M - 1/2}^{N + 1/2}
%-----------------------------------
k3 = A_mm_np_half;

% Define k_2 = y(X_m,T_n)\times \Delta t, which is the same as y_m^k at k=0
%--------------------------------------------------------------------------
k2 = a_Admittances(1)*delK;

% get nodal coordinate X_m, and R0_m = R0(X_m):
%----------------------------------------------
Xm = Xvec(end);
if ( isR0function )
    R0m = getR0function(Xm, Rtop, Rbottom, LenVessel);
else
    R0m = getR0data(Xm); %%%% NEED TO IMPLEMENT THIS FUNCTION
end

% calculate \frac{k_3 + x_2}{2}
%------------------------------
k3_p_x2_half    = 0.5*(k3 + a_SolX(2));

% calculate pressure P(X_m, \frac{k_3+x_2}{2}):
%----------------------------------------------
Pm_k3_p_x2_half = getPressure(R0m, k3_p_x2_half, P0);

% NOTE BY Debanjan: Revised till here on 05/06/2014 6:50 PM

Qtms_m_np_half  = 0.0;

for countK = 2:nsteps
    if (a_N - countK + 0.5 >= 0.0)
        % calculate A_{m+1/2}^{n-k+1/2} from current period:
        %---------------------------------------------------
        %A_mp_nmk_half = solA_mphalf(a_N-countK,2);     %%%% CHECK INDEXING
        A_mp_nmkp_half = solA_mp_np_half(countK,2);     %%%% CHECK INDEXING
        
        % calculate Q_{m+1/2}^{n-k+1/2} from current period:
        %---------------------------------------------------
        %Q_mp_nmk_half = solQ_mphalf(a_N-countK,2);     %%%% CHECK INDEXING
        Q_mp_nmkp_half = solQ_mp_np_half(countK,2);     %%%% CHECK INDEXING
    else
        % calculate A_{m+1/2}^{n-k+1/2} from previous period:
        %----------------------------------------------------
        %A_mp_nmk_half = solA_mphalf(a_N-countK,1);     %%%% CHECK INDEXING
        A_mp_nmkp_half = solA_mp_np_half(countK,1);     %%%% CHECK INDEXING
        
        % calculate A_{m+1/2}^{n-k+1/2} from previous period:
        %----------------------------------------------------
        %Q_mp_nmk_half = solQ_mphalf(a_N-countK,1);     %%%% CHECK INDEXING
        Q_mp_nmkp_half = solQ_mp_np_half(countK,1);     %%%% CHECK INDEXING
    end
    
    if (a_N - countK >= 0.0) %%% OLDER VERSION
        % calculate A_m^{n-k} from current period:
        %-----------------------------------------
        %A_m_nmk = solA_m(a_N-countK,2);                %%%% CHECK INDEXING
        A_m_nmk = solA_m_n(countK,2);                   %%%% CHECK INDEXING
        
        % calculate A_{m-1}^{n-k} from current period:
        %---------------------------------------------
        %A_mm1_nmk = solA_mm1(a_N-countK,2);            %%%% CHECK INDEXING
        A_mm1_nmk = solA_mm1_n(countK,2);               %%%% CHECK INDEXING
        
        % calculate Q_{m}^{n-k} from current period:
        %-------------------------------------------
        %Q_m_nmk = solQ_m(a_N-countK,2);                %%%% CHECK INDEXING
        Q_m_nmk = solQ_m_n(countK,2);                   %%%% CHECK INDEXING
        
        % calculate Q_{m-1}^{n-k} from current period:
        %---------------------------------------------
        %Q_mm1_nmk = solQ_mm1(a_N-countK,2);            %%%% CHECK INDEXING
        Q_mm1_nmk = solQ_mm1_n(countK,2);               %%%% CHECK INDEXING
        
    else
        % calculate A_m^{n-k} from previous period:
        %------------------------------------------
        %A_m_nmk = solA_m(a_N-countK,1);                %%%% CHECK INDEXING
        A_m_nmk = solA_m_n(countK,1);                   %%%% CHECK INDEXING
        
        % calculate A_{m-1}^{n-k} from previous period:
        %----------------------------------------------
        %A_mm1_nmk = solA_mm1(a_N-countK,1);            %%%% CHECK INDEXING
        A_mm1_nmk = solA_mm1_n(countK,1);               %%%% CHECK INDEXING
        
        % calculate Q_{m}^{n-k} from previous period:
        %--------------------------------------------
        %Q_m_nmk = solQ_m(a_N-countK,1);                %%%% CHECK INDEXING
        Q_m_nmk = solQ_m_n(countK,1);                   %%%% CHECK INDEXING
        
        % calculate Q_{m-1}^{n-k} from previous period:
        %----------------------------------------------
        %Q_mm1_nmk = solQ_mm1(a_N-countK,1); %%%% CHECK INDEXING
        Q_mm1_nmk = solQ_mm1_n(countK,1);      %%%% CHECK INDEXING
    end
    
    % calculate A_{m-1/2}^{n-k+1/2} using these values:
    %--------------------------------------------------
    A_mm_nmkp_half = 0.5*(A_m_nmk + A_mm1_nmk) - ...
        (delK/(2.0*delH))*(Q_m_nmk - Q_mm1_nmk);
    
    % calculate A_{m}^{n-k+1/2} using A_{m-1/2}^{n-k+1/2} and
    % A_{m+1/2}^{n-k+1/2} values:
    %--------------------------------------------------------
    A_m_nmkp_half    = 0.5*(A_mp_nmkp_half + A_mm_nmkp_half);
    P_m_nmkp_half    = getPressure(A_m_nmkp_half, R0m, P0);
    Qtms_m_np_half  = Qtms_m_np_half + ...
        P_m_nmkp_half*a_Admittances(countK)*delK;
    
    %%%% CHECK THE INDEXING OF a_Admittances
end

% calculate Q_{m-1/2}^{n+1/2} using discretization scheme:
%--------------------------------------------------------------------------
% get Q_m^n from existing solution:
%----------------------------------
Q_m_n   = solQ(a_N, end);

% get A_m^n from  existing solution:
%-----------------------------------
A_m_n   = solA(a_N, end);

% get R0_m_n:
%------------
X_m     = Xvec(end);
if ( isR0function )
    R0_m    = getR0function(X_m, Rtop, Rbottom, LenVessel);
    dR0dX_m = diffR0function(X_m, Rtop, Rbottom, LenVessel);
else
    R0_m    = getR0data(X_m);       %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_m = 0;                    %%%% HAVE NOT ENCODED THIS CASE
end

% get R(2)_m^n using these calculated and stored values:
%-------------------------------------------------------
R_m_n   = getFluxVector(A_m_n,Q_m_n,R0_m,P0);
R2_m_n  = R_m_n(2);

% get S(2)_m^n using these calculated and stored values:
%-------------------------------------------------------
S_m_n   = getSourceVector(A_m_n,Q_m_n,R0_m,dR0dX_m,P0);
S2_m_n  = S_m_n(2);

% get Q_{m-1}^n from existing solution:
%--------------------------------------
Q_mm1_n = solQ(a_N, end-1);

% get A_{m-1}^n from existing solution:
%--------------------------------------
A_mm1_n = solA(a_N, end-1);

% get R0_{m-21}^n:
%-----------------
X_mm1   = Xvec(end-1);
if ( isR0function )
    R0_mm1    = getR0function(X_mm1, Rtop, Rbottom, LenVessel);
    dR0dX_mm1 = diffR0function(X_mm1, Rtop, Rbottom, LenVessel);
else
    R0_mm1      = getR0data(X_mm1);    %%%% NEED TO IMPLEMENT THIS FUNCTION
    dR0dX_mm1   = 0.0;                 %%%% HAVE NOT ENCODED THIS CASE
end

% get R(2)_{m-1}^n using these calculated and stored values:
%-----------------------------------------------------------
R_mm1_n     = getFluxVector(A_mm1_n,Q_mm1_n,R0_mm1,P0);
R2_mm1_n    = R_mm1_n(2);

% get S(2)_{m-1}^n using these calculated and stored values:
%-----------------------------------------------------------
S_mm1_n     = getSourceVector(A_mm1_n,Q_mm1_n,R0_mm1,dR0dX_mm1,P0);
S2_mm1_n    = S_mm1_n(2);

% plug all these values into the discretization of Q_{m-1/2}^{n+1/2}:
% Q_{m-1/2}^{n+1/2} = \frac{Q_m^n + Q_{m-1}^n}{2} 
%           - \frac{k}{2}\frac{R(2)_m^n - R(2)_{m-1}^n}{h} 
%           + \frac{k}{2}\frac{S(2)_m^n + S(2)_{m-1}^n}{2}
%--------------------------------------------------------------------
Q_mm_np_half = 0.5*(Q_m_n + Q_mm1_n) - ...
    (delK/(2.0*delH))*(R2_m_n - R2_mm1_n) + ...
    (delK/4.0)*(S2_m_n + S2_mm1_n);

% use this value of (Qtms)_m^{n+1/2} and Q_{m-1/2}^{n+1/2} to calculate k1:
%--------------------------------------------------------------------------
%Qtms_m_np_half
%Q_mm_np_half
k1 = Qtms_m_np_half - 0.5*Q_mm_np_half;

% use the calculated values of k1, k2, k3 to calculate f1:
%---------------------------------------------------------
%k1
%Pm_k3_p_x2_half
%a_SolX(1)
f1 = k1 + Pm_k3_p_x2_half*k2 - a_SolX(1)/2.0;
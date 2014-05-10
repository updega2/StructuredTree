function f2 = getNonLinearEq2(a_SolX, a_N, a_Admittances)
%--------------------------------------------------------------------------
% Function to calculate the non-linear equation corresponding to the
% outflow boundary condition for variable A_{m+1/2}^{n+1/2} as defined by
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
global delK delH nsteps         % time-step size and mesh grid spacing
global isR0function             % property of the geometry

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
global P0

% Define the length of the vessel from X1 and X0:
%------------------------------------------------
LenVessel = X1 - X0;

% k_2 = y(X_m, t_n) which is the same as y_m^0:
%----------------------------------------------
k2 = a_Admittances(1)*delK;

% get the value of initial radius at X_m:
%----------------------------------------
Xm = Xvec(end);
if ( isR0function )
    R0_m = getR0function(Xm, Rtop, Rbottom, LenVessel);
else
    R0_m = getR0data(Xm);    %%%% NEED TO IMPLEMENT THIS FUNCTION
end

% calculate P(X_m, x(4)):
%------------------------
P_m         = getPressure(a_SolX(4),R0_m,P0);

Qtms_m_np1  = 0;

for countK = 2:nsteps
    
    %if ( a_N - countK + 0.5 >= 0 ) %%%% OLD VERSION
    if ( a_N - countK + 1.0 >= 0 ) %%%% TESTING SOMETHING
        
%         % calculate Q_{m+1/2}^{n-k+1/2} from current period:
%         %---------------------------------------------------
%         % OLD CODE:
%         %----------
%         %Q_mp_nmkp_half = solQ_mphalf(a_N-countK,2); %%%% 1. CHECK INDEXING
%         Q_mp_nmkp_half = solQ_mp_np_half(countK,2);  %%%% 2. CHECK INDEXING
%         
%         % calculate A_m^{n-k} from current period:
%         %-----------------------------------------
%         % OLD CODE:
%         %----------
%         %A_m_nmk = solA_m(a_N-countK,2);            %%%% 1. CHECK INDEXING
%         A_m_nmk = solA_m_n(countK,2);               %%%% 2. CHECK INDEXING
        
        % TESTING SOMETHING:
        %-------------------
        A_m_nmkp1 = solA_m_n(countK,2); % from current period
    else
        
%         % calculate Q_{m+1/2}^{n-k+1/2} from previous period:
%         %----------------------------------------------------
%         % OLD CODE:
%         %----------
%         %Q_mp_nmkp_half = solQ_mphalf(a_N-countK,1);%%%% 1. CHECK INDEXING
%         Q_mp_nmkp_half = solQ_mp_np_half(countK,1); %%%% 2. CHECK INDEXING
%         
%         % calculate A_m^{n-k} from previous period:
%         %------------------------------------------
%         % OLD CODE
%         %A_m_nmk = solA_m(a_N-countK,1);            %%%% 1. CHECK INDEXING
%         A_m_nmk = solA_m_n(countK,1);               %%%% 2. CHECK INDEXING
        
        % TESTING SOMETHING:
        %-------------------
        A_m_nmkp1 = solA_m_n(countK,1); % from previous period
    end
    
    if ( a_N - countK >= 0 )
        % calculate Q_m^{n-k} from current period:
        %-----------------------------------------
        %Q_m_nmk = solQ_m(a_N-countK,2);                %%%% 1. CHECK INDEXING
        Q_m_nmk = solQ_m_n(countK,2);                  %%%% 2. CHECK INDEXING
        
        % calculate Q_{m-1}^{n-k} from current period:
        %---------------------------------------------
        %Q_mm1_nmk = solQ_mm1(a_N-countK,2);            %%%% CHECK INDEXING
        Q_mm1_nmk = solQ_mm1_n(countK,2);               %%%% CHECK INDEXING
        
        % calculate A_m^{n-k} from current period:
        %-----------------------------------------
        %A_m_nmk = solA_m(a_N-countK,2);                %%%% CHECK INDEXING
        A_m_nmk = solA_m_n(countK,2);                   %%%% CHECK INDEXING
        
        % calculate A_{m-1}^{n-k} from current period:
        %---------------------------------------------
        %A_mm1_nmk = solA_mm1(a_N-countK,2);            %%%% CHECK INDEXING
        A_mm1_nmk = solA_mm1_n(countK,2);               %%%% CHECK INDEXING
    else
        % calculate Q_m^{n-k} from previous period:
        %------------------------------------------
        %Q_m_nmk = solQ_m(a_N-countK,1);                %%%% CHECK INDEXING
        Q_m_nmk = solQ_m_n(countK,1);                   %%%% CHECK INDEXING
        
        % calculate Q_{m-1}^{n-k} from previous period:
        %----------------------------------------------
        %Q_mm1_nmk = solQ_mm1(a_N-countK,1);            %%%% CHECK INDEXING
        Q_mm1_nmk = solQ_mm1_n(countK,1);               %%%% CHECK INDEXING
        
        % calculate A_m^{n-k} from previous period:
        %------------------------------------------
        %A_m_nmk = solA_m(a_N-countK,1);                %%%% CHECK INDEXING
        A_m_nmk = solA_m_n(countK,1);                   %%%% CHECK INDEXING
        
        % calculate A_{m-1}^{n-k} from previous period:
        %----------------------------------------------
        %A_mm1_nmk = solA_mm1(a_N-countK,1);            %%%% CHECK INDEXING
        A_mm1_nmk = solA_mm1_n(countK,1);               %%%% CHECK INDEXING
    end
    
%     % calculate R0(X_m):
%     %-------------------
%     X_m = Xvec(end);
%     if ( isR0function )
%         R0_m    = getR0function(X_m,...
%             Rtop, Rbottom, LenVessel);
%         dR0dX_m = diffR0function(X_m,...
%             Rtop, Rbottom, LenVessel);
%     else
%         R0_m    = getR0data(X_m);    %%%% NEED TO IMPLEMENT THIS FUNCTION
%         dR0dX_m = 0.0;               %%%% HAVE NOT ENCODED THIS
%     end
%     
%     % calculate R(2)_{m}^{n-k}:
%     %--------------------------
%     R_m_nmk     = getFluxVector(A_m_nmk,Q_m_nmk,R0_m,P0);
%     R2_m_nmk    = R_m_nmk(2);
%     
%     % calculate S(2)_{m}^{n-k}:
%     %--------------------------
%     S_m_nmk     = getSourceVector(A_m_nmk,Q_m_nmk,R0_m,dR0dX_m,P0);
%     S2_m_nmk    = S_m_nmk(2);
%     
%     % calculate R0(X_{m-1}):
%     %-----------------------
%     X_mm1 = Xvec(end-1);
%     if ( isR0function )
%         R0_mm1    = getR0function(X_mm1,Rtop, Rbottom, LenVessel);
%         dR0dX_mm1 = diffR0function(X_mm1,Rtop, Rbottom, LenVessel);
%     else
%         R0_mm1    = getR0data(X_mm1);  %%%% NEED TO IMPLEMENT THIS FUNCTION
%         dR0dX_mm1 = 0.0;               %%%% HAVE NOT ENCODED THIS
%     end
%     
%     % calculate R(2)_{m-1}^{n-k}:
%     %----------------------------
%     R_mm1_nmk     = getFluxVector(A_mm1_nmk,Q_mm1_nmk,R0_mm1,P0);
%     R2_mm1_nmk    = R_mm1_nmk(2);
%     
%     % calculate S(2)_{m-1}^{n-k}:
%     %----------------------------
%     S_mm1_nmk     = getSourceVector(A_mm1_nmk,Q_mm1_nmk,R0_mm1,dR0dX_mm1,P0);
%     S2_mm1_nmk    = S_mm1_nmk(2);
%     
%     % use these calculated values to plug into the discretisation for
%     % Q_{m-1/2}^{n-k+1/2}:
%     %----------------------------------------------------------------
%     Q_mm_nmkp_half = 0.5*(Q_m_nmk + Q_mm1_nmk) - ...
%         (delK/(2*delH))*(R2_m_nmk - R2_mm1_nmk) + ...
%         (delK/4)*(S2_m_nmk + S2_mm1_nmk);
%     
%     % use the value of Q_{m-1/2}^{n-k+1/2} to get A_m^{n-k+1}:
%     %---------------------------------------------------------
%     A_m_nmkp1 = A_m_nmk - (delK/4)*(Q_mp_nmkp_half - Q_mm_nmkp_half);
    
    % use this to get P(X_m, A_m^{n-k+1}):
    %-------------------------------------
    P_nmkp1 = getPressure(A_m_nmkp1, R0_m, P0);
    
    % use this to get (Q_{tms})_{m}^{n+1} = sum P(X_m,
    % A_m^{n-k+1})y_m^k\Delta t:
    %-------------------------------------------------
    Qtms_m_np1 = Qtms_m_np1 + P_nmkp1*a_Admittances(countK)*delK;
    
    %%%% CHECK THE INDEXING OF a_Admittances
end

k4 = Qtms_m_np1;

f2 = k4 + P_m*k2 - a_SolX(3);
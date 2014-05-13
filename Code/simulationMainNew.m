
%% One Dimensional Blood Flow simulation
%
% Main script for performing One-dimensional flow simulations using a
% structured tree outflow boundary condition.
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.

%% Configuring simulation inputs
%
% Ask the user to enter the input file where all global simulation data
% have been configured:
%--------------------------------------------------------------------------
% a. We start by defining all appropriate global variables required for
% the simulation
%--------------------------------------------------------------------------

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
global Rmin                     % minimum tree radius
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
global isMethodCharacteristics
global isCopyOlufsen
global isMaterialExpModel
global expModelK1 expModelK2 expModelK3

isMaterialExpModel = 1;
isCopyOlufsen = 0;
isMatlabNonLinSolve = 1;
isMethodCharacteristics = 0;

expModelK1 = 2.0*10^6;
expModelK2 = -2253.0;
expModelK3 = 8.65*10^4;

Imat = [1.0,0.0;0.0,1.0];

% NOTE: Data should be written into global solution arrays only from the
% main script. Data can be read from global solution arrays from any sub
% routine.

%--------------------------------------------------------------------------
% b. Then we ask user to input the simulation input dcleclclcccata file on the
% command window:
%--------------------------------------------------------------------------
inputFile = input('Enter the simulation input file \n','s');

%--------------------------------------------------------------------------
% c. Then from the input data file, all global data values are parsed:
%--------------------------------------------------------------------------
fileID = fopen(inputFile,'r');
if ( fileID == -1 )
    fprintf('Cannot read file\n')
end
while 1
    dataKey = fgetl(fileID);
    if (strcmpi(dataKey,'number-of-spatial-nodes'))
        nnodes = str2double(fgetl(fileID));
        fprintf('nnodes is %i \n', nnodes)
    elseif (strcmpi(dataKey,'number-of-temporal-nodes'))
        nsteps = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'cardiac-period'))
        Tcard = str2double(fgetl(fileID));
    elseif ( strcmpi(dataKey,'number-of-periods') )
        numPeriods = str2double(fgetl(fileID));
    elseif ( strcmpi(dataKey,'total-cardiac-output') )
        CardiacOut = str2double(fgetl(fileID));
    elseif ( strcmpi(dataKey,'instant-of-peak-cardiac-output') )
        CardiacPeak = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'domain-starting-point'))
        X0 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'domain-end-point'))
        X1 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'lower-frequency-value'))
        omega0 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'upper-frequency-value'))
        omega1 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'operating-pressure'))
        P0 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'artery-wall-elasticity'))
        Elast = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'artery-wall-thickness'))
        wallH = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'assumed-boundary-layer'))
        wallDelta = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'blood-viscosity'))
        Nu = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'blood-density'))
        rho = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'root-radius'))
        rootR = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'vessel-start-radius'))
        Rtop = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'vessel-end-radius'))
        Rbottom = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'initial-radius-is-a-function'))
        isR0function = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'iterative-error-tolerance'))
        errTol = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'minimum-radius'))
        Rmin = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'END-OF-FILE'))
        break
    end
end

%--------------------------------------------------------------------------
% d. Get Mu from density and viscosity:
%--------------------------------------------------------------------------
Mu = Nu * rho;

%% Create Spatial and Temporal Discretisations:
%-------------------------------------------------------------------------
Xvec = linspace(X0, X1, nnodes);    % The vector for x coordinates
Tvec = linspace(0, Tcard, nsteps);  % The vector for time in cardiac period
delK = Tvec(2) - Tvec(1);
delH = Xvec(2) - Xvec(1);

% NOTE: need to encode a function to check CFL condition for these

%% Generate a distribution of initial radius values
%--------------------------------------------------------------------------
LenVessel = X1-X0;
if ( isR0function )
    R0      = getR0function(Xvec, Rtop, Rbottom, LenVessel);
    dR0dX   = diffR0function(Xvec, Rtop, Rbottom, LenVessel);
else
    R0 = getR0data(Xvec);       % NEED TO IMPLEMENT THIS FUNCTION
end

A0 = pi*R0.*R0;

if ( isCopyOlufsen == 1 )
    m1     = 7.3693;
    m2     = 1.2122;
    m3     = 5.6517;
    m4     = 0.21763;
    for i = 1:nnodes
        womersley(i) = R0(i)*sqrt((2.0*pi/Tcard)/Nu);
        Cu(i)        = atan(m4*(-womersley(i)+m3))/m1 + m2;
    end
end

%% Fix a chosen value for flow-rate to initialize (let it be 1 for now):
%--------------------------------------------------------------------------
Qstart = 1*10^-6;

%% Initialize the solution vectors:
%--------------------------------------------------------------------------
Qvec = Qstart*ones(size(Xvec));     % The vector for flow rate solutions
Avec = A0;                          % The vector for area solutions

% a. Initialize the arrays for solution storage. Storage should be for one
% cardiac period's worth of spatial data
%--------------------------------------------------------------------------
% store _A_, _Q_ as nodes along columns, and time-steps along rows:
%--------------------------------------------------------------------------
solA = zeros(nsteps, nnodes);
solQ = zeros(nsteps, nnodes);
for timecounter = 1:nsteps
    solA(timecounter,:) = A0;
    solQ(timecounter,:) = Qstart;
end

solQ_mp_np_half     = zeros(nsteps, 2);
solA_mp_np_half     = zeros(nsteps, 2);
solQ_m_n            = zeros(nsteps, 2);
solA_m_n            = zeros(nsteps, 2);
solQ_mm1_n          = zeros(nsteps, 2);
solA_mm1_n          = zeros(nsteps, 2);
Aold                = zeros(1, nnodes);
Qold                = zeros(1, nnodes);

%% Build a tree for the outlet:
%--------------------------------------------------------------------------
flowTree = Tree(rootR);
flowTree.BuildTree(Rmin,10,0.5,0.5);

k = [(-nsteps/2):1:(nsteps/2)];
freqVec = ((2.*pi*k)/Tcard);

zTree_F = ones(size(k));


for i=1:length(k)
    indPlus = indexMap(k(i));
    if k(i) > 0
        zTree_F(indPlus-1) = flowTree.CalculateImpedance(freqVec(indPlus),1,0,50);
    end
end

for i=1:length(k)
    if k(i) <= 0
        indPlus = indexMap(-k(i));
        indMinus = indexMap(k(i));
        zTree_F(indMinus) = conj(zTree_F(indPlus));
    end
end

temp = zTree_F(nsteps/2);
tempTree = zTree_F;
for i = (nsteps/2):nsteps
    zTree_F(i+1) = tempTree(i);
end
zTree_F(nsteps/2) = temp;

if isCopyOlufsen
    zTree_T = real(ifft(fftshift(zTree_F)/Tcard,'symmetric'));
else
    zTree_T = real(ifft(zTree_F,'symmetric'));
end

yTree_T = 1./zTree_T;
%yTree_T = yTree_T;


%% Loop through the number of cardiac cycles:
%-------------------------------------------
simTime = 0;
for periodCount = 1:numPeriods
    
    % b. Store values of previous period for appropriate variables before
    % updating values for the current period:
    %---------------------------------------------------------------------
    if ( periodCount == 1 )
        % Simply initialize all arrays to a standard starting point:
        %-----------------------------------------------------------
        solA_m_n(:,:)        = A0(end);
        solQ_m_n(:,:)        = Qstart;
        solA_mm1_n(:,:)      = A0(end);
        solQ_mm1_n(:,:)      = Qstart;
        solA_mp_np_half(:,:) = A0(end);
        solQ_mp_np_half(:,:) = Qstart;
    else
        % Store all previous period solutions (from column 2) into column 1
        %------------------------------------------------------------------
        solA_m_n(:,1)           = solA_m_n(:,2);
        solQ_m_n(:,1)           = solQ_m_n(:,2);
        solA_mm1_n(:,1)         = solA_mm1_n(:,2);
        solQ_mm1_n(:,1)         = solQ_mm1_n(:,2);
        solA_mp_np_half(:,1)    = solA_mp_np_half(:,2);
        solQ_mp_np_half(:,1)    = solQ_mp_np_half(:,2);
    end
    
    % Add Initial Conditions For The Cardiac Period Solutions:
    %---------------------------------------------------------
    if ( periodCount == 1 )
        solQ(1,1) = getFlowRateIn(Tvec(1), CardiacOut, CardiacPeak, Tcard);
        
        Resistance = 33330500;    
        %Resistance = 333305; 
        %solQ(1,end) = getPressure(solA(1,end),R0(end),P0)/Resistance;
        solQ(1,end) = getPressure(solA(1,end),R0(end),P0)*yTree_T(1)*delK;
        solQ(1,end);
    else
        %------------------------------------------------------------------
        % get the first time instant (T=0) for the N+1'th period to match
        % the last time-instant (T=nsteps) for the N'th period to maintain
        % periodicity in time counting
        %------------------------------------------------------------------
        solA(1,:) = solA(nsteps,:);
        solQ(1,:) = solQ(nsteps,:);
    end
    
    
    % loop through the time-steps in a cardiac cycle. note that we are
    % going to calculate the solution at time = timeCounter, from the
    % solution at time = timeCounter - 1, in an explicit manner
    %----------------------------------------------------------------------
    for timeCounter = 2:nsteps
        fprintf('Time %i Period %i \n',timeCounter,periodCount);
   
        simTimeP = Tvec(timeCounter);
        simTime  = simTime + simTimeP;


        Aold = solA(timeCounter-1,:);
        Qold = solQ(timeCounter-1,:);
        
        % b. apply boundary condition at starting point:
        %-----------------------------------------------
        % b.1. calculate flux and source terms at X = X(2):
        %--------------------------------------------------
        R_1_n   = getFluxVector(Aold(2),Qold(2),R0(2),P0);
        S_1_n   = getSourceVector(Aold(2),Qold(2),R0(2),dR0dX(2),P0);
        
        % b.2. calculate flux and source terms at X = X(1):
        %--------------------------------------------------
        R_0_n   = getFluxVector(Aold(1),Qold(1),R0(1),P0);
        S_0_n   = getSourceVector(Aold(1),Qold(1),R0(1),dR0dX(1),P0);
        
        % b.3. update the solutions at the inlet:
        %----------------------------------------
        Qvec(1)         = getFlowRateIn(simTimeP,...
            CardiacOut, CardiacPeak, Tcard);
        Q_zero_nphalf   = getFlowRateIn(simTimeP + 0.5*delK,...
            CardiacOut, CardiacPeak, Tcard);
        
        if ( isCopyOlufsen == 1 )
            Q_half_nphalf   = 0.5*(Qold(2) + Qold(1)) - ...
                (delK/(2.0*Cu(1)*delH))*(R_1_n(2) - R_0_n(2)) + ...
                (delK/(4.0*Cu(1)))*(S_1_n(2) + S_0_n(2));
        else
            Q_half_nphalf   = 0.5*(Qold(2) + Qold(1)) - ...
                (delK/(2.0*delH))*(R_1_n(2) - R_0_n(2)) + ...
                (delK/4.0)*(S_1_n(2) + S_0_n(2));
        end
        
        Avec(1) = Aold(1) - (2.0*delK/delH)*(Q_half_nphalf - Q_zero_nphalf);
        
        
        % c. loop through the interior nodes:
        %------------------------------------
        
        for x = 2:(nnodes-1)
            
            %fprintf('Node %i Time %i Period %i \n',x,timeCounter,periodCount);
            
            U_X_T   = transpose([Aold(x), Qold(x)]);
            U_Xp1_T = transpose([Aold(x+1), Qold(x+1)]);
            U_Xm1_T = transpose([Aold(x-1), Qold(x-1)]);
            U_X_Tp1 = transpose([0,0]);
            
            % c.1. calculate R_{x+1}^t
            %---------------------------
            R_Xp1_T = getFluxVector(Aold(x+1),Qold(x+1),R0(x+1),P0);
            
            % c.2. calculate R_{x-1}^t:
            %--------------------------
            R_Xm1_T = getFluxVector(Aold(x-1),Qold(x-1),R0(x-1),P0);
            
            % c.3. calculate R_{x}^t:
            %------------------------
            R_X_T = getFluxVector(Aold(x),Qold(x),R0(x),P0);
            
            % c.4. calculate (dR0/dX)_{x+1}^t:
            %-----------------------------------
            if ( isR0function == 1 )
                dR0dX_Xp1_T = diffR0function(Xvec(x+1),...
                    Rtop, Rbottom, LenVessel);
            else
%                 if ( x == nnodes-1 )
%                     dR0dX_Xp1_T   = (R0(x+1) - R0(x))/delH;
%                 else
%                    dR0dX_Xp1_T   = (R0(x+2) - R0(x))/(2.0*delH);
%                 end
            end
            
            % c.5. calculate S_{x+1}^t:
            %--------------------------
            S_Xp1_T = getSourceVector(Aold(x+1),Qold(x+1),...
                R0(x+1),dR0dX_Xp1_T,P0);
            
            % c.6. calculate (dR0/dX)_{x-1/2}^t:
            %-----------------------------------
            if ( isR0function == 1 )
                dR0dX_Xm1_T = diffR0function(Xvec(x-1),...
                    Rtop, Rbottom, LenVessel);
            else
%                 if ( x == 2 )
%                     dR0dX_Xm1_T = (R0(2) - R0(1))/delH;
%                 else
%                     dR0dX_Xm1_T = (R0(x) - R0(x-2))/(2.0*delH);
%                 end
            end
            
            % c.7. calculate S_{x-1}^t:
            %--------------------------
            S_Xm1_T = getSourceVector(Aold(x-1),Qold(x-1),...
                R0(x-1),dR0dX_Xm1_T,P0);
            
            % c.9. calculate S_x^t:
            %----------------------
            if ( isR0function == 1 )
                dR0dX_X_T = diffR0function(Xvec(x),...
                    Rtop, Rbottom, LenVessel);
            else
%                 dR0dX_X_T = (R0(x+1) - R0(x-1))/(2.0*delH);
            end
            S_X_T = getSourceVector(Aold(x),Qold(x),R0(x),dR0dX_X_T,P0);
            
            % c.10. calculate U_{x+1/2}^{t+1/2}:
            %-----------------------------------
            if ( isCopyOlufsen == 1 )
                U_Xp_Tp_half(1) = 0.5*(U_Xp1_T(1) + U_X_T(1)) - ...
                    (0.5*delK/delH)*(R_Xp1_T(1) - R_X_T(1)) + ...
                    (delK/4.0)*(S_Xp1_T(1) + S_X_T(1));
                
                U_Xp_Tp_half(2) = 0.5*(U_Xp1_T(2) + U_X_T(2)) - ...
                    (0.5*delK/(delH*Cu(x)))*(R_Xp1_T(2) - R_X_T(2)) + ...
                    (delK/(4.0*Cu(x)))*(S_Xp1_T(2) + S_X_T(2));
            else
                
                U_Xp_Tp_half = 0.5*(U_Xp1_T + U_X_T) - ...
                    (0.5*delK/delH)*(R_Xp1_T - R_X_T) + ...
                    (delK/4.0)*(S_Xp1_T + S_X_T);
            end
            
            % c.11. calculate U_{x-1/2}^{t+1/2}:
            %-----------------------------------
            if ( isCopyOlufsen == 1 )
                U_Xm_Tp_half(1) = 0.5*(U_X_T(1) + U_Xm1_T(1)) - ...
                    (0.5*delK/delH)*(R_X_T(1) - R_Xm1_T(1)) + ...
                    (delK/4.0)*(S_X_T(1) + S_Xm1_T(1));
                
                U_Xm_Tp_half(2) = 0.5*(U_X_T(2) + U_Xm1_T(2)) - ...
                    (0.5*delK/(delH*Cu(x)))*(R_X_T(2) - R_Xm1_T(2)) + ...
                    (delK/(4.0*Cu(x)))*(S_X_T(2) + S_Xm1_T(2));
            else
                U_Xm_Tp_half = 0.5*(U_X_T+ U_Xm1_T) - ...
                    (0.5*delK/delH)*(R_X_T - R_Xm1_T) + ...
                    (delK/4.0)*(S_X_T + S_Xm1_T);
            end
            
            % assign $ Q_{x+1/2}^{t+1/2} $:
            %------------------------------
            Q_Xp_Tp_half = U_Xp_Tp_half(2);
            
            % assign A_{x+1/2}^{t+1/2}:
            %--------------------------
            A_Xp_Tp_half = U_Xp_Tp_half(1);
            
            % assign Q_{x-1/2}^{t+1/2}:
            %--------------------------
            Q_Xm_Tp_half = U_Xm_Tp_half(2);
            
            % assign A_{x-1/2}^{t+1/2}:
            %--------------------------
            A_Xm_Tp_half = U_Xm_Tp_half(1);
            
            % calculate R0_{x+1/2}:
            %----------------------
            if (isR0function)
                R0_Xp_half = getR0function(Xvec(x)+0.5*delH,Rtop,...
                    Rbottom,LenVessel);
            else
                R0_Xp_half = 0.5*(R0(x+1) + R0(x));
            end
            
            % calculate R0_{x-1/2}:
            %----------------------
            if (isR0function)
                R0_Xm_half = getR0function(Xvec(x)-0.5*delH,Rtop,...
                    Rbottom,LenVessel);
            else
                R0_Xm_half = 0.5*(R0(x-1) + R0(x));
            end
            % calculate R_{x+1/2}^{t+1/2}:
            %-----------------------------
            R_Xp_Tp_half = getFluxVector(A_Xp_Tp_half, Q_Xp_Tp_half, ...
                R0_Xp_half, P0);
            
            % calculate R_{x-1/2}^{t+1/2}:
            %-----------------------------
            R_Xm_Tp_half = getFluxVector(A_Xm_Tp_half, Q_Xm_Tp_half, ...
                R0_Xm_half, P0);
            
            % calculate (dR0/dX)_{x+1/2}:
            %----------------------------
            if (isR0function)
                dR0dX_Xp_half = diffR0function(Xvec(x)+0.5*delH,Rtop,...
                    Rbottom,LenVessel);
            else
                dR0dX_Xp_half = 0.5*(dR0dX_Xp1_T + dR0dX_X_T);
            end
            
            % calculate (dR0/dX)_{x-1/2}:
            %----------------------------
            if (isR0function)
                dR0dX_Xm_half = diffR0function(Xvec(x)-0.5*delH,Rtop,...
                    Rbottom,LenVessel);
            else
                dR0dX_Xm_half = 0.5*(dR0dX_Xm1_T + dR0dX_X_T);
            end

            
            % calculate S_{x+1/2}^{t+1/2}:
            %-----------------------------
            S_Xp_Tp_half = getSourceVector(A_Xp_Tp_half, Q_Xp_Tp_half, ...
                R0_Xp_half, dR0dX_Xp_half, P0);
            
            % calculate S_{x-1/2}^{t+1/2}:
            %-----------------------------
            S_Xm_Tp_half = getSourceVector(A_Xm_Tp_half, Q_Xm_Tp_half, ...
                R0_Xm_half, dR0dX_Xm_half, P0);
            
            % calculate the solution update:
            %-------------------------------
            if ( isCopyOlufsen == 1 )
                U_X_Tp1(1) = U_X_T(1) - ...
                    (delK/delH)*(R_Xp_Tp_half(1) - R_Xm_Tp_half(1)) + ...
                    (delK/2.0)*(S_Xp_Tp_half(1) + S_Xm_Tp_half(1));
                
                U_X_Tp1(2) = U_X_T(2) - ...
                    (delK/(delH*Cu(x)))*(R_Xp_Tp_half(2) - R_Xm_Tp_half(2)) + ...
                    (delK/(2.0*Cu(x)))*(S_Xp_Tp_half(2) + S_Xm_Tp_half(2));
            else
                U_X_Tp1 = U_X_T - ...
                    (delK/delH)*(R_Xp_Tp_half - R_Xm_Tp_half) + ...
                    (delK/2.0)*(S_Xp_Tp_half + S_Xm_Tp_half);
            end
            
            % store this solution update into solA and solQ (CHECK IMPLEMENTATION)
            
            Avec(x) = U_X_Tp1(1);
            
            Qvec(x) = U_X_Tp1(2);
            
            if (Avec(x) <0)
                Avec(x)
                error('Avec is negative');
            end
            
        end
        
        % d. apply boundary condition at the ending point
        %------------------------------------------------
        
        % d.3. continue the iterative solution:
        %--------------------------------------
        %plot(linspace(1,nnodes,nnodes),Qvec);
        nPeriod = timeCounter;
        if ( isMethodCharacteristics == 1 )
            [A,Q] = characteristicOutflow(nPeriod,yTree_T);
            %[A,Q] = characteristicResistance(nPeriod);
            Qvec(nnodes) = Q;
            Avec(nnodes) = A;
            Qvec;
            Avec;
            % e. update all nodal solution values in the global solution arrays
            %------------------------------------------------------------------
            solA(timeCounter,:) = Avec(:);
            solQ(timeCounter,:) = Qvec(:);
            
            solA_m_n(timeCounter,2)         = Avec(nnodes);
            solQ_m_n(timeCounter,2)         = Qvec(nnodes);
            solA_mm1_n(timeCounter,2)       = Avec(nnodes-1);
            solQ_mm1_n(timeCounter,2)       = Qvec(nnodes-1);
        else
            % d.2. create the vector for iterative solution guess:
            %-----------------------------------------------------
            Q_m_n   = Qold(end);
            A_m_n   = Aold(end);
            Q_mm1_n = Qold(end-1);
            A_mm1_n = Aold(end-1);
            if ( isR0function )
                R0_m_n      = getR0function(Xvec(end),...
                    Rtop, Rbottom, LenVessel);
                dR0dX_m_n   = diffR0function(Xvec(end),...
                    Rtop, Rbottom, LenVessel);
            else
                R0_m_n = getR0data(Xvec(end));
            end
            if ( isR0function )
                R0_mm1_n    = getR0function(Xvec(end-1),...
                    Rtop, Rbottom, LenVessel);
                dR0dX_mm1_n = diffR0function(Xvec(end-1),...
                    Rtop, Rbottom, LenVessel);
            else
                R0_mm1_n = getR0data(Xvec(end-1));
            end
            
            R_m_n   = getFluxVector(A_m_n, Q_m_n, R0_m_n, P0);
            S_m_n   = getSourceVector(A_m_n, Q_m_n, R0_m_n, dR0dX_m_n, P0);
            R1_m_n  = R_m_n(1);
            R2_m_n  = R_m_n(2);
            S1_m_n  = S_m_n(1);
            S2_m_n  = S_m_n(2);
            R_mm1_n = getFluxVector(A_mm1_n, Q_mm1_n,R0_mm1_n, P0);
            S_mm1_n = getSourceVector(A_mm1_n, Q_mm1_n, R0_mm1_n, dR0dX_mm1_n, P0);
            R1_mm1_n = R_mm1_n(1);
            R2_mm1_n = R_mm1_n(2);
            S1_mm1_n = S_mm1_n(1);
            S2_mm1_n = S_mm1_n(2);
            Q_mm_np_half = 0.5*(Q_m_n + Q_mm1_n) - ...
                (delK/(2*delH))*(R2_m_n - R2_mm1_n) + ...
                (delK/4)*(S2_m_n + S2_mm1_n);
            A_mm_np_half = 0.5*(A_m_n + A_mm1_n) - ...
                (delK/(2*delH))*(R1_m_n - R1_mm1_n) + ...
                (delK/4)*(S1_m_n + S1_mm1_n);
            %xIter0  = [Q_m_n, A_m_n, Q_mm_np_half, A_mm_np_half]';
            xIter0  = [Q_mm_np_half, A_mm_np_half, Q_m_n, A_m_n]';
            fX      = [0, 0, 0, 0]';
            xIter   = xIter0;
            
            errIter = 10*errTol;
            if ( isMatlabNonLinSolve == 1 )
                f1 = @(x) getNonLinearEq1(x, nPeriod, yTree_T);
                f2 = @(x) getNonLinearEq2(x, nPeriod, yTree_T);
                f3 = @(x) getNonLinearEq3(x, nPeriod, yTree_T);
                f4 = @(x) getNonLinearEq4(x, nPeriod, yTree_T);
                options = optimset('Display','off');
                [xNew, fVal] = fsolve(@(x)[f1(x), f2(x), f3(x), f4(x)], xIter0,options);
                fVal;
            else
                while (errIter > errTol)
                    fX(1)   = getNonLinearEq1(xIter, nPeriod, yTree_T);
                    fX(2)   = getNonLinearEq2(xIter, nPeriod, yTree_T);
                    fX(3)   = getNonLinearEq3(xIter, nPeriod, yTree_T);
                    fX(4)   = getNonLinearEq4(xIter, nPeriod, yTree_T);
                    Df      = getJacobianX(xIter, nPeriod, yTree_T);
                    xNew    = xIter - Df\fX;
                    errIter = norm(xNew - xIter,2)/norm(xNew - xIter0);
                    xIter   = xNew;
                end
            end
            % d.4. use the converged solution to update outlet values:
            %---------------------------------------------------------
            Qvec(end) = xIter(3);
            Avec(end) = xIter(4);
            % e. update all nodal solution values in the global solution arrays
            %------------------------------------------------------------------
            solA(timeCounter,:) = Avec(:);
            solQ(timeCounter,:) = Qvec(:);
            
            solA_m_n(timeCounter,2)         = Avec(end);
            solQ_m_n(timeCounter,2)         = Qvec(end);
            solA_mm1_n(timeCounter,2)       = Avec(end-1);
            solQ_mm1_n(timeCounter,2)       = Qvec(end-1);
            solA_mp_np_half(timeCounter,2)  = xIter(2);
            solQ_mp_np_half(timeCounter,2)  = xIter(1);
        end
        
        kcfl = checkCFL(solA(timeCounter,:),solQ(timeCounter,:));
        %plot(linspace(1,nnodes,nnodes),solA(timeCounter,:));
        %hold on;
        %plot(linspace(1,nnodes,nnodes),A0,'r')
    end
    
end
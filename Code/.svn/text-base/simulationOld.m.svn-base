%-------------------------------------------------------------------------%
% Main script for performing One-dimensional flow simulations using a
% structured tree outflow boundary condition.
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%-------------------------------------------------------------------------%
% 1. Ask the user to enter the input file where all global simulation data
% have been configured:
%--------------------------------------------------------------------------
global pi
pi = 22.0/7.0;

global nnodes
global nsteps
global T0 T1
global X0 X1
global omega0 omega1
global delK delH
global rootR isR0function
global Elast wallH wallDelta
global Nu rho

inputFile = input('Enter the simulation input file \n','s');

fileID = fopen(inputFile,'r');
if ( fileID == -1 )
    fprintf('Cannot read file\n')
end
while 1
    dataKey = fgetl(fileID);
    if (strcmpi(dataKey,'NUM-SPATIAL-NODES'))
        nnodes = str2double(fgetl(fileID));
        fprintf('nnodes is %i \n', nnodes)
    elseif (strcmpi(dataKey,'NUM-TEMPORAL-NODES'))
        nsteps = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'SIM-START-TIME'))
        T0 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'SIM-STOP-TIME'))
        T1 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'domain-starting-point'))
        X0 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'domain-end-point'))
        X1 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'lower-frequency-value'))
        omega0 = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'upper-frequency-value'))
        omega1 = str2double(fgetl(fileID));
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
    elseif (strcmpi(dataKey,'initial-radius-is-a-function'))
        isR0function = str2double(fgetl(fileID));
    elseif (strcmpi(dataKey,'END-OF-FILE'))
        break
    end
end

% 2. Create Spatial and Temporal Discretisations:
%--------------------------------------------------------------------------
delK = (T1 - T0)/nsteps;
delH = (X1 - X0)/nnodes;
% encode a function to check CFL condition for these

% 3. Initialize the solution vectors: (we have let them be 0 for now)
%--------------------------------------------------------------------------
Xvec = linspace(X1, X0, nnodes);   % The vector for x coordinates
Qvec = zeros('like',Xvec);          % The vector for flow rate solutions
Avec = zeros('like',Xvec);          % The vector for area solutions
%R0   = zeros('like',Xvec);          % The vector for initial radius values

% 4. Generate a distribution of initial radius values
%--------------------------------------------------------------------------
if ( isR0function )
    R0 = getR0function(Xvec);   % NEED TO IMPLEMENT THIS FUNCTION
else
    R0 = getR0data(Xvec);       % NEED TO IMPLEMENT THIS FUNCTION
end

timeCounter = 1;

Imat = [1.0,0.0;0.0,1.0];

% 5. Build a tree for the outlet:
%--------------------------------------------------------------------------
flowTree = Tree(rootR);
flowTree.BuildTree(0.3*rootR,10,0.5,0.5);
freqVec = linspace(omega0,omega1,100);

zTree_F = flowTree.CalculateImpedance(freqVec,1,0,1.055,0.049,0.046,50);
zTree_T = zeros('like',zTree_F);

while (simTime <= T1)
    
    % a. update simulation time
    simTime = T0 + delK*timeCounter;
    
    % b. apply boundary condition at starting point
    Avec(1) = pi*R0(1)*R0(1);
    Qvec(1) = getFlowrate(simTime); % NEED TO IMPLEMENT THIS FUNCTION
    
    % c. loop through the interior nodes
    for x = 2:(nnodes-1)
        U_X_T = [Avec(x),Qvec(x)]';
        
        % c.1. calculate R_{x+1}^t:
        R_Xp1_T = getFluxVector(Avec(x+1),Qvec(x+1),R0(x+1),P0);
        
        % c.2. calculate R_{x-1}^t:
        R_Xm1_T = getFluxVector(Avec(x-1),Qvec(x-1),R0(x-1),P0);
        
        % c.3. calculate R_{x}^t:
        R_X_T = getFluxVector(Avec(x),Qvec(x),R0(x),P0);
        
        % c.4. calculate M_{x+1}^t:
        M_Xp1_T = JacobianDRDU(Avec(x+1),Qvec(x+1),R0(x+1),P0);
        
        % c.5. calculate M_{x-1}^t:
        M_Xm1_T = JacobianDRDU(Avec(x-1),Qvec(x-1),R0(x-1),P0);
        
        % c.6. calculate M_{x}^t:
        M_X_T = JacobianDRDU(Avec(x),Qvec(x),R0(x),P0);
        
        % c.7. calculate S_{x+1}^t:
        if ( isR0function == 1 ) then
            dR0dX = diffR0function(Xvec(x+1));
        else
            if ( x == nnodes-1 )
                dR0dX   = (R0(x+1) - R0(x))/delH;
            else
                dR0dX   = (R0(x+2) - R0(x))/(2.0*delH);
            end
        end
        S_Xp1_T = getSourceVector(Avec(x+1),Qvec(x+1),R0(x+1),dR0dX,P0);
        
        % c.8. calculate S_{x-1}^t:
        if ( isR0function == 1 )
            dR0dX = diffR0function(Xvec(x-1));
        else
            if ( x == 2 )
                dR0dX = (R0(2) - R0(1))/delH;
            else
                dR0dX = (R0(x) - R0(x-2))/(2.0*delH);
            end
        end
        S_Xm1_T = getSourceVector(Avec(x-1),Qvec(x-1),R0(x-1),dR0dX,P0);
        
        % c.9. calculate S_x^t:
        if ( isR0function == 1 )
            dR0dX = diffR0function(Xvec(x));
        else
            dR0dX = (R0(x+1) - R0(x-1))/(2.0*delH);
        end
        S_X_T = getSourceVector(Avec(x),Qvec(x),R0(x),dR0dX,P0);
        
        % c.10. calculate the discretization update:
        term1 = (Imat - ...
            jacobianSource(Avec(x),Qvec(x),R0(x),dR0dX)*(delK/2.0))*Uvec;
        
        term2 = (delK - ((delK^2)/(4.0*delH))*(M_Xp1_T - M_Xm1_T))*...
            ((R_Xp1_T - R_Xm1_T)/(2.0*delH));
        
        term3 = (delK - ((delK^2)/(4.0*delH))*(M_Xp1_T - M_Xm1_T))*S_X_T;
        
        term4 = ((delK^2)/(8.0*delH^2))*M_X_T*...
            (R_Xp1_T - 2.0*R_X_T + R_Xm1_T);
        
        term5 = ((delK^2)/(4.0*delH))*M_X_T*(S_Xp1_T - S_Xm1_T);
        
        rhs = term1 - term2 + term3 + term4 - term5;
        Uvec = (Imat - jacobianSource(Avec(x),Qvec(x),R0(x),dR0dX,P0)\rhs;
        
        Avec(x) = Uvec(1);
        Qvec(x) = Uvec(2);
    end
    
    % d. apply boundary condition at the ending point
    % NEED TO IMPLEMENT THIS
    zTree_T = ifft(zTree_F);
    
end
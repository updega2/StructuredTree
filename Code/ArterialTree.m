% ---------------------------------------------- %
%                                                %
%  This is the Class Definition for each node    %
%  of a Arterial Tree defining the geometry      %
%  to be solved for the problem                  %
%  where the tree branches, we will need         %
%  information about radius and the two          %
%  child branches.                               %
%  Date: 5/9/2014                                %
%  Developers: Adam Updegrove, Alex Baelde       %
%              Debanjan Mukherjee                %
% ---------------------------------------------- %
classdef ArterialTree < handle
    %Define the member data of this structure
    properties
        GeomFile = 'geom.dat';
        Nodes = [0,0,0,0,0,0]; 
        %[daugther 1,daughter2,parent,radius start,radius end,length]
        
    end
    %Define the functions of this structure
    methods
        
        % This is a function to set all the necessary 
        % parameters for a model with a geom file (i.e. you create an
        % areterial tree by calling arteries =
        % ArterialTree(radius)
        function ReadGeomFile(obj)
            if nargin > 0
                fileID = fopen(obj.GeomFile,'r');
                if ( fileID == -1 )
                    fprintf('Geometry file called geom.dat not found\n')
                end

                while 1
                    dataKey = fgetl(fileID);
                    if (strncmpi(dataKey,'root',4))
                        obj.Nodes(1,3) = 0; % This is the root and doesn't
                                            % have a parent vessel
                        while 1
                            newLine = fgetl(fileID);
                            if (strncmpi(newLine,'radius',6))
                                rLine = sscanf(fgetl(fileID),'%f %f');
                                obj.Nodes(1,4) = rLine(1);
                                obj.Nodes(1,5) = rLine(2);
                            elseif (strncmpi(newLine,'length',6))
                                obj.Nodes(1,6) = str2double(fgetl(fileID));
                            elseif (strncmpi(newLine,'daughter',8))
                                dLine = sscanf(fgetl(fileID),'%d %d');
                                obj.Nodes(1,1) = dLine(1);
                                obj.Nodes(1,2) = dLine(2);
                            elseif (strncmpi(newLine,'terminal',8))
                                obj.Nodes(1,1) = 0;
                                obj.Nodes(1,2) = 0;
                            elseif (strncmpi(newLine,'end',3));
                                break
                            end
                        end
                    elseif (strncmpi(dataKey,'branch',5));
                        obj.Nodes(end+1,1) = 0; % There is a new row
                        while 1
                            newLine = fgetl(fileID);
                            if (strncmpi(newLine,'radius',6))
                                rLine = sscanf(fgetl(fileID),'%f %f');
                                obj.Nodes(end,4) = rLine(1);
                                obj.Nodes(end,5) = rLine(2);
                            elseif (strncmpi(newLine,'parent',6))
                                obj.Nodes(end,3) = ...
                                    str2double(fgetl(fileID));
                            elseif (strncmpi(newLine,'length',6))
                                obj.Nodes(end,6) = ...
                                    str2double(fgetl(fileID));
                            elseif (strncmpi(newLine,'daughter',8))
                                dLine = sscanf(fgetl(fileID),'%d %d');
                                obj.Nodes(end,1) = dLine(1);
                                obj.Nodes(end,2) = dLine(2);
                            elseif (strncmpi(newLine,'terminal',8))
                                obj.Nodes(end,1) = 0;
                                obj.Nodes(end,2) = 0;
                            elseif (strncmpi(newLine,'end',3));
                                break
                            end
                        end
                    elseif (strcmpi(dataKey,'END-OF-FILE'))
                        break
                    end
                end
            end
        end
        
        function BuildTree(obj,rmin,generations,alpha,beta)
            if nargin > 0
               iterations = 2^generations;
       
               obj.Nodes(1,3) = 0;
               obj.Nodes(1,4) = obj.R0;
            
               for i=1:iterations
                   
                   dims = size(obj.Nodes);
                   x = dims(1);
                   
                   if (x == i && i ~= 1)
                       return;
                   end
                   
                   r1 = obj.Nodes(i,4)*alpha;
                   r2 = obj.Nodes(i,4)*beta;
                   if r1 < rmin
                       %disp('Terminal Branch');
                   else
                       obj.Nodes(i,1) = x+1;
                       obj.Nodes(x+1,3) = i;
                       obj.Nodes(x+1,4) = r1;
                   end
                       
                   if r2 < rmin
                       %disp('Terminal Branch');
                   else
                       obj.Nodes(i,2) = x+2;
                       obj.Nodes(x+2,3) = i;
                       obj.Nodes(x+2,4) = r2;
                   end
                   if r2< rmin && r1< rmin
                       %disp('Both Branches are Terminal');
                   end
               end    
                   
            else
                error('Must provide input radius,minimum radius,alpha,and beta');
            end
        end

        %Function to Calculate the Effective Impedance as a function of
        %frequency
        function [Ztotal] = CalculateImpedance(obj,omega,x,Rterm,lratio)
            if nargin > 0
                %Define global variables
                global Mu Nu rho            % physical properties of the fluid
                
                %This is a recursive function to find the total impedance
                %of a tree. It is based on the impedance of the two child
                %branches. The function returns if both child vessles are smaller
                %than the prescribed minimum radius
                
                %These are all the values needed by the Impedance
                %computation
                r = obj.Nodes(x,4);
                %Eh is a function of root radius and some constants (Olufssen)
                Eh = obj.R0*((2*10^7)*exp(-1*obj.R0*22.53)+8.65*10^5);
                %Root area
                A = pi*(r)^2;
                %Wormersley number
                w0 = sqrt(((1i^3)*r^2*omega)/Nu);
                Fj = (2*besselj(1,w0))/(w0*besselj(0,w0));
                %Compliance
                C = (3*A*r)/(2*Eh);
                %Wave speed
                c = sqrt((A*(1-Fj))/(rho*C));
                %Length of vessel
                L = r*lratio;
                g=C*c;
                c1 = obj.Nodes(x,1);
                c2 = obj.Nodes(x,2);
                
                %If both children are zero, make them terminal impedance
                %and return
                if (c1 == 0 && c2 == 0)

                    Z1 = Rterm;
                    Z2 = Rterm;
                    
                    ZL = 1/((1/Z1)+(1/Z2));
                    
                    if omega ~= 0
                        
                        Zwtop = (1i*(1/g)*sin(omega*L/c))+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+1i*g*ZL*sin(omega*L/c);

                        Ztotal = Zwtop/Zwbottom;
                    else
                        Ztotal = (8*Mu*lratio)/(pi*obj.Nodes(x,4)^3)+ZL;
                    end
                %If one child is zero, make it terminal impedance
                elseif (c1 ~= 0 && c2 == 0)

                    Z1 = obj.CalculateImpedance(omega,obj.Nodes(x,1),Rterm,lratio);
                    Z2 = Rterm;
                    
                    ZL = 1/((1/Z1)+(1/Z2));
                    
                    if omega ~= 0
                        
                        Zwtop = (1i*(1/g)*sin(omega*L/c))+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+1i*g*ZL*sin(omega*L/c);

                        Ztotal = Zwtop/Zwbottom;
                    else
                        Ztotal = (8*Mu*lratio)/(pi*obj.Nodes(x,4)^3)+ZL;
                    end
                    
                %If other child is zero, make it terminal impedance
                elseif (c2 ~= 0 && c1 == 0)
                   
                    Z2 = obj.CalculateImpedance(omega,obj.Nodes(x,2),Rterm,lratio);
                    Z1 = Rterm;
                    
                    ZL = 1/((1/Z1)+(1/Z2));
                    
                    if omega ~= 0
                        
                        Zwtop = (1i*(1/g)*sin(omega*L/c))+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+1i*g*ZL*sin(omega*L/c);

                        Ztotal = Zwtop/Zwbottom;
                    else
                        Ztotal = (8*Mu*lratio)/(pi*obj.Nodes(x,4)^3)+ZL;
                    end
                    
                %Else if both children are not zero, calculate their
                %impedances
                else

                    Z1 = obj.CalculateImpedance(omega,obj.Nodes(x,1),Rterm,lratio);
                    Z2 = obj.CalculateImpedance(omega,obj.Nodes(x,2),Rterm,lratio);
                    
                    ZL = 1/((1/Z1)+(1/Z2));
                    
                    if omega ~= 0
                        
                        Zwtop = (1i*(1/g)*sin(omega*L/c))+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+1i*g*ZL*sin(omega*L/c);

                        Ztotal = Zwtop/Zwbottom;
                    else
                        Ztotal = (8*Mu*lratio)/(pi*obj.Nodes(x,4)^3)+ZL;
                    end
                end
            else
                error('Must provide frequency')
            end
        end
    end
end

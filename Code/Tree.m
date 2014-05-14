% ---------------------------------------------- %
%                                                %
%  This is the Class Definition for each node    %
%  of a Structured Tree Bed; at each point       %
%  where the tree branches, we will need         %
%  information about radius and the two          %
%  child branches.                               %
%  Date: 5/9/2014                                %
%  Developers: Adam Updegrove, Alex Baelde       %
%              Debanjan Mukherjee                %
% ---------------------------------------------- %
classdef Tree < handle
    %Define the member data of this structure
    properties
        R0        %Units(m)
        Nodes = [0 0 0 0] 
        Eh
    end
    %Define the functions of this structure
    methods
        
        % This is a constructor for the tree class (i.e. you create a
        % tree by calling tree0 =
        % Tree(radius)
        function obj = Tree(InputRadius)
            if nargin > 0
                obj.R0 = InputRadius;
            else
                error('Must provide Radius');
            end
        end
        
        function BuildTree(obj,rmin,generations,beta,alpha)
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
                       if r1 < rmin
                           obj.Nodes(i,2) = x+1;
                           obj.Nodes(x+1,3) = i;
                           obj.Nodes(x+1,4) = r2;
                       else
                           obj.Nodes(i,2) = x+2;
                           obj.Nodes(x+2,3) = i;
                           obj.Nodes(x+2,4) = r2;
                       end
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
                global Elast wallH
                global isMaterialExpSmall
                global expModelK1 expModelK2 expModelK3
      
                
                %This is a recursive function to find the total impedance
                %of a tree. It is based on the impedance of the two child
                %branches. The function returns if both child vessles are smaller
                %than the prescribed minimum radius
                
                %These are all the values needed by the Impedance
                %computation
                r = obj.Nodes(x,4);
          
                %Root area
                A = pi*(r)^2;
                %Wormersley number
                w0 = sqrt(((1i^3)*r^2*omega)/Nu);
                Fj = (2*besselj(1,w0))/(w0*besselj(0,w0));
                hold on;
                %Compliance
                if isMaterialExpSmall
                    Eh = r*(expModelK1*exp(expModelK2*r)...
                    + expModelK3);
                    C = (3*A*r)/(2*Eh);
                else
                    C = (3*A*r)/(2*Elast*wallH);
                end
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
                        Ztotal = (8*Mu*lratio)/(pi*r^3)+ZL;
                    end
                %If one child is zero, make it terminal impedance
                elseif (c1 ~= 0 && c2 == 0)

                    Z1 = obj.CalculateImpedance(omega,c1,Rterm,lratio);
                    Z2 = Rterm;
                    
                    ZL = 1/((1/Z1)+(1/Z2));
                    
                    if omega ~= 0
                        
                        Zwtop = (1i*(1/g)*sin(omega*L/c))+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+1i*g*ZL*sin(omega*L/c);

                        Ztotal = Zwtop/Zwbottom;
                    else
                        Ztotal = (8*Mu*lratio)/(pi*r^3)+ZL;
                    end
                    
                %If other child is zero, make it terminal impedance
                elseif (c2 ~= 0 && c1 == 0)
                   
                    Z2 = obj.CalculateImpedance(omega,c2,Rterm,lratio);
                    Z1 = Rterm;
                    
                    ZL = 1/((1/Z1)+(1/Z2));
                    
                    if omega ~= 0
                        
                        Zwtop = (1i*(1/g)*sin(omega*L/c))+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+1i*g*ZL*sin(omega*L/c);

                        Ztotal = Zwtop/Zwbottom;
                    else
                        Ztotal = (8*Mu*lratio)/(pi*r^3)+ZL;
                    end
                    
                %Else if both children are not zero, calculate their
                %impedances
                else

                    Z1 = obj.CalculateImpedance(omega,c1,Rterm,lratio);
                    Z2 = obj.CalculateImpedance(omega,c2,Rterm,lratio);
                    
                    ZL = 1/((1/Z1)+(1/Z2));
                    
                    if omega ~= 0
                        
                        Zwtop = (1i*(1/g)*sin(omega*L/c))+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+1i*g*ZL*sin(omega*L/c);

                        Ztotal = Zwtop/Zwbottom;
                    else
                        Ztotal = (8*Mu*lratio)/(pi*r^3)+ZL;
                    end
                end
            else
                error('Must provide frequency')
            end
        end
    end
end

        
                
            
        
    
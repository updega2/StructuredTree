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
classdef Node < handle
    %Define the member data of this structure
    properties
        Radius        %Units(m)
        Impedance=1
        Child1
        Child2
    end
    %Define the functions of this structure
    methods
        
        % This is a constructor for the branch class (i.e. you create a
        % branch by calling node0 =
        % Node(radius,child1,child2)
        function obj = Node(radius,child1,child2)
            if nargin >
                obj.Radius = radius;
                obj.Child1 = child1;
                obj.Child2 = child2;
            else
                error('Must provide Radius, and Child Branches');
            end
        end
        
        %Function to Calculate the Effective Impedance as a function of
        %frequency
        function CalculateImpedance(node,frequency,rmin,Rterm,rho,nu,E,h,lratio)
            if nargin > 0
                %This is a recursive function to find the total impedance
                %of a tree. It is based on the impedance of the two child
                %branches. The function returns if the vessel is smaller
                %than the prescribed minimum radius
                
                if (node.Radius < rmin)
                    node.Impedance = Rterm;
                    return;
                else
                    %Omega
                    omega = 2*pi*frequency;
                    %Root area
                    A = pi*(node.Radius^2);
                    %Wormersley number
                    w0 = sqrt(((i^3)*node.Impedance^2)/nu);
                    Fj = (2*besselj(1,w0))/(w0*besselj(0,w0));
                    %Compliance
                    C = (3*A*node.Radius)/(2*E*h);
                    %Wave speed
                    c = sqrt(A*(1-Fj))/(rho*C));
                    %Length of vessel
                    L = node.Radius*lratio;
                    g=C*c;

                    Y1 = node.Child1.CalculateImpedance(frequency,rmin,Rterm,rho,nu,E,h,lratio);
                    Y2 = node.Child2.CalculateImpedance(frequency,rmin,Rterm,rho,nu,E,h,lratio);
                    
                    ZL = 1/((1/node.Child1.Impedance)+(1/node.Child2.Impedance));
                    
                    if omega ~= 0
                        
                        Zwtop = (i*(1/g)*sin(omega*L/c)+ ZL*cos(omega*L/c);
                        Zwbottom = cos(omega*L/c)+i*(1.g)*ZL*sin(omega*L/c);

                        node.Impedance = Zwtop/Zwbottom;
                    else
                        node.Impedance = (8*mu*lratio)/(pi*node.Radius^3)+ZL;
                    end
                end
                   
            else
                error('Must provide frequency')
            end
        end
    end
end

        
                
            
        
    
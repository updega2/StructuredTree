% ---------------------------------------------- %
%                                                %
%  This is the Class Definition for each branch  %
%  of a Structured Tree Bed; at each point       %
%  where the tree branches, we will need         %
%  information about generation, index,          %
%  impedance, etc.                               %
%  Date: 5/9/2014                                %
%  Developers: Adam Updegrove, Alex Baelde       %
%              Debanjan Mukherjee                %
% ---------------------------------------------- %
classdef Branch < handle
    %Define the member data of this structure
    properties
        Generation    %This is the level in the tree
        Index         %This is the index in the generation
        Radius        %Units(m)
        LengthRatio 
        isEnd=false
        Impedance=1
    end
    %Define the functions of this structure
    methods
        
        % This is a constructor for the branch class (i.e. you create a
        % branch by calling branch0 =
        % Branch(generationvalue,indexvalue,radiusvalue,lengthratio,true/false)
        function obj = Branch(generation,index,radius,lengthratio,isend)
            if nargin > 0
                obj.Generation = generation;
                obj.Index = index;
                obj.Radius = radius;
                obj.LengthRatio = lengthratio;
                obj.isEnd = isend;
            else
                error('Must provide, Generation#, Index#, Radius, and Length Ratio');
            end
        end
        
        %Function to Calculate the Effective Impedance as a function of
        %frequency
        function CalculateImpedance(branch,frequency)
            if nargin > 0
                %Here will be the definition for this function
                branch.Impedance = 'NewValue';
                %
                %
                %
            else
                error('Must provide frequency')
            end
        end
        
        %Function to create a branch based on radius tolerance
        function MakeBranch(branch,tolerance)
            if nargin > 0
                %Here will be the definition for this function
                %
                %
                %
            else
                error('Must provide tolerance')
            end
        end
    end
end

        
                
            
        
    
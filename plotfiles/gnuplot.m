function [X,Y,Z] = gnuplot (A)
% GNUPLOT
%
%  [X,Y,Z] = gnuplot (A)
%
% Create 3D mesh from program output (3 columns) gnuplot style
% The first column must be x-values repeated for each different y.
% The second column must be the y-values.
% The third must be the corresponding z value.
% It is assumed that each block is separated with a blank line
% A can be input with load fname -ascii
% 
% Output can be used with mesh(X,Y,Z) and similar routines.

[height,width] = size(A);
if(width ~= 3)
  fprintf('Input matrix must have three columns\n');
  return;
end;

x = A(:,1);
y = A(:,2);
z = A(:,3); 

xx = x(1);               % Determine length of y-axis
ylength = 2;
while (xx == x(ylength)) 
  ylength = ylength + 1;
end;
ylength = ylength - 1;

xlength = height/ylength;          % Determine length of x-axis

X = zeros (xlength,ylength); 
Y = zeros (xlength,ylength);
Z = zeros (xlength,ylength);

for i=1:xlength 
  a = (i-1)*ylength+1;
  b = i*ylength;   
  X(i, :) = x(a:b)';
  Y(i, :) = y(a:b)';
  Z(i, :) = z(a:b)';
end


% GNUPLOT
%
%  [X,Y,Z1,Z2,Z3, Z4] = gnuplot (A)
%
% Create 3D mesh from program output (3 columns) gnuplot style
% The first column must be x-values repeated for each different y.
% The second column must be the y-values.
% The third must be the corresponding z value.
% It is assumed that each block is separated with a blank line
% A can be input with load fname -ascii
% 
% Output can be used with mesh(X,Y,Z) and similar routines.


function [X,Y,Z1,Z2,Z3,Z4] = gnuplot (A)

[height,width] = size(A);
if(width ~= 6)
  fprintf('Input matrix must have three columns\n');
  return;
end;

x = A(:,1);
y = A(:,2);
z1 = A(:,3);
z2 = A(:,4); 
z3 = A(:,5);
z4 = A(:,6);

xx = x(1);               % Determine length of y-axis
ylength = 2;
while (xx == x(ylength)) 
  ylength = ylength + 1;
end;
ylength = ylength - 1;

xlength = height/ylength;          % Determine length of x-axis

X = zeros (xlength,ylength); 
Y = zeros (xlength,ylength);
Z1 = zeros (xlength,ylength);
Z2 = zeros (xlength,ylength);
Z3 = zeros (xlength,ylength);
Z4 = zeros (xlength,ylength);

for i=1:xlength 
  a = (i-1)*ylength+1;
  b = i*ylength;   
  X(i, :) = x(a:b)';
  Y(i, :) = y(a:b)';
  Z1(i, :) = z1(a:b)';
  Z2(i, :) = z2(a:b)';
  Z3(i, :) = z3(a:b)';
  Z4(i, :) = z4(a:b)';
end


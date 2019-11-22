function out = plus(A,B) 
% This function overloads the default plus function, but only if the 
% default plus function can't solve it. Usefull to fix the function 
% magnifyOnFigure.m, in which the '+' operator is used on figure handles, 
% but after matlab 2014b these aren't seen as integers anymore. 
try 
out = double(A) + double(B); 
catch 
%Fail with the built-in 
builtin('plus',A,B); 
end
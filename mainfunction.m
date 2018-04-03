clear 
clc
close all


x=load('X.txt');  % import variables X, each column is one varible

y=load('y.txt');  %  import output y

index=SCR(x,y)    % index of identified contributing variables using contour variable selection algorithm

% index=SIR(x,y)  %  index of identified contributing variables using inverse variable selection algorithm








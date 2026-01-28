function [TH,Th,P1,P0] = getCoarseFineTriangulation(x_a,x_b,H_level,h_level)
%GETCOARSEFINETRIANGULATION Summary of this function goes here
%   Detailed explanation goes here

%% initial mesh
p = [x_a, x_a; x_b, x_a; x_a, x_b; x_b, x_b];
t = [1 2 4; 1 4 3];
Dnodes = [1;2;3;4];
T0 = struct('t',t,'p',p,'Dnodes',Dnodes);
[TH] = refineMesh(T0,H_level);
[Th,P1,~,P0,~] = refineMesh(TH,h_level-H_level);

end


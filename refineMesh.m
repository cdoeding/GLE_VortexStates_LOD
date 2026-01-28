function [Th,P1,R1,P0,P1dg] = refineMesh(Th,nref)
%REFINE_MESH Summary of this function goes here
%   Detailed explanation goes here

% number of vertices
NH = size(Th.p,1);
% compute fine reference mesh
P1 = speye(NH);
P0 = speye(size(Th.t,1));
P1dg = speye(numel(Th.t));
for k=1:nref
    [Th,p,p0,p1dg] = refine(Th);
    P1dg = p1dg*P1dg;
    P1 = p*P1;
    P0 = p0*P0;
end
% compute coarse mesh nodal interpolant
Nh = size(Th.p,1);
[i,j] = find(P1>1-2*eps);
R1 = sparse(j,i,ones(length(i),1),NH,Nh);
end


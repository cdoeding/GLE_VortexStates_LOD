function [e,nmbe] = computeEdges(T)
%COMPUTE_EDGES Summary of this function goes here
%   Detailed explanation goes here

np = size(T.p,1);
e = [T.t(:,[1,2]); T.t(:,[1,3]); T.t(:,[2,3])];
e = sort(e,2);
d2p = sparse(e(:,1),e(:,2),1,np,np);
[e1,e2,nmbe] = find(d2p);
e = [e1,e2];
end


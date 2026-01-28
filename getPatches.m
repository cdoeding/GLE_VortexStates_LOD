function patch = getPatches(TH,ell)
%GETPATCHES Summary of this function goes here
%   Detailed explanation goes here

d = 2;
NH = size(TH.p,1);
NTH = size(TH.t,1);

% incidence matrix of coarse elements and vertices
Ivt = sparse(TH.t,repmat((1:NTH).',1,d+1),1,NH,NTH);
% incidence matrix of coarse elements
Itt = spones(Ivt.'*Ivt);
% initialize patch matrix as identity
patch = speye(NTH);
% succesively apply Itt ell times
for k = 1:ell
    patch = Itt*patch;
end
% set non-zero entries to 1 and diagonal entries to 2
patch = spones(patch)+speye(NTH);
end


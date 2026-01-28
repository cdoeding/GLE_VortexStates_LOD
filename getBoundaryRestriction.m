function B = getBoundaryRestriction(B_H)
%GETBOUNDARYRESTRICTION Summary of this function goes here
%   Detailed explanation goes here

N_H = size(B_H,1);
active_nodes = diag(B_H);
N_H_active = sum(active_nodes);

% generate restriction operator from active_nodes
b_i = zeros(N_H_active,1);
b_j = zeros(N_H_active,1);
index = 1;
for i = 1:N_H
    if active_nodes(i)~=0
        b_i(index) = index;
        b_j(index) = i;
        index = index + 1;
    end
end

B = sparse(b_i,b_j,ones(N_H_active,1),N_H_active,N_H);

end


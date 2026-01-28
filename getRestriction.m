function [Rl_H,Rl_h] = getRestriction(T_H,T_h,l,patches,P0,B_H,B_h)
%GETRESTRICTION This function computes the
% coarse restriction operator Rl_H and the
% fine restriction operator Rl_h
% for the l-th local patch

N_h = size(T_h.p,1); % number of fine nodes
N_H = size(T_H.p,1); % number of coarse nodes

%% compute coarse restriction operator Rl_H
active_triangles = T_H.t(find(logical(patches(l,:))),:); % hold the triangles of the current patch (active triangles)

% mark active coarse nodes in vector active_nodes
active_nodes = false(N_H,1);
for i = 1:3
    for j = 1:size(active_triangles,1)
        active_nodes(active_triangles(j,i)) = true;
    end
end

active_nodes = logical(B_H*active_nodes); % unmark Dirichlet points in active coarse nodes
N_H_active = sum(active_nodes); % number of active coarse nodes

% generate restriction operator from active_nodes
r_i = zeros(N_H_active,1);
r_j = zeros(N_H_active,1);
index = 1;
for i = 1:N_H
    if active_nodes(i)~=0
        r_i(index) = index;
        r_j(index) = i;
        index = index + 1;
    end
end
Rl_H = sparse(r_i,r_j,ones(N_H_active,1),N_H_active,N_H);

%% compute Rl_h
active_fine_triangles = T_h.t(find(logical(P0*patches(l,:)')),:); % compute fien triangles in currect patch (fine active triangles)

% mark active fine nodes in vector active_fine_nodes
active_fine_nodes = zeros(N_h,1);
for i = 1:3
    for j = 1:size(active_fine_triangles,1)
        active_fine_nodes(active_fine_triangles(j,i)) = 1;
    end
end

active_fine_nodes = B_h*active_fine_nodes; % unmark Dirichlet points in active fine nodes
N_h_active = sum(active_fine_nodes);  % number of active fine nodes

% generate restriction operator from active_nodes
r_i = zeros(N_h_active,1);
r_j = zeros(N_h_active,1);
index = 1;
for i = 1:N_h
    if active_fine_nodes(i)~=0
        r_i(index) = index;
        r_j(index) = i;
        index = index + 1;
    end
end
Rl_h = sparse(r_i,r_j,ones(N_h_active,1),N_h_active,N_h);

end


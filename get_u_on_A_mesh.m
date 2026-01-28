function u_on_A = get_u_on_A_mesh(u,h_level,ht_level,P,M_h,M_ht)
%GET_U_ON_A_MESH Summary of this function goes here
%   Detailed explanation goes here

if h_level <= ht_level
    % polongation
    u_on_A = P'*u;

else
    % L2-projection
    u_on_A = M_ht'\(P*M_h'*u);

end


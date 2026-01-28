function U = evaluate_initial_value_LOD(u0,T,Nd,M_LOD,Phi,nodes2mesh_h)
%EVALUATE_INITAL_VALUE returns vector of initial value
% projected to finite element space

r = assemble_rhs(u0,T,Nd,nodes2mesh_h);
U = M_LOD'\(Phi*r);

end




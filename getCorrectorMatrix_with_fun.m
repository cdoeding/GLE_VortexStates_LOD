function Q = getCorrectorMatrix_with_fun(T_H,T_h,patches,A,kappa,beta,S_h,M_h,P1,P0,B_H,B_h,parallel)
%GETCORRECTORMA Summary of this function goes here
%   Detailed explanation goes here

if parallel
    Q = getCorrectorMatrixParallel_with_fun(T_H,T_h,patches,A,kappa,beta,S_h,M_h,P1,P0,B_H,B_h);
else
    Q = getCorrectorMatrixSequential_with_fun(T_H,T_h,patches,A,kappa,beta,S_h,M_h,P1,P0,B_H,B_h);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Q = getCorrectorMatrixSequential_with_fun(T_H,T_h,patches,A,kappa,beta,S_h,M_h,P1,P0,B_H,B_h)
N_h = size(T_h.p,1);
N_H = size(T_H.p,1);
NT_H = size(T_H.t,1);

C_h = P1'*M_h;
Q = sparse(N_H,N_h);

for l = 1:NT_H
    [Rl_H,Rl_h] = getRestriction(T_H,T_h,l,patches,P0,B_H,B_h);
    T = getNode2MeshMatrix(T_H,l);
    S_loc = assembleLocalBilinearForm_with_fun(A,kappa,T_h,P0,l);
    M_loc = assembleLocalMassMatrix(T_h,P0,l);
    
    Nl_h = size(Rl_h,1);
    Nl_H = size(Rl_H,1);
    
    S = Rl_h*(S_h + beta^2*M_h)*Rl_h';
    C = Rl_H*C_h*Rl_h';
    
    rhs = -T*B_H*P1'*(S_loc + beta^2*M_loc)*Rl_h';
    
    Y = S\C';
    S_inv = (C*Y)\speye(Nl_H,Nl_H);
    
    w = zeros(3,Nl_h);
    for i = 1:3
        q = S\rhs(i,:)';
        lambda = S_inv*C*q;
        w(i,:) = (q - Y*lambda)';
    end
    
    Q = Q + T'*sparse(w)*Rl_h;
    
    if mod(l,500)==0
        fprintf(strjoin([{num2str(l,'%i')},{'patches of'},{num2str(NT_H,'%i')},{'computed \n'}]))
    end
end

end



%Compute minimizers u of the (simplified) Ginzburg-Landau energy in 2d
% E(u) = \int_\Omega |(i/kappa)*grad(u)+Au|^2+(1/2)*(1-|u|^2)^2 dx
%
% with
% LOD discretization in order parameter u
% conjugate Sobolev gradient (with Polak–Ribière damping) flow discretization for minimization

%% initialize parameters
save = false;
plot = true;
parallel = false;

kappa = 10; % GL parameter
A = @(x) sqrt(2)*[sin(pi*x(1)).*cos(pi*x(2)); -sin(pi*x(2)).*cos(pi*x(1))];

iv_case = 1; % choose initial value here
alpha = (1/sqrt(2))*(1 + 1i);

switch iv_case
    case 1
        u0_ref = @(x) alpha*exp(-(x(1)^2+x(2)^2));
    case 2
        u0_ref = @(x) (x(1) + 1i*x(2))*alpha*exp(-(x(1)^2+x(2)^2));
    case 3
        u0_ref = @(x) 1 - alpha*exp(-5*(x(1)^2+x(2)^2));
    case 4
        u0_ref = @(x) 0.75*((x(1)+1i*x(2))/(sqrt(pi)))*exp(-(x(1)^2+x(2)^2)) + (0.25/(sqrt(pi)))*exp(-(x(1)^2+x(2)^2));
    case 5
        u0_ref = @(x) alpha*exp(10*1i*(x(1)^2+x(2)^2));
    case 6
        u0_ref = @(x) (x(1) + 1i*x(2))*alpha*exp(10*1i*(x(1)^2+x(2)^2));
    case 7
        u0_ref = @(x) 1 - alpha*sin(4*pi*x(1))^2*sin(4*pi*x(2))^2;
    case 8
        u0_ref = @(x) ((x(1) - 0.5) + 1i*(x(2) - 0.5))*alpha*exp(10*1i*((x(1) - 0.5)^2+(x(2) - 0.5)^2));
    case 9
        u0_ref = @(x) alpha*(1i + x(1) - 0.5);
    case 10
        u0_ref = @(x) alpha;
end

u0 = @(x) u0_ref(2*(x - 0.5*ones(size(x))));

x_a = 0; % domain left/bottom end point
x_b = 1; % domain right/top end point
area = 1; % area of rectangle

H_level = 7; % coarse mesh size level u (LOD)
h_level = 10; % fine mesh size level u (LOD)
ell = 10; % oversampling of ell-layers for u (LOD)

tol = 10^(-15); % tolerance for termination
i_max = 100000; % maximum number of iterations

save_path = "solution_kappa"+kappa+"_iv"+iv_case+".mat"; % path for save file

%% Eigenvalue parameter
lambda_tol = 10^(-10);
lambda_ev = zeros(10,10);
counter_ev = 1;
restart = true;

%% coarse, fine mesh and patches
boundary_u = 'Neumann';
[T_H,T_h,P1,P0] = getCoarseFineTriangulation(x_a,x_b,H_level,h_level);
T = T_H.t;
Nd = T_H.p;
P1 = P1';

% boundaries for H mesh (u equation)
B_H = getBoundaryNodes(T_H.p,x_a,x_b,'Neumann');
[nodes2mesh_H,~,~] = getNodes2Mesh(T_H.p,x_a,x_b,'Neumann','non-natural');

% boundaries for h mesh (u equation)
B_h = getBoundaryNodes(T_h.p,x_a,x_b,'Neumann');
[nodes2mesh_h,~,~] = getNodes2Mesh(T_h.p,x_a,x_b,'Neumann','non-natural');

% patches
patches = getPatches(T_H,ell); % patches_ij non-zero iff jth triangle is in patch of ith triangle

%% global mass and stiffness matricies
disp("start assemble global standard matricies")
S_h = assemble_bilinear_form_with_fun(A,kappa,T_h.t,T_h.p,parallel);
M_h = assemble_mass_matrix(T_h.t,T_h.p,nodes2mesh_h,parallel);
G_h = assemble_density_with_fun(A,T_h.t,T_h.p,parallel);

%% compute LOD basis (corrector)
disp("start computing corrector")
tic;
beta = 0;
Q = getCorrectorMatrix_with_fun(T_H,T_h,patches,A,kappa,beta,S_h,M_h,P1',P0,B_H,B_h,parallel);
time_Q = toc;

%% compute LOD matricies
S_LOD = (P1 + Q)*S_h*(P1 + Q)';
M_LOD = (P1 + Q)*M_h*(P1 + Q)';
G_LOD = (P1 + Q)*G_h*(P1 + Q)';

%% Matrices for eigenvalue computation
Phi = P1 + real(Q) + imag(Q);
N = size(M_LOD,1);
M_LOD_rc = sparse(2*N,2*N);
M_LOD_rc(1:N,1:N) = Phi*M_h*Phi';
M_LOD_rc(N+1:2*N,N+1:2*N) = Phi*M_h*Phi';

%% get initial value
u_LOD = evaluate_initial_value_LOD(u0,T_h.t,T_h.p,M_LOD,P1+Q,nodes2mesh_h);

% initial value in LOD space (L2-projection) and representation
%u_LOD = M_LOD'\((P1+Q)*M_h'*u_h);
u_h = (P1+Q)'*u_LOD;

%% corrector computation and assembling
Fu_h = assemble_nonlinear_term_impl(u_h,T_h.t,T_h.p,nodes2mesh_h,parallel);
Fu_LOD = (P1 + Q)*Fu_h*(P1 + Q)';

%% compute energy
E = 0.5*real(u_LOD'*S_LOD*u_LOD + 0.5*area - u_LOD'*M_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*u_LOD);

%% computation of Sobolev gradient flow
delta = tol;
counter = 0;
stab = 0;

Mat = S_LOD + Fu_LOD + stab*M_LOD + G_LOD;
grad_E = zeros(size(u_LOD));
d_u = grad_E;

while abs(delta) >= tol && counter < i_max

    E_old = E;
    Mat_old = Mat;
    grad_E_old = grad_E;
    d_u_old = d_u;

    % compute Sobolev gradient and descent direction
    Mat = S_LOD + Fu_LOD + stab*M_LOD + G_LOD;
    g_u = Mat\(((1+stab)*M_LOD + G_LOD)*u_LOD);
    grad_E = u_LOD - g_u;

    % compute Polak-Riebiere parameter
    if counter == 0 || delta == tol
        beta_PR_parameter = 0;
    else
        beta_PR_parameter =  real( ((grad_E-grad_E_old)'*Mat*grad_E)/(grad_E_old'*Mat_old*grad_E_old) );
        beta_PR_parameter = max( [ 0 , beta_PR_parameter ] );
    end

    % find optimal step size tau
    d_u = -grad_E + beta_PR_parameter*d_u_old;
    d_uh = (P1 + Q)'*d_u;

    [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1(u_h,d_uh,T_h.t,T_h.p,parallel);
    F_un_dn_LOD = (P1 + Q)*F_un_dn*(P1 + Q)';
    F_dn_dn_LOD = (P1 + Q)*F_dn_dn*(P1 + Q)';

    c0 = 0.5*real( u_LOD'*S_LOD*u_LOD ...
        + 0.5*area ...
        - u_LOD'*M_LOD*u_LOD  ...
        + 0.5*u_LOD'*Fu_LOD*u_LOD );

    c1 = 0.5*real( u_LOD'*S_LOD*d_u + d_u'*S_LOD*u_LOD...
        - d_u'*M_LOD*u_LOD - u_LOD'*M_LOD*d_u ...
        + 0.5*d_u'*Fu_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*d_u ...
        + u_LOD'*F_un_dn_LOD*u_LOD );

    c2 = 0.5*real( d_u'*S_LOD*d_u ...
        - d_u'*M_LOD*d_u ...
        + 0.5*u_LOD'*F_dn_dn_LOD*u_LOD ...
        + u_LOD'*F_un_dn_LOD*d_u ...
        + d_u'*F_un_dn_LOD*u_LOD ...
        + 0.5*d_u'*Fu_LOD*d_u );

    c3 = 0.5*real( d_u'*F_un_dn_LOD*d_u ...
        + 0.5*u_LOD'*F_dn_dn_LOD*d_u + 0.5*d_u'*F_dn_dn_LOD*u_LOD );

    c4 = 0.5*real( 0.5*d_u'*F_dn_dn_LOD*d_u );

    g_tau = @(t) ( c0 + c1*t + c2*t.^2 + c3*t.^3 + c4*t.^4 );

    upper_search_bound = 30;
    lower_search_bound = 0.1;
    [tau,g_tau_value] = golden_search_section(lower_search_bound,upper_search_bound,g_tau);

    % Sobolev gradient descent step
    u_LOD = u_LOD + tau*d_u;
    u_h = (P1 + Q)'*u_LOD;

    % update matrices and energy
    Fu_LOD = Fu_LOD + 2*tau*F_un_dn_LOD + (tau^2)*F_dn_dn_LOD;
    E = 0.5*real(u_LOD'*S_LOD*u_LOD + 0.5*area - u_LOD'*M_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*u_LOD);

    delta = E_old - E;
    counter = counter + 1;
    disp(counter)
    disp('step size')
    disp(tau)
    disp('energy difference')
    disp(delta)
    disp('energy')
    disp(E)

    %% check if limit is minimizer
    if abs(delta) < tol
        E_primeprime_LOD = assembleEnergyPrimePrimeMatrix_LOD(A,kappa,u_h,T_h.t,T_h.p,Phi);
        [V_LOD,lambda_LOD] = eigs(E_primeprime_LOD,M_LOD_rc,10,'smallestabs');
        lambda_ev(:,counter_ev) = diag(lambda_LOD);

        for k = 10:-1:2
            ev = lambda_ev(k,counter_ev);
            if (ev < 0 && abs(ev) >= lambda_tol) || (abs(ev) < lambda_tol)

                %% trigger restart & one step in EV direction
                d_uh = Phi'*(V_LOD(1:N,k) + 1i*V_LOD(N+1:2*N,k));
                d_u = M_LOD'\((P1+Q)*M_h'*d_uh);

                [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1(u_h,d_uh,T_h.t,T_h.p,parallel);
                F_un_dn_LOD = (P1 + Q)*F_un_dn*(P1 + Q)';
                F_dn_dn_LOD = (P1 + Q)*F_dn_dn*(P1 + Q)';

                c0 = 0.5*real( u_LOD'*S_LOD*u_LOD ...
                    + 0.5*area ...
                    - u_LOD'*M_LOD*u_LOD  ...
                    + 0.5*u_LOD'*Fu_LOD*u_LOD );

                c1 = 0.5*real( u_LOD'*S_LOD*d_u + d_u'*S_LOD*u_LOD...
                    - d_u'*M_LOD*u_LOD - u_LOD'*M_LOD*d_u ...
                    + 0.5*d_u'*Fu_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*d_u ...
                    + u_LOD'*F_un_dn_LOD*u_LOD );

                c2 = 0.5*real( d_u'*S_LOD*d_u ...
                    - d_u'*M_LOD*d_u ...
                    + 0.5*u_LOD'*F_dn_dn_LOD*u_LOD ...
                    + u_LOD'*F_un_dn_LOD*d_u ...
                    + d_u'*F_un_dn_LOD*u_LOD ...
                    + 0.5*d_u'*Fu_LOD*d_u );

                c3 = 0.5*real( d_u'*F_un_dn_LOD*d_u ...
                    + 0.5*u_LOD'*F_dn_dn_LOD*d_u + 0.5*d_u'*F_dn_dn_LOD*u_LOD );

                c4 = 0.5*real( 0.5*d_u'*F_dn_dn_LOD*d_u );

                g_tau = @(t) ( c0 + c1*t + c2*t.^2 + c3*t.^3 + c4*t.^4 );

                upper_search_bound = 30;
                lower_search_bound = 0.1;
                [tau,g_tau_value] = golden_search_section(lower_search_bound,upper_search_bound,g_tau);

                % Sobolev gradient descent step
                u_LOD = u_LOD + tau*d_u;
                u_h = (P1 + Q)'*u_LOD;

                % update matrices and energy
                Fu_LOD = Fu_LOD + 2*tau*F_un_dn_LOD + (tau^2)*F_dn_dn_LOD;
                E = 0.5*real(u_LOD'*S_LOD*u_LOD + 0.5*area - u_LOD'*M_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*u_LOD);

                counter_ev = counter_ev + 1;
                delta = tol;
                disp('--------------restart------------------')
                break
            end
        end
    end

end

%% plot
if plot == true
    model = createpde();
    geometryFromMesh(model,T_h.p',T_h.t');
    font = 28;
    pdeplot(model,"XYData",abs(u_h).^2, "ZData",abs(u_h).^2,'colormap',colormap(flipud(copper)),'Mesh','off','FaceAlpha', 1)
    view(0,90)
    clim([0 1])
    xlabel("$x$",'Interpreter', 'Latex', 'Fontsize', font)
    ylabel("$y$",'Interpreter', 'Latex', 'Fontsize', font)
    title("$|u_h|^2$", 'Interpreter', 'Latex', 'Fontsize', 20)
    ax = gca;
    ax.FontSize = font;
end

%% save
if save == true
    clearvars -except save_path H_level delta counter counter_ev tol E Q u_LOD u_h M_h T_h nodes2mesh_h parallel lambda_ev iv_case
    save(save_path,'-v7.3')
end

if parallel == true
    delete(gcp())
end





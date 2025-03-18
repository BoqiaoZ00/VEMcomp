function [u, v, t, uprime_norm, vprime_norm, u_average, v_average] = solver_parabolic_bulk_surf_block(dOmega, dGamma, f, g, h, P, M, MS, K, KS, R, T, tau, u0, v0)

n = length(dOmega);
m = length(dGamma);
PGamma = R'*P;
[omega_rows, omega_cols] = size(M);
[gamma_rows, gamma_cols] = size(MS);

% Precompute LHS matrices and preconditioners
LHS_Omega = cell(1, n);
for i = 1:n
    LHS_Omega{i} = M + tau * dOmega(i) * K; 
end
LHS_Gamma = cell(1, m);
for i = 1:m
    LHS_Gamma{i} = MS + tau * dGamma(i) * KS; 
end
LHS_Omega_block = blkdiag(LHS_Omega{:});
LHS_Gamma_block = blkdiag(LHS_Gamma{:});

NT = ceil(T / tau);
u = zeros(size(u0,1), size(u0,2));
u(:) = u0(:);
u_block = zeros(size(u,1)*size(u,2), 1);
u_block(:)=u(:);
v = zeros(size(v0,1), size(v0,2));
v(:) = v0(:);
v_block = zeros(size(v,1)*size(v,2), 1);
v_block(:)=v(:);
if nargout >= 4
    % spatial L2 norm of time derivative of fist component of u and v
    uprime_norm = zeros(1, NT);
    vprime_norm = zeros(1, NT);
end
if nargout >= 6
    % spatial average of first component of u and v
    u_average = zeros(1, NT);
    v_average = zeros(1, NT);
end
progress_handle = waitbar(0);
axes_handle = findobj(progress_handle, 'type','axes');
title_handle = get(axes_handle,'title');
set(title_handle, 'FontSize', 18);
waitbar(0, progress_handle, 'Timestepping in progress: 0 %')
percent_prev = 0;

tol = 1e-10;
maxit = 200;
precond_Omega = ichol(LHS_Omega_block, struct('type', 'ict', 'droptol', 1e-3)); 
precond_Gamma = ichol(LHS_Gamma_block, struct('type', 'ict', 'droptol', 1e-3)); 

RHS_Omega_block = zeros(size(M, 1) * n, 1);
RHS_Gamma_block = zeros(size(MS, 1) * m, 1);
RHS_Omega = zeros(size(M, 1), n);
RHS_Gamma = zeros(size(MS, 1), m);

for i=0:NT-1
   percent_new = round(i*100/NT);
   if percent_new > percent_prev
       waitbar(i/NT,progress_handle, sprintf('Timestepping in progress: %d%%', percent_new));
       percent_prev = percent_new;
   end

   u(:) = u_block(:);
   for j=1:n
       RHS_Omega(:,j)= M*(u(:,j)); % + tau*f{j}(u,P,i*tau)) + tau*R*MS*h{j}(R'*u,v, PGamma, i*tau); 
   end
   RHS_Omega_block(:) = RHS_Omega(:);

   v(:) = v_block(:);
   for j=1:m
       RHS_Gamma(:,j)= MS*(v(:,j)); % + tau*g{j}(R'*u,v,PGamma,i*tau));
   end
   RHS_Gamma_block(:) = RHS_Gamma(:);

    [unew_block, flag] = pcg(LHS_Omega_block, RHS_Omega_block, tol, maxit, precond_Omega, precond_Omega', u_block); 
    if flag ~= 0
        warning('PCG did not converge for Omega at time step %d', i);
    end
    [vnew_block, flag] = pcg(LHS_Gamma_block, RHS_Gamma_block, tol, maxit, precond_Gamma, precond_Gamma', v_block); 
    if flag ~= 0
        warning('PCG did not converge for Gamma at time step %d', i);
    end

   % if nargout >= 4
   %     incr_u = unew(:,1)-u(:,1);
   %     incr_v = vnew(:,1)-v(:,1);
   %     uprime_norm(i+1) = incr_u'*M*incr_u;
   %     vprime_norm(i+1) = incr_v'*MS*incr_v;
   % end
   % if nargout >= 6
   %     u_average(i+1) = sum(M*unew(:,1));
   %     v_average(i+1) = sum(MS*vnew(:,1));
   % end

   u_block = unew_block;
   v_block = vnew_block;
end

u(:) = u_block(:);
v(:) = v_block(:);

close(progress_handle);

if nargout >= 3
    t = linspace(tau, NT*tau, NT);
end
if nargout >= 4
    uprime_norm = sqrt(uprime_norm);
    vprime_norm = sqrt(vprime_norm);
end
if nargout >= 6
    u_average = u_average/sum(sum(M));
    v_average = v_average/sum(sum(MS));
end

end
% Solves Laplace equation on the unit cube through VEM with k=1
%
%   - Delta(u) + u = (3*pi^2+1)*cos(pi*x)*cos(pi*y)*cos(pi*z)
%
% Exact solution: u(x,y,z) = cos(pi*x)*cos(pi*y)*cos(pi*z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
% - Nx: number of nodes along each direction
% - plotflag: indicates whether generate plots
%
% OUTPUTS:
%
% - l2err: absolute error in L2 norm
% - l2err_rel: relative error in L2 norm
% - h: meshsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l2err,l2err_rel, h] = laplace3dcube(Nx,plotflag)

hx = 1/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
x = linspace(0,1,Nx);
P = zeros(Nx^3,3);
for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(i+1) x(j+1) x(k+1)];
        end
    end
end
  
P1 = [0 0 0; 0 1 0; 1 1 0; 1 0 0]; % bottom face
P2 = [0 0 1; 0 1 1; 1 1 1; 1 0 1]; % top face
P3 = [0 0 0; 0 1 0; 0 1 1; 0 0 1]; % back face
P4 = [1 0 0; 1 1 0; 1 1 1; 1 0 1]; % front face
P5 = [0 0 0; 1 0 0; 1 0 1; 0 0 1]; % left face
P6 = [0 1 0; 1 1 0; 1 1 1; 0 1 1]; % right face

E1S = element2dsquare(P1);
E2S = element2dsquare(P2);
E3S = element2dsquare(P3);
E4S = element2dsquare(P4);
E5S = element2dsquare(P5);
E6S = element2dsquare(P6);

PS = unique([P1; P2; P3; P4; P5; P6],'rows');
ES = element3dcube([E1S;E2S;E3S;E4S;E5S;E6S], PS);
KE = ES.K;
ME = ES.M;

K = spalloc(Nx^3,Nx^3,57*Nx^3);
M = spalloc(Nx^3,Nx^3,57*Nx^3);
for i=0:Nx-2
    for j=0:Nx-2
        for k=0:Nx-2
            indexes = [Nx^2*i+Nx*j+k+1
                       Nx^2*i+Nx*j+k+2
                       Nx^2*i+Nx*(j+1)+k+1
                       Nx^2*i+Nx*(j+1)+k+2
                       Nx^2*(i+1)+Nx*j+k+1
                       Nx^2*(i+1)+Nx*j+k+2
                       Nx^2*(i+1)+Nx*(j+1)+k+1
                       Nx^2*(i+1)+Nx*(j+1)+k+2];
            K(indexes, indexes) = K(indexes, indexes) + KE; %#ok
            M(indexes, indexes) = M(indexes, indexes) + ME; %#ok
        end
    end
end
K = K*hx;
M = M*hx^3;

f = @(x,y,z) (3*pi^2+1)*cos(pi*x).*cos(pi*y).*cos(pi*z);
esol = @(x,y,z) cos(pi*x).*cos(pi*y).*cos(pi*z);

% Solving linear system (K+M)u = Mf
tic
u = (K+M)\(M*f(P(:,1),P(:,2),P(:,3)));
es = esol(P(:,1),P(:,2),P(:,3));
toc

% Computing L2 error
err = es - u;
l2err = sqrt(err'*M*err);
l2solnorm = sqrt(es'*M*es);
l2err_rel = l2err/l2solnorm;

if plotflag
    
% Plotting Exact solution
figure
set(gcf,'Color','white')
scatter3(P(:,1), P(:,2), P(:,3), 30, esol(P(:,1), P(:,2), P(:,3)),'filled')
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
title('Exact solution')
colorbar

% Plotting Numerical solution
figure
set(gcf,'Color','white')
scatter3(P(:,1), P(:,2), P(:,3), 30, u,'filled')
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
title('Numerical solution')
colorbar

end

end
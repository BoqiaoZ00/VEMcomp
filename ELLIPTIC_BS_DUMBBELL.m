% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all

load('dumbbell.mat')
N = length(P); % Overall amount of nodes
NGamma = length(MS); % Amount of boundary nodes

% TOY PROBLEM 2 - 2 EIGENMODES
esol_u = @(P) P(:,1).^4 - P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 0.1;
esol_v = @(P) 0*P(:,1) + 1;
f_u = @(P) P(:,1).^4 - 13*P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 2.1;
f_v = @(P) 0*P(:,1) + 1;
flux_u = @(P) sqrt((4*P(:,1).^3 - 2*P(:,1)).^2 +4*(P(:,2).^2 + P(:,3).^2));

RM = spalloc(N, NGamma, NGamma);
boundarynode = unique(EGamma(:));
RM(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

tic
numsol = [K+M, 0*RM*MS; 0*MS*RM', KS+MS]\[M*f_u(P) + RM*MS*flux_u(P(boundarynode,:)); MS*f_v(P(boundarynode,:))];
u = numsol(1:N,1);
v = numsol(N+1:end,1);
toc
es_u = esol_u(P);
es_v = esol_v(P(boundarynode,:));
err_u = u - es_u;
err_v = v - es_v;
% L2err_u = sqrt(err_u'*M*err_u);
% L2err_v = sqrt(err_v'*MS*err_v);
% H1err_u = sqrt(err_u'*K*err_u);
% H1err_v = sqrt(err_v'*KS*err_v);

L2err_product_abs = sqrt(err_u'*M*err_u + err_v'*MS*err_v);
H1err_product_abs = sqrt(err_u'*(K+M)*err_u + err_v'*(KS+MS)*err_v);
L2norm_product = sqrt(es_u'*M*es_u + es_v'*MS*es_v);
H1norm_product = sqrt(es_u'*(K+M)*es_u + es_v'*(KS+MS)*es_v);
L2err_product_rel = L2err_product_abs/L2norm_product;
H1err_product_rel = H1err_product_abs/H1norm_product;


%xcut = range(1,1) + h/sqrt(3)*2;
xcut = 0;
Rcirc = sqrt(.1);      % Radius of circle 

% Plotting Numerical Solution - Bulk Component u
figure
indsol = P(:,1) >= xcut;
chull = convhull(P(indsol,:));
set(gcf, 'Color','white')
trisurf(chull, P(indsol,1), P(indsol,2), P(indsol,3), u(indsol,1), 'FaceColor', 'interp', 'EdgeColor','none')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('$U$', 'interpreter', 'latex')
colorbar
colormap jet
hold on
Ccirc = [xcut, 0, 0];   % Center of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud

% Plotting Numerical Solution - Surface Component v
figure
set(gcf,'Color','white')
trisurf(EGamma, P(:,1), P(:,2), P(:,3), RM*v, 'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
xlim([0,1])
title('$V$', 'interpreter', 'latex')
colorbar
colormap jet
hold on
Ccirc = [xcut, 0, 0];   % Center of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud

% Plotting Numerical Error - Bulk Component u
figure
indsol = P(:,1) >= xcut;
chull = convhull(P(indsol,1), P(indsol,2), P(indsol,3));
set(gcf, 'Color','white')
trisurf(chull, P(indsol,1), P(indsol,2), P(indsol,3), es_u(indsol,1)-u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
caxis([min(es_u-u), max(es_u-u)])
title('$u^{-\ell}-U$', 'interpreter', 'latex')
colorbar
colormap jet
hold on
Ccirc = [xcut, 0, 0];   % Center of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud

% Plotting Numerical Error - Surface Component v
figure
set(gcf,'Color','white')
trisurf(EGamma, P(:,1), P(:,2), P(:,3), RM*(es_v-v), 'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
caxis([min(es_v-v), max(es_v-v)])
xlim([0,1])
title('$v^{-\ell}-V$','interpreter', 'latex')
colorbar
colormap jet
hold on
Ccirc = [xcut, 0, 0];   % Center of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud
% calculate buttressing on BEDMAP2 grid
% using a constant viscosity parameter or
% a perturbation on the inverted viscosity
% INPUT:        see scripts in ..
%               meta_setup.m
%
% OUTPUT:       NONE
%
% by Johannes J. Fuerst

% =====================================================================
% VARIABLES
% =====================================================================

% =====================================================================
% CALCULATE VELOCITY GRADIENTS
% =====================================================================

small = 1E-10;

input.dudx      = input.x*NaN;
input.dvdx      = input.x*NaN;
input.dudy      = input.x*NaN;
input.dvdy      = input.x*NaN;
input.KN        = input.x*NaN;
input.KT        = input.x*NaN;
rf              = input.x*0;

rf           = rf + AA;

nx = size(input.x,1);
ny = size(input.x,2);
dx  = (input.x(1,2)-input.x(1,1));
dy  = (input.y(2,1)-input.y(1,1));
input.dudx(2:nx-1,2:ny-1) = (input.vx(2:nx-1,3:ny)-input.vx(2:nx-1,1:ny-2))./2./(input.x(1,2)-input.x(1,1));
input.dvdx(2:nx-1,2:ny-1) = (input.vy(2:nx-1,3:ny)-input.vy(2:nx-1,1:ny-2))./2./(input.x(1,2)-input.x(1,1));

input.dudy(2:nx-1,2:ny-1) = (input.vx(3:nx,2:ny-1)-input.vx(1:nx-2,2:ny-1))./2./(input.x(1,2)-input.x(1,1));
input.dvdy(2:nx-1,2:ny-1) = (input.vy(3:nx,2:ny-1)-input.vy(1:nx-2,2:ny-1))./2./(input.x(1,2)-input.x(1,1));

exx          = input.dudx;
eyy          = input.dvdy;
exy          = 0.5*(input.dudy+input.dvdx);
ezz          = -exx-eyy;
eff          = 0.5*(exx.*exx+eyy.*eyy+ezz.*ezz)+exy.*exy+small;
%eff          = 1.*(exx.*exx+eyy.*eyy+ezz.*ezz)+small;
%eff          = 1.*(exx.*exx+eyy.*eyy+exx.*eyy+exy.*exy)+small;

ttxx = (input.dudx).*rf.^(-1./nn).*eff.^((1.-nn)/(2.*nn));
ttyy = (input.dvdy).*rf.^(-1./nn).*eff.^((1.-nn)/(2.*nn));
ttxy = 0.5*(input.dudy+input.dvdx).*rf.^(-1./nn).*eff.^((1.-nn)/(2.*nn));
ttyx = ttxy;


% the following tensor represents the full stress distribution
% being a linear combination of the deviatoric stresses
%
%% THIS TENSOR IS DECISIVE FOR HORIZONTAL EXTENSION AND COMPRESSION
%
% om_xx = 2 tau_xx + 1 tau_yy
% om_yy = 2 tau_yy + 1 tau_xx
% om_xy = tau_xy
% om_yx = tau_yx
%
% the additional pressure is an offset to the eigenvalues but does not
% change the direction of the horizontal eigenvectors



t11    = 2.*ttxx+1.*ttyy;
t12    = ttxy;
t21    = ttxy;
t22    = 1.*ttxx+2.*ttyy;


% DETERMINE EIGENVALUES
% characteristic polynom
%
% aa * lambda^2 + bb*lambda + cc = 0;

aa           = ones(size(t11));
bb           = (-1)*(t11+t22);
cc           = t11.*t22-t21.*t12;

lambda1      = (-bb+sqrt(bb.^2-4*aa.*cc))./(2.*aa);
lambda2      = (-bb-sqrt(bb.^2-4*aa.*cc))./(2.*aa);

eig1   = lambda1;
eig2   = lambda2;

% DATERMINE EIGEN VECTORS v1 v2
% t11 * v11 + t12 * v12 = lambda1 * v11
% t21 * v21 + t22 * v22 = lambda2 * v21
%
% v12 = (lambda1 - t11)/t12 v11
% v21 = (lambda2 - t22)/t21 v22

% ? ERROR FROM HERE ON 
%
% ? EIGENVALUES ARE GOOD

eigv11 = ones(size(t11));
eigv12 = (lambda1 - t11).*eigv11./t12;

eigv11 = eigv11./sqrt(eigv11.^2+eigv12.^2);
eigv12 = eigv12./sqrt(ones(size(input.x)).^2+eigv12.^2);

eigv21 = ones(size(t11));
eigv22 = (lambda2 - t11).*eigv21./t12;

eigv21 = eigv21./sqrt(eigv21.^2+eigv22.^2);
eigv22 = eigv22./sqrt(ones(size(input.x)).^2+eigv22.^2);


%% CHOICE OF BUTTRESSING DIRECTION
% 1. FLOW LINES
% 2. EIGENVECTOR v1
% 3. EIGENVECTOR v2

n1          = input.vx./sqrt(input.vx.^2+input.vy.^2);
n2          = input.vy./sqrt(input.vx.^2+input.vy.^2);

tau    = ( (t11.*n1+t12.*n2).*n1 + (t21.*n1+t22.*n2).*n2 );
psi    = ( (t11.*n2-t12.*n1).*n2 - (t21.*n2-t22.*n1).*n1 );
  
tau0  = 0.5*rhoice*grav*thiampl*input.thi*(1-rhoice/rhowater);

KT    = psi./tau0;
KN1   = 1-tau./tau0;
KN2   = 1-eig1./tau0;
KN3   = 1-eig2./tau0;

clear  n1 n2 t11 t21 t12 t22 aa bb cc lambda1 lambda2 psi tau ttxx ttyy ttxy ttyx eigv11 eigv12 eigv21 eigv22 eig1 eig2 lambda1 lambda2
clear dx dy exx eyy exy eyx ezz kk nx ny rf tau0



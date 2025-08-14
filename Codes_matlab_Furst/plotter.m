% Matlab Plotting Script for PIG
%
% script by Johannes J. Fuerst, 2013
%
% INPUT:        run first meta_setup.m
%
% OUTPUT:       NONE
%

path(path,'../bin/')

% =====================================================================
% CONSTANTS
% =====================================================================

nexp   = 24;
scrsz  = get(0,'screensize');

% =====================================================================
% ELMER SHELF GRID + CALVING EXTENTS
% =====================================================================

%{
figure
set(gcf,'name','calving extents','position',[1 scrsz(4)/1.5 scrsz(3)/2.5 scrsz(4)/1.5]);

colourise = [1 1 0; 1 0 0; 0 1 0; 0.2 0.5 1; 0 0 0];
% assemble
  for n = 1:5
    plot(c(n).x.*c(n).shelfmask/1e3,c(n).y.*c(n).shelfmask/1e3,'o','MarkerFaceColor',squeeze(colourise(n,:)),'MarkerEdgeColor','none')
    axis equal
    hold on
  end
set(gca,'FontSize',20,'FontName','Myriad Pro');
legend
xlabel('zonal extent [km]')
ylabel('meridional extent [km] ')
title('calving extents [m] ')
%}

% =====================================================================
% ELMER SHELF GRID + CALVING EXTENTS
% =====================================================================


for n=2:1:5%1:nexp

figure
set(gcf,'name','buttressing','position',[1 scrsz(4)/1.5 scrsz(3)/2.11 scrsz(4)/1.5]);
axes_font

%% LAYER 1
% define plotted fields
cvec2 = [-10 -2 -0.5 0:0.1:1 2 4 100];
xx = input.x/1e3;
yy = input.y/1e3;
zz = squeeze(input.KN3(n,:,:)).*squeeze(input.shelfELM(n,:,:))./squeeze(sign(input.vmELM(n,:,:)));
if(n>2)
zeta = squeeze(input.shelfELM(n-1,:,:));
zeta(isnan(zeta)==1) = 0;
end
zeta2 = squeeze(input.shelfELM(n+23,:,:));
zeta2(isnan(zeta2)==1) = 0;

% define colorbar
clear cmap
cmap=bluebrown(length(cvec2),3);
colormap(cmap)

% plot
[a2,ha2] = contourf(xx,yy,zz,cvec2,'LineColor','none');

%put colorbar
cbarf(cvec2(1:end),cvec2(1+3:end-3),'vertical','nonlinear',1);


%tcmap=jet(length(cvec2)-4);
%tcmap=[squeeze(tcmap(1,:));squeeze(tcmap(1,:));squeeze(tcmap(1,:));squeeze(tcmap(1,:));squeeze(tcmap(1,:));tcmap;squeeze(tcmap(end,:));tcmap;squeeze(tcmap(end,:));tcmap;squeeze(tcmap(end,:));tcmap;squeeze(tcmap(end,:));squeeze(tcmap(end,:))];
%[a2,ha2] = contourf(input.x/1000,input.y/1000,((squeeze(input.KN2(n,:,:)-1)/1+1)).*squeeze(input.shelfELM(n,:,:))./squeeze(sign(input.vmELM(n,:,:))),cvec2,'LineColor','none');%,'LineStyle','none','LineColor','none');

% fix patch colors
  p=get(ha2,'children');
  thechild=get(p,'CData');   
  cdat=cell2mat(thechild);
  %loop through and manually set facecolor of each patch to the colormap you  made:
  for i=1:length(cvec2)
   set(p(cdat==cvec2(i)),'Facecolor',cmap(i,:));
  end
clear zz i p thechild cdat


hold on
%% LAYER 2
% background velocity field
zz = nthroot(input.vm_obs,7).*isnan(squeeze(input.shelfELM(n,:,:)));
cvec1 = nthroot([0 1 10 20 50 100 200 500 1000 2000 5000],7);
[a1,ha1]=contourf(xx,yy,zz,cvec1,'LineStyle','none');
caxis([0 nthroot(max(max(input.vm_obs)),7)])
zz = input.vm_obs./double(input.mainland);
contour(xx,yy,zz,[-5000 100 10000],'LineWidth',1.2,'LineColor','k','LineStyle','--')

if(n<=6 || n==11 || n==16 || n==21)
else
contour(xx,yy,zeta,[-10000 0.5 10000],'LineColor',[1 1 1])
end
contour(xx,yy,zeta2,[-10000 0.5 10000],'LineColor',[1 0 0])

axis equal
cmap=colormap(strcat('gray(',num2str(length(cvec1)),')'));
%cbarf(cvec1,cvec1(1:end-1),'vertical','nonlinear',7);
xlabel('zonal extent [km]')
ylabel('meridional extent [km] ')
title('buttressing [-] ')

print('-dpng','-r600',strcat('./charts/butter_n',num2str(n),'_eig2_dir.png'))
%print('-dpng','-r600',strcat('./charts/buttressing_n',num2str(n),'_eig2_dir.png'))
%close

end








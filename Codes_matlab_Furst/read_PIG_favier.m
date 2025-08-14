% Matlab Script to load PIG calving experiments
% from Favier et al. (2014)
%
% INPUT:        bedmap2_surface.nc
%               bedmap2_bed.nc
%               bedmap2_thickness.nc
%
% scrupt by Johannes J. Fuerst, 2014



%close all
%clear all

cc  = [0 0 0 0];
mm1 = [1 11 20 30 40];
cm1 = [0 0 1 1 0];
mm2 = [1 10 20 30 40];
cm2 = [0 0 1 1 1];
mm3 = [1 10 20 30 40];
cm3 = [0 0 0 1 1];
mm4 = [1 10 20 30 40];
cm4 = [0 1 1 1];



for n=1:nexp

% =====================================================================
% DATA SPECIFICATIONS
% =====================================================================
if(n==1)
twister    = 1;
dir_name   = strcat('/../../../input/favier2014/relaxation/');
file_name  = strcat('export_relax_yr68_v05.csv');
elseif(n>1&&n<=5)
twister    = cc(n-1);
dir_name   = strcat('/../../../input/favier2014/calving/c',num2str(n-1),'/');
file_name  = strcat('export_c',num2str(n-1),'_v05.csv');
elseif(n>5&&n<=10)
twister    = cm1(n-5);
dir_name   = strcat('/../../../input/favier2014/melting/m1/');
if(mm1(n-5)<10)
file_name  = strcat('export_m1_yr0',num2str(mm1(n-5)),'_v05.csv');
else
file_name  = strcat('export_m1_yr',num2str(mm1(n-5)),'_v05.csv');
end
elseif(n>10&&n<=15)
twister    = cm2(n-10);
dir_name   = strcat('/../../../input/favier2014/melting/m2/');
if(mm2(n-10)<10)
file_name  = strcat('export_m2_yr0',num2str(mm2(n-10)),'_v05.csv');
else
file_name  = strcat('export_m2_yr',num2str(mm2(n-10)),'_v05.csv');
end
elseif(n>15&&n<=20)
twister    = cm3(n-15);
dir_name   = strcat('/../../../input/favier2014/melting/m3/');
if(mm3(n-15)<10)
file_name  = strcat('export_m3_yr0',num2str(mm3(n-15)),'_v05.csv');
else
file_name  = strcat('export_m3_yr',num2str(mm3(n-15)),'_v05.csv');
end
elseif(n>20&&n<=24)
twister    = cm4(n-20);
dir_name   = strcat('/../../../input/favier2014/melting/m4/');
if(mm4(n-20)<10)
file_name  = strcat('export_m4_yr0',num2str(mm4(n-20)),'_v05.csv');
else
file_name  = strcat('export_m4_yr',num2str(mm4(n-20)),'_v05.csv');
end
elseif(n>=25)
if(n==26) %calving experiments after 50yr
twister = 1;
else
twister = 0;
end
dir_name   = strcat('/../../../input/favier2014/calving/c',num2str(n-24),'/');
file_name  = strcat('export_50yr_c',num2str(n-24),'_v05.csv');
end

AAA = load(strcat(dir_name,file_name));
c(n).x    = squeeze(AAA(:,15));
c(n).y    = squeeze(AAA(:,16));
c(n).z    = squeeze(AAA(:,17));
c(n).vx   = squeeze(AAA(:,12));
c(n).vy   = squeeze(AAA(:,13));
c(n).vz   = squeeze(AAA(:,14));
c(n).txx  = squeeze(AAA(:, 6))*MPa;
c(n).tyy  = squeeze(AAA(:, 7))*MPa;
c(n).tzz  = squeeze(AAA(:, 8))*MPa;
c(n).txy  = squeeze(AAA(:, 9))*MPa;
c(n).txz  = squeeze(AAA(:,10))*MPa;
c(n).tyz  = squeeze(AAA(:,11))*MPa;
if(twister==1)
c(n).bot  = squeeze(AAA(:, 4));
c(n).sur  = squeeze(AAA(:, 3));
c(n).groundedmask  = squeeze(AAA(:, 5));
else
c(n).bot  = squeeze(AAA(:, 5));
c(n).sur  = squeeze(AAA(:, 4));
c(n).groundedmask  = squeeze(AAA(:, 3));
end
c(n).thi  = squeeze(AAA(:, 2));
c(n).bed  = squeeze(AAA(:, 1));
display(strcat('3D field size:',num2str(size(c(n).x,1))))
%keep surface data only (1 out of 11 layers; each 25577 nodes)
c(n).groundedmask(c(n).bot~=c(n).z) = [];
x0 = c(n).x;
y0 = c(n).y;
x0(c(n).bot~=c(n).z) = [];
y0(c(n).bot~=c(n).z) = [];
c(n).bot(c(n).bot~=c(n).z) = [];
c(n).shelfmask             = c(n).groundedmask;
c(n).shelfmask(c(n).shelfmask == 1) = 0;

c(n).x(c(n).sur~=c(n).z)   = [];
c(n).y(c(n).sur~=c(n).z)   = [];
c(n).vx(c(n).sur~=c(n).z)  = [];
c(n).vy(c(n).sur~=c(n).z)  = [];
c(n).txx(c(n).sur~=c(n).z) = [];
c(n).tyy(c(n).sur~=c(n).z) = [];
c(n).tzz(c(n).sur~=c(n).z) = [];
c(n).txy(c(n).sur~=c(n).z) = [];
c(n).txz(c(n).sur~=c(n).z) = [];
c(n).tyz(c(n).sur~=c(n).z) = [];
c(n).thi(c(n).sur~=c(n).z) = [];
c(n).bed(c(n).sur~=c(n).z) = [];
c(n).sur(c(n).sur~=c(n).z) = [];
c(n).z                     = c(n).sur;
c(n).bott                  = zeros(size(c(n).x));
for kk=1:length(c(n).x)
    kindex = find(c(n).x(kk).^2+c(n).y(kk).^2==x0.^2+y0.^2,1);
    kkindex(kk) = kindex;
    c(n).bott(kk) = c(n).bot(kindex);
%    c(n).shelfmask(kk) = NaN;
end
display(strcat('2D field size:',num2str(size(c(n).x,1))))
clear AAA y0 x0 kk kindex kkindex
c(n).domain_mask           = ones(size(c(n).x));
c(n).shelfmask             = ones(size(c(n).x));
c(n).shelfmask(abs(c(n).bott-c(n).bed) < 1e-6) = NaN;
%}
end

clear n

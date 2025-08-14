% regrid PIG data on BEDMAP2
% INPUT:        see scripts in ..
%               meta_setup.m
%
% OUTPUT:       input.KN1
%               input.KN2
%               input.KN3
%               input.KT
%               input.vxELM
%               input.vyELM
%               input.domain_mask
%
% by Johannes J. Fuerst

% =====================================================================
% CONSTANTS
% =====================================================================

% =====================================================================
% VARIABLES
% =====================================================================

% =====================================================================
% PROCESSING
% =====================================================================

input.KN1      = zeros(nexp,size(input.x,1),size(input.x,2));
input.KN2      = zeros(nexp,size(input.x,1),size(input.x,2));
input.KN3      = zeros(nexp,size(input.x,1),size(input.x,2));
input.KT       = zeros(nexp,size(input.x,1),size(input.x,2));
input.vxELM    = zeros(nexp,size(input.x,1),size(input.x,2));
input.vyELM    = zeros(nexp,size(input.x,1),size(input.x,2));
input.thiELM   = zeros(nexp,size(input.x,1),size(input.x,2));
input.shelfELM = zeros(nexp,size(input.x,1),size(input.x,2))*NaN;
input.domainELM = zeros(nexp,size(input.x,1),size(input.x,2))*NaN;

%XX = reshape(input.x,1,size(input.x,1)*size(input.x,2));
%YY = reshape(input.y,1,size(input.x,1)*size(input.x,2));

for n = 1:nexp
    display(strcat('INTERPOLATE SHELF step :',num2str(n)))
for i = 1:size(input.x,1)
for j = 1:size(input.x,2)
    distit    = sqrt((c(n).x-input.x(i,j)).^2+(c(n).y-input.y(i,j)).^2);
    mindist = min(distit);
    kindex  = find(distit==mindist(1),1);
    if(mindist<3e3)%necessary in the middle of the coarse gridded shelf
    input.shelfELM(n,i,j)  = c(n).shelfmask(kindex);
    input.domainELM(n,i,j) = c(n).domain_mask(kindex);
    end
end
end
%{
    %for k = 1:length(c(n).x)
    %dist = sqrt((c(n).x(k)-XX).^2+(c(n).y(k)-YY).^2);
    %mindist = min(dist);
    %kindex = find(dist==mindist(1),1);
    %[i,j] = ind2sub(size(input.x),kindex);
    %input.shelfELM(n,i,j) = c(n).shelfmask(k);
    %end
%}
end

for n = 1:nexp

KN1       = griddata(c(n).x,c(n).y,c(n).KN1,input.x,input.y);
input.KN1(n,:,:) = KN1.*squeeze(input.shelfELM(n,:,:));

KN2       = griddata(c(n).x,c(n).y,c(n).KN2,input.x,input.y);
input.KN2(n,:,:) = KN2.*squeeze(input.shelfELM(n,:,:));

KN3       = griddata(c(n).x,c(n).y,c(n).KN3,input.x,input.y);
input.KN3(n,:,:) = KN3.*squeeze(input.shelfELM(n,:,:));

KT        = griddata(c(n).x,c(n).y,c(n).KT,input.x,input.y);
input.KT(n,:,:)  = KT.*squeeze(input.shelfELM(n,:,:));

vx        = griddata(c(n).x,c(n).y,c(n).vx,input.x,input.y);
input.vxELM(n,:,:) = vx.*squeeze(input.shelfELM(n,:,:));

vy        = griddata(c(n).x,c(n).y,c(n).vy,input.x,input.y);
input.vyELM(n,:,:) = vy.*squeeze(input.shelfELM(n,:,:));

thi       = griddata(c(n).x,c(n).y,c(n).thi,input.x,input.y);
input.thiELM(n,:,:) = thi.*squeeze(input.shelfELM(n,:,:));

clear KN1 KN2 KN3 KT vx vy DM distit mindist kindex thi

end
clear i j k n

input.vmELM = sqrt(input.vxELM.^2+input.vyELM.^2);


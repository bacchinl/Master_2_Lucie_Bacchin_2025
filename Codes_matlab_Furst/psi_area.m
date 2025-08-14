% calculate relative passive shelf area
% INPUT:        see scripts in ..
%               meta_setup.m
%
% OUTPUT:       NONE
%
% by Johannes J. Fuerst

clear shelfarea deadshelfarea deadshelfarea2
% TABLE OF BUTTRESSING
for kk = 1:nanmax(nanmax(input.shelfID))
%    display(strcat('Processing step: ',num2str(kk),' on :',shelf_names(kk),' ice shelf'));
    shelfID                    = input.shelfID;
    shelfID(input.shelfID~=kk) = NaN;
    shelfID(shelfID==kk)       = 1;
    shelfarea(kk)              = nansum(nansum(shelfID));
    butter                     = kn3kn3.*input.shelfmask;
    butter(input.shelfID~=kk)  = NaN;
    butter(butter<butthreshold) = NaN;
    butter(isnan(butter)~=1)   = 1;
    deadshelfarea(kk)          = shelfarea(kk)-nansum(nansum(butter));
    butter2                    = kn1kn1.*input.shelfmask;
    butter2(input.shelfID~=kk) = NaN;
    butter2(butter2<butthreshold) = NaN;
    butter2(isnan(butter2)~=1) = 1;
    deadshelfarea2(kk)         = shelfarea(kk)-nansum(nansum(butter2));
    divider                    = (kn3kn3-kn1kn1).*input.shelfmask;
    divider(input.shelfID~=kk) = NaN;
    division(kk)               = nanmean(nanmean(divider));
%    display(strcat('Shelf grid cells :',num2str(shelfarea(kk)),'    Dead shelf-ice cells :',num2str(deadshelfarea(kk)),'   ratio:',num2str(floor(100*deadshelfarea(kk)/shelfarea(kk))),'.',num2str(round((100*deadshelfarea(kk)/shelfarea(kk)-floor(100*deadshelfarea(kk)/shelfarea(kk)))*10)),' %'));
    display(strcat(shelf_names(kk),' ice shelf : ','    ----   ratio :',num2str(floor(100*deadshelfarea(kk)/shelfarea(kk))),'.',num2str(round((100*deadshelfarea(kk)/shelfarea(kk)-floor(100*deadshelfarea(kk)/shelfarea(kk)))*10)),' %','    ----   ratio (flow):',num2str(floor(100*deadshelfarea2(kk)/shelfarea(kk))),'.',num2str(round((100*deadshelfarea2(kk)/shelfarea(kk)-floor(100*deadshelfarea2(kk)/shelfarea(kk)))*10)),' %'));
 
    clear butter shelfID
end

shelfarea(end+1)      = sum(shelfarea);
deadshelfarea(end+1)  = sum(deadshelfarea);
deadshelfarea2(end+1) = sum(deadshelfarea2);
display('ANTARCTICA')
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells :',num2str(deadshelfarea(end)),'   ratio:',num2str(floor(100*deadshelfarea(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea(end)/shelfarea(end)-floor(100*deadshelfarea(end)/shelfarea(end)))*10)),' %'));
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells (flow) :',num2str(deadshelfarea2(end)),'   ratio (flow):',num2str(floor(100*deadshelfarea2(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea2(end)/shelfarea(end)-floor(100*deadshelfarea2(end)/shelfarea(end)))*10)),' %'));
shelfarea(end)         = [];
deadshelfarea(end)     = [];
deadshelfarea2(end)    = [];

shelfarea(end+1) = sum(shelfarea([2:8]));
deadshelfarea(end+1) = sum(deadshelfarea([2:8]));
deadshelfarea2(end+1) = sum(deadshelfarea2([2:8]));
display('WEDDELL SEA')
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells :',num2str(deadshelfarea(end)),'   ratio:',num2str(floor(100*deadshelfarea(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea(end)/shelfarea(end)-floor(100*deadshelfarea(end)/shelfarea(end)))*10)),' %'));
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells (flow) :',num2str(deadshelfarea2(end)),'   ratio (flow):',num2str(floor(100*deadshelfarea2(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea2(end)/shelfarea(end)-floor(100*deadshelfarea2(end)/shelfarea(end)))*10)),' %'));
shelfarea(end)         = [];
deadshelfarea(end)     = [];
deadshelfarea2(end)    = [];

shelfarea(end+1) = sum(shelfarea([9:12]));
deadshelfarea(end+1) = sum(deadshelfarea([9:12]));
deadshelfarea2(end+1) = sum(deadshelfarea2([9:12]));
display('WEST INDIAN OCEAN')
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells :',num2str(deadshelfarea(end)),'   ratio:',num2str(floor(100*deadshelfarea(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea(end)/shelfarea(end)-floor(100*deadshelfarea(end)/shelfarea(end)))*10)),' %'));
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells (flow) :',num2str(deadshelfarea2(end)),'   ratio (flow):',num2str(floor(100*deadshelfarea2(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea2(end)/shelfarea(end)-floor(100*deadshelfarea2(end)/shelfarea(end)))*10)),' %'));
shelfarea(end)         = [];
deadshelfarea(end)     = [];
deadshelfarea2(end)    = [];

shelfarea(end+1) = sum(shelfarea([13:22]));
deadshelfarea(end+1) = sum(deadshelfarea([13:22]));
deadshelfarea2(end+1) = sum(deadshelfarea2([13:22]));
display('EAST INDIAN OCEAN')
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells :',num2str(deadshelfarea(end)),'   ratio:',num2str(floor(100*deadshelfarea(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea(end)/shelfarea(end)-floor(100*deadshelfarea(end)/shelfarea(end)))*10)),' %'));
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells (flow) :',num2str(deadshelfarea2(end)),'   ratio (flow):',num2str(floor(100*deadshelfarea2(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea2(end)/shelfarea(end)-floor(100*deadshelfarea2(end)/shelfarea(end)))*10)),' %'));
shelfarea(end)         = [];
deadshelfarea(end)     = [];
deadshelfarea2(end)    = [];

shelfarea(end+1) = sum(shelfarea([23:25]));
deadshelfarea(end+1) = sum(deadshelfarea([23:25]));
deadshelfarea2(end+1) = sum(deadshelfarea2([23:25]));
display('ROSS SEA')
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells :',num2str(deadshelfarea(end)),'   ratio:',num2str(floor(100*deadshelfarea(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea(end)/shelfarea(end)-floor(100*deadshelfarea(end)/shelfarea(end)))*10)),' %'));
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells (flow) :',num2str(deadshelfarea2(end)),'   ratio (flow):',num2str(floor(100*deadshelfarea2(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea2(end)/shelfarea(end)-floor(100*deadshelfarea2(end)/shelfarea(end)))*10)),' %'));
shelfarea(end)         = [];
deadshelfarea(end)     = [];
deadshelfarea2(end)    = [];

shelfarea(end+1) = sum(shelfarea([26:32]));
deadshelfarea(end+1) = sum(deadshelfarea([26:32]));
deadshelfarea2(end+1) = sum(deadshelfarea2([26:32]));
display('AMUNDSEN SEA')
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells :',num2str(deadshelfarea(end)),'   ratio:',num2str(floor(100*deadshelfarea(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea(end)/shelfarea(end)-floor(100*deadshelfarea(end)/shelfarea(end)))*10)),' %'));
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells (flow) :',num2str(deadshelfarea2(end)),'   ratio (flow):',num2str(floor(100*deadshelfarea2(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea2(end)/shelfarea(end)-floor(100*deadshelfarea2(end)/shelfarea(end)))*10)),' %'));
shelfarea(end)         = [];
deadshelfarea(end)     = [];
deadshelfarea2(end)    = [];

shelfarea(end+1) = sum(shelfarea([33:39]));
deadshelfarea(end+1) = sum(deadshelfarea([33:39]));
deadshelfarea2(end+1) = sum(deadshelfarea2([33:39]));
display('BELLINGHAUSEN SEA')
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells :',num2str(deadshelfarea(end)),'   ratio:',num2str(floor(100*deadshelfarea(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea(end)/shelfarea(end)-floor(100*deadshelfarea(end)/shelfarea(end)))*10)),' %'));
display(strcat('Shelf grid cells :',num2str(shelfarea(end)),'    Dead shelf-ice cells (flow) :',num2str(deadshelfarea2(end)),'   ratio (flow):',num2str(floor(100*deadshelfarea2(end)/shelfarea(end))),'.',num2str(round((100*deadshelfarea2(end)/shelfarea(end)-floor(100*deadshelfarea2(end)/shelfarea(end)))*10)),' %'));
shelfarea(end)         = [];
deadshelfarea(end)     = [];
deadshelfarea2(end)    = [];

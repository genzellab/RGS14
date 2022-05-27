%% Process data for ripple analysis (in terms of shape)

data=rat3_sd14_c1_bp_rgs; % manually enter the data to be changed

x=cell(1,9);
for i=1:9
    x{1,i}=cell(1,3);
end

for i=1:9
    x{1,i}{1,1} = [];
    x{1,i}{1,3} = [];
    x{1,i}{1,2} = [];
end

for i=1:size(data,1)
    if data(i,9)=="ps"
        a=data(i,2);
        x{1,1}{1,1} = [x{1,1}{1,1},a];
        x{1,1}{1,3} = [x{1,1}{1,3},data(i,3)];
        x{1,1}{1,2} = [x{1,1}{1,2},data(i,4)];
    elseif data(i,9)=="pt1"
        a=data(i,2);
        x{1,2}{1,1} = [x{1,2}{1,1},a];
        x{1,2}{1,3} = [x{1,2}{1,3},data(i,3)];
        x{1,2}{1,2} = [x{1,2}{1,2},data(i,4)];
    elseif data(i,9)=="pt2"
        a=data(i,2);
        x{1,3}{1,1} = [x{1,3}{1,1},a];
        x{1,3}{1,3} = [x{1,3}{1,3},data(i,3)];
        x{1,3}{1,2} = [x{1,3}{1,2},data(i,4)];
    elseif data(i,9)=="pt3"
        a=data(i,2);
        x{1,4}{1,1} = [x{1,4}{1,1},a];
        x{1,4}{1,3} = [x{1,4}{1,3},data(i,3)];
        x{1,4}{1,2} = [x{1,4}{1,2},data(i,4)];
    elseif data(i,9)=="pt4"
        a=data(i,2);
        x{1,5}{1,1} = [x{1,5}{1,1},a];
        x{1,5}{1,3} = [x{1,5}{1,3},data(i,3)];
        x{1,5}{1,2} = [x{1,5}{1,2},data(i,4)];
    elseif data(i,9)=="pt5.1"
        a=data(i,2);
        x{1,6}{1,1} = [x{1,6}{1,1},a];
        x{1,6}{1,3} = [x{1,6}{1,3},data(i,3)];
        x{1,6}{1,2} = [x{1,6}{1,2},data(i,4)];
    elseif data(i,9)=="pt5.2"
        a=data(i,2);
        x{1,7}{1,1} = [x{1,7}{1,1},a];
        x{1,7}{1,3} = [x{1,7}{1,3},data(i,3)];
        x{1,7}{1,2} = [x{1,7}{1,2},data(i,4)];
    elseif data(i,9)=="pt5.3"
        a=data(i,2);
        x{1,8}{1,1} = [x{1,8}{1,1},a];
        x{1,8}{1,3} = [x{1,8}{1,3},data(i,3)];
        x{1,8}{1,2} = [x{1,8}{1,2},data(i,4)];
    elseif data(i,9)=="pt5.4"
        a=data(i,2);
        x{1,9}{1,1} = [x{1,9}{1,1},a];
        x{1,9}{1,3} = [x{1,9}{1,3},data(i,3)];
        x{1,9}{1,2} = [x{1,9}{1,2},data(i,4)];
    end
end

% Do not forgt to change names here
rat3_sd14_c1_bp_rgs=x;
save rat3_sd14_c1_bp_rgs.mat rat3_sd14_c1_bp_rgs


%save rat9_veh.mat












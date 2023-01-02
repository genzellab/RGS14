clear variables
addpath(genpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/CorticoHippocampal'))
addpath ('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/ADRITOOLS')
cd('/Volumes/Samsung_T5/Milan_DA/GCsd')
rat_folder = getfolder;
k = 8; 
cd(rat_folder{k});
path = cd;
dinfo = dir(path);
dinfo = {dinfo.name};
dinfo(1) = [];
dinfo(1) = [];
idx = find(cellfun(@(x)~isempty(strfind(x,'compilation')),dinfo));
dinfo(idx) = [];
for i = 1: length(dinfo)
    
    filename = dinfo{i};
    filename = strsplit(filename, '_');
    
    idx = find(contains(filename, 'rat', 'IgnoreCase',true));
        if length(idx)>1
         ratn = str2double(filename{idx(end)}(4));
        else
         ratn = str2double(filename{idx}(4));
        end
    idx = find(contains(filename, 'sd', 'IgnoreCase',true));
        SDn = str2double(filename{idx}(1,3:end));
    
 OScon = [{'HC'}, {'CON'},{'OD'},{'OR'}];
 
   for c  = 1:length(OScon)
        idx = find(contains(filename, OScon{c}, 'IgnoreCase',true));
        if ~isempty(idx)
            condition = filename{idx};
            break
        else 
            continue 
        end 
   end 
trial = [{'ps'},{'pt1'},{'pt2'},{'pt3'},{'pt4'},{'pt5.1'},{'pt5.2'},{'pt5.3'},{'pt5.4'}];

GC = load(dinfo{i});
fname  = fieldnames(GC);
GC = GC.(fname{1});

for j = 1:length(GC)
    GC{j}(cellfun('isempty', GC{j}(:,1)),:) = [];
    l = size(GC{j},1);
    n = string(1:l');
    uid = cell([l,1]);
    uid(:,1) =  {'UID'};
    for d = 1:l
        uid{d,1} = strcat(uid{d,1},n(d));
    end
    GC{j}(:,5) = dinfo(i);
    GC{j}(:,6) = {ratn};
    GC{j}(:,7) = {SDn};
    GC{j}(:,8) = {condition};
    GC{j}(:,9) = trial(j); 
    GC{j}(:,10) = uid;
end
    if  contains(fname{1},'broadband')
    GC_window_broadband_ripples_total = GC;
    save (dinfo{i} ,'GC_window_broadband_ripples_total')
    else
    GC_window_ripples_total = GC;
    save (dinfo{i} ,'GC_window_ripples_total')
    end
end







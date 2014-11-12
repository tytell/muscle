function plot_duty_data

filenames = getfilenames('dutydata*.mat');
fields = {'x','lc','vc','Pc','fexp'};

phitest2 = 0:0.2:0.8;
dutyvals = [0.1 0.36 0.5 0.55];
k3test = 0.6:0.1:1.4;
k4test = [0.8 1 1.2];
[phivals2,dutyvals2, k3vals2, k4vals2] = ndgrid(phitest2,dutyvals,k3test,k4test);
assert(numel(phivals2) == length(filenames));

%start with "normal" k3, k4 values
k3ind = find(k3test == 1);
k4ind = find(k4test == 1);

[i1,i2,i3,i4] = ndgrid(1:length(phitest2),1:length(dutyvals),1:length(k3test),1:length(k4test));
ind = sub2ind(size(phivals2), i1,i2,i3,i4);

timedWaitBar(0,'Loading data');
showind = ind(:,:,k3ind,k4ind);
for i = 1:numel(showind)
    F = load(filenames{showind(i)});
    
    for j = 1:length(fields)
        dutydata.(fields{j}) = F.data1.(fields{j});
    end
    timedWaitBar(i/numel(showind),'Loading data');
end
timedWaitBar(1);

putvar dutydata;



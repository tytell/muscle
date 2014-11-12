function plot_muscle_duty_data

filename = 'test_muscle_duty.h5';
loadvars = {'Ca','Caf','m','lc','vc','Pc','L','V',...
    'fexp','fx','phi','t','k3','k4','duty'};

par.k30 = 51.3537;              % 1/s
par.k40 = 19.3801;              % 1/s

phi1 = 0:0.2:0.8;
duty1 = [0.1 0.2 0.3 0.36 0.4 0.5];
k31 = [0.8 1 1.2] * par.k30;
k41 = [0.8 1 1.2] * par.k40;

%order of dimensions
% phi, duty, k3, k4

for i = 1:length(loadvars)
    data.(loadvars{i}) = h5read(filename,['/' loadvars{i}]);
end

t12 = log(0.5) ./ real(squeeze(data.fexp(:,:,:,:,:)));

clf;
plot(data.t, data.Pc(:,:,1,1,1));
plotgroups(flatten(data.duty(:,:,2,:)), flatten(t12(2,:,:,2,:)),{flatten(data.k4(:,:,2,:))},{'cm'},'means','error','std','xoffset',0.01);


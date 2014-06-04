function merge_duty_data(filenames,outfile)

phitest2 = 0:0.2:0.8;
dutyvals = [0.1 0.36 0.5 0.55];
k3test = 0.6:0.1:1.4;
k4test = [0.8 1 1.2];
sz = [length(phitest2) length(dutyvals) length(k3test) length(k4test)];
assert(prod(sz) == length(filenames));

F = load(filenames{1});

if (exist(outfile,'file') && ...
        inputyn(sprintf('Overwrite output file %s',outfile),'default',false))
    delete(outfile);
end

fields = fieldnames(F.data1);
sizes = cell(size(fields));
for k = 1:length(fields)
    if (isnumeric(F.data1.(fields{k})))
        if (numel(F.data1.(fields{k})) == 1)
            sz1 = [];
        else
            sz1 = size(F.data1.(fields{k}));
            if (sz1(2) == 1)
                sz1 = sz1(1);
            end
        end
        sizes{k} = sz1;
        h5create(outfile,['/' fields{k}], [sz1 sz]);
    end
end

timedWaitBar(0,'Merging data...');
for i = 1:length(filenames)
    F = load(filenames{i},'data1');
    data1 = struct2cell(F.data1);

    [a,b,c,d] = ind2sub(sz,i);
    
    for k = 1:length(data1)
        if (isnumeric(data1{k}))
            start = [ones(1,length(sizes{k})) a b c d];
            if (numel(data1{k}) == 1)
                count = ones(1,length(sz));
            else
                sz1 = size(data1{k});
                if (sz1(2) == 1)
                    count = [sz1(1) ones(1,length(sz))];
                else
                    count = [sz1 ones(1,length(sz))];
                end
            end
            h5write(outfile,['/' fields{k}], data1{k}, start,count);
        end
    end
    timedWaitBar(i/length(filenames),'Merging data...');
end
timedWaitBar(1);


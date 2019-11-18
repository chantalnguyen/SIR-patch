times = 0:0.5:200;
data = zeros(length(times),Npatch,length(XX));
for i = 1:length(XX)
    weights = [0.25, 0.50, 0.75, 1];
    xxx = XX{i}(:,2:end,:);
    xxxx = zeros(size(xxx,1),Npatch);
    for k = 1:size(xxx,3)
        xxxx(:,k) = (weights * xxx(:,:,k)')';
    end
%     xxx = squeeze(sum(XX{i}(:,2:end,:),2))./9;
    
    for j = 1:Npatch
        data(:,j,i) = interp1(t{i}(1:end-1),xxxx(1:end-1,j),times);
        Ind = find(isnan(data(:,j,i)),1,'first');
        if ~isempty(Ind)
            data(Ind:end,j,i) = data(Ind-1,j,i);
            
        end
    end
end

dataa = mean(data,3);
%%
figure()
for i = 1:size(dataa,1)
    bcolor(reshape(dataa(i,:),5,5));
    colorbar;
    caxis([0 9]);
    title(['t = ' num2str(times(i))])
    pause(0.01)
end
%%
times = 0:1:250;
XXX = cell(size(XX));
for i = 1:length(XX)
    XXX{i} = squeeze(sum(XX{i},3));
end
data = zeros(length(times),Npatch,length(XXX));
for i = 1:length(XXX)
%     weights = [0.25, 0.50, 0.75, 1];
    weights = [1 1 1 1];
    xxx = XXX{i}(:,2:end,:);
    xxxx = zeros(size(xxx,1),Npatch);
    for k = 1:size(xxx,3)
        xxxx(:,k) = (weights * xxx(:,:,k)')';
    end
    
    for j = 1:Npatch
        data(:,j,i) = interp1(t{i}(1:end-1),xxxx(1:end-1,j),times);
        Ind = find(isnan(data(:,j,i)),1,'first');
        if ~isempty(Ind)
            data(Ind:end,j,i) = data(Ind-1,j,i);
            
        end
    end
end

dataa = mean(data,3);
%%
figure()
for i = 1:size(dataa,1)
    bcolor(reshape(dataa(i,:),5,5));
    colorbar;
    caxis([0 9]);
    title(['t = ' num2str(times(i))])
    pause(0.01)
end

%% number of infected per patch over time
times = 0:1:250;
YYY = cell(size(YY));
for i = 1:length(YY)
    YYY{i} = squeeze(sum(YY{i},3));
end
data = zeros(length(times),Npatch,length(YYY));
for i = 1:length(YYY)
%     weights = [0.25, 0.50, 0.75, 1];
    weights = [1];
    yyy = YYY{i}(:,2:end,:);
    yyyy = zeros(size(yyy,1),Npatch);
    for k = 1:size(yyy,3)
        yyyy(:,k) = (weights * yyy(:,:,k)')';
    end
    
    for j = 1:Npatch
        data(:,j,i) = interp1(t{i}(1:end-1),yyyy(1:end-1,j),times);
        Ind = find(isnan(data(:,j,i)),1,'first');
        if ~isempty(Ind)
            data(Ind:end,j,i) = data(Ind-1,j,i);
            
        end
    end
end

dataa = mean(data,3);
%%
figure()
for i = 1:size(dataa,1)
    bcolor(reshape(dataa(i,:),5,5));
    colorbar;
    caxis([0 3]);
    title(['t = ' num2str(times(i))])
    pause(0.01)
end
%%
figure()
for i = 1:size(dataa,2)
    plot(dataa(:,i)); hold on;
end
%% plot final epidemic size
figure()
for i = 1:size(yy)
    plot(t{i},225-yy{i}(:,1)); hold on;
end
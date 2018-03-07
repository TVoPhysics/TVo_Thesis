%load in Data_in (time series { t,DCPD[?],PZT[Volts] })
data_in = load('InputData.txt');

[pks,locs] = findpeaks(data_in(:,2),'MinPeakDistance',500);

%load in model
model_in = load('InModel.txt');

semilogy(data_in(:,3),data_in(:,2),'--',...
        data_in(locs,3),pks,'o');
grid on
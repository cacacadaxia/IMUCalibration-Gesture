clear all;
% load('caldata.mat')
cal_data = dlmread('text.txt');
dt = 1/200;
time = dt:dt:(length(cal_data))*dt;time = time';
data = [time, cal_data];


m=size(data,1);
for i=1:m
%      data(i,2:4)=(Ta*Ka*(data(i,2:4)'+Ba))';
%      data(i,5:7)=(Tg*Kg*(data(i,5:7)'+Bg))';
     
     
   
end

acc = data(:,2:4);
gyro = data(:,5:7);

figure;plot(gyro);
figure;plot(acc);


% =========================================================================
%
%                 IMU的标定
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 5月12日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.标定陀螺仪系统误差
%        2.为什么选择旋转的量呢？这一点让我很费解
%        3.
%        4.
%
%--------------------------------------------------------------------------

clear all;
close all;

cal_data = dlmread('text.txt');

cal_data = cal_data(1800:end,:);
dt = 1/200;             %% 这是采样时间间隔
time = dt:dt:(length(cal_data))*dt;time = time';
data = [time, cal_data];

[~,fix_point,rotation] = FindFixDataxiaomi(data,20);%%fix_point是什么？

[Ta,Ka,Ba]=ICRA2014_acc(fix_point);

Bg = -mean(fix_point(:,4:6),1)';%%allan方差的方法吗？

n = size(rotation,1);

rotation{n+1}=Ta;
rotation{n+2}=Ka;
rotation{n+3}=Ba;
% rotation{n+4}=Bg;
rotation{n+4} = [0,0,0];

[ Tg , Kg ] = ICRA_2014_gyro(rotation);

%[Tm2a,Bm,Vm]=mag2acc_matrix(fix_point,Ta,Ka,Ba);

% [Tm2a,Bm,Vm,mag_strength] = Cal_mag4acc_frame(rotation,fix_point,Tg,Kg);

% Set_Bias_Gyro = [0.1,-0.4,1.5];

% See_Gesture( data,Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm,Set_Bias_Gyro);


% [Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm,mag_strength ] = ImuCalibration_Gesture(cal_data);
% 算法mag2acc_matrix假设重力与磁向量的夹角不变，算法Cal_mag4acc_frame利用不同姿态下传感器感受的磁通向量的变化与姿态变化的相关性，计算参数。

tpdata = data; %%保存之间的数据
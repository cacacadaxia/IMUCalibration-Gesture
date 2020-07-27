
% =========================================================================
%
%                 IMU�ı궨
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 5��12��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�궨������ϵͳ���
%        2.Ϊʲôѡ����ת�����أ���һ�����Һܷѽ�
%        3.
%        4.
%
%--------------------------------------------------------------------------

clear all;
close all;

cal_data = dlmread('text.txt');

cal_data = cal_data(1800:end,:);
dt = 1/200;             %% ���ǲ���ʱ����
time = dt:dt:(length(cal_data))*dt;time = time';
data = [time, cal_data];

[~,fix_point,rotation] = FindFixDataxiaomi(data,20);%%fix_point��ʲô��

[Ta,Ka,Ba]=ICRA2014_acc(fix_point);

Bg = -mean(fix_point(:,4:6),1)';%%allan����ķ�����

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
% �㷨mag2acc_matrix����������������ļнǲ��䣬�㷨Cal_mag4acc_frame���ò�ͬ��̬�´��������ܵĴ�ͨ�����ı仯����̬�仯������ԣ����������

tpdata = data; %%����֮�������
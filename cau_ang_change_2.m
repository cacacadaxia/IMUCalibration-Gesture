
% =========================================================================
%
%                 利用标定值计算角度
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 5月12日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.标定陀螺仪系统误差
%        2.
%        3.
%        4.
%
%--------------------------------------------------------------------------
cal_data = dlmread('Stationary_1.txt');
dt = 1/200;             %% 这是采样时间间隔
time = dt:dt:(length(cal_data))*dt;time = time';
data = [time, cal_data];
% data = tpdata;
m = size(data,1);

%====================================================
t(1)=0;
Pk_EKF=[];
eInt=[];
Pk_EKF_bias=[];
Bias_Ekf=zeros(1,9);
Pk_ESKF=[];
delta_X=zeros(1,6);
Bg_ESKF=zeros(1,3);
clear Q;
for i=1:m
    cabil = 0;
    if cabil == 1
        data(i,2:4) = (Ta*Ka*(data(i,2:4)' + Ba))';
        data(i,5:7) = (Tg*Kg*(data(i,5:7)' + Bg))';
    else
        data(i,2:4) = (data(i,2:4) + Ba') * 9.8;
        data(i,5:7) = (data(i,5:7) + Bg')* 0.0175;
    end
     norm_a(i)=norm(data(i,2:4));
     norm_g(i)=norm(data(i,5:7));
     
     if i>1
         t(i)=data(i,1)-data(i-1,1);
     end
     %%静止检测
    if norm_g(i)<0.0873    %5*pi/180
        q(i,:)=[1,0,0,0];
    else
        q(i,:) = axisAngle2quatern(data(i,5:7)/norm_g(i), norm_g(i)*t(i));%%真实角度
    end
    if i==1 
%         size(accMeg2qRichard(data(i,:)));
%         Q(i,:)  = accMeg2qRichard(data(i,:));  % only gyro

        Q(i,:) = cau_init(data');
        Q_RK4(i,:) = Q(i,:);              % only gyro in using RK4 updata
        %%进行初始化，刚开始的时候是随机的
        Q_random= randn(1,4);
        Q_random= Q_random/norm(Q_random);
        QfuseHL(i,:)=Q_random;            % high low pass filter
        QfuseEKF(i,:)=Q_random;           % EKF filter
        QfuseMahony(i,:)=Q_random;        % Mahony filter 
        QfuseEKF_bias(i,:)=Q_random;      % EKF filter to gyro bias
        QfuseESKF(i,:)=Q_random;          % ESKF filter
       %%问题是这样加噪后，bias
    else
        Q(i,:) = quaternProd(Q(i-1,:),q(i,:));    %Q(i-1,:)*q(i,:)
        
        %%直接积分计算角度
        Q_RK4(i,:) = attitude_update_RK4(Q_RK4(i-1,:)',t(i),data(i-1,5:7)',data(i,5:7)')';
        
        angsave(i,:) = qtpoula(Q_RK4(i,:));
        
        %% 下面3个都需要地磁计
        %%高通滤波器
%         QfuseHL(i,:) = HighLowPassFilter(QfuseHL(i-1,:),data(i,:),t(i));
        
        %%EKF滤波器
%         [QfuseEKF(i,:),Pk_EKF] = EkfFilter(QfuseEKF(i-1,:),data(i,:),t(i),Vm,Pk_EKF);
        
        %%直接滤波
%         [QfuseMahony(i,:),eInt]=MahonyFilter(QfuseMahony(i-1,:),data(i,:),t(i),Vm,eInt);  
        
        
        %% 利用卡尔曼滤波进行bias的估计过程
        %% 都需要地磁计
        %% 估计陀螺仪的固定偏差，比较准
        Set_Bias_Gyro = [0.1,-0.4,1.5]; %%虽然之前的gyro的bias去除掉了，但是这里又加偏差，也相对合理
        data_bias = [0,0,0,0,Set_Bias_Gyro,0,0,0];%%为什么要对陀螺仪的bias进行添加呢？还是固定不变？
        
%         [QfuseEKF_bias(i,:),Bias_Ekf(i,:),Pk_EKF_bias]  = EKF_Gyro_bias(QfuseEKF_bias(i-1,:),Bias_Ekf(i-1,:),data(i,:)+data_bias,t(i),Vm,Pk_EKF_bias);
       
%         [QfuseESKF(i,:),Bg_ESKF(i,:),delta_X(i,:),Pk_ESKF]  = ESKF(QfuseESKF(i-1,:),Bg_ESKF(i-1,:),delta_X(i-1,:),data(i,:)+data_bias,t(i),Vm,Pk_ESKF);
        
    end
end

figure;plot(data(:,1),Q_RK4(:,1:3));
figure;plot(data(:,1),angsave(:,1:3)/pi*180);legend 1 2 3

% figure;plot(1:m,QfuseEKF(:,1),'g', 1:m,QfuseMahony(:,1) ,1:m,QfuseHL(:,1) , 1:m,Q_RK4(:,1));legend 1 2 3 4

PLOT = 0;
if PLOT ==1
%%
figure(1)
p1(1)=subplot(6,1,1);
plot(1:m,QfuseHL(:,1),'r',1:m,QfuseEKF(:,1),'g', 1:m,QfuseMahony(:,1),'b',1:m,Q(:,1),'m',1:m,Q_RK4(:,1),'c');
ylim([-1,1]);
legend('HLP','EKF','Mahony','only gyro','gyro in RK4');

p1(2)=subplot(6,1,2);
plot(1:m,QfuseHL(:,2),'r',1:m,QfuseEKF(:,2),'g', 1:m,QfuseMahony(:,2),'b',1:m,Q(:,2),'m',1:m,Q_RK4(:,2),'c');
ylim([-1,1]);

p1(3)=subplot(6,1,3);
plot(1:m,QfuseHL(:,3),'r',1:m,QfuseEKF(:,3),'g', 1:m,QfuseMahony(:,3),'b',1:m,Q(:,3),'m',1:m,Q_RK4(:,3),'c');
ylim([-1,1]);

p1(4)=subplot(6,1,4);
plot(1:m,QfuseHL(:,4),'r',1:m,QfuseEKF(:,4),'g', 1:m,QfuseMahony(:,4),'b',1:m,Q(:,4),'m',1:m,Q_RK4(:,4),'c');
ylim([-1,1]);

p1(5)=subplot(6,1,5);
plot(1:m,norm_g);

p1(6)=subplot(6,1,6);
plot(1:m,norm_a);
linkaxes(p1,'x');

end

%%
function init_quat = cau_init(imuacc)
init_a = mean(imuacc(1:3,1:200),2);
init_a = init_a / norm(init_a);

init_psi =  0;
init_theta = -asin(init_a(1));
init_phi = atan2(init_a(2), init_a(3));

init_quat = angle2quat(init_psi, init_theta, init_phi);

end
function q = axisAngle2quatern(axis, angle)
    q0 = cos(angle./2);
    q1 = axis(:,1)*sin(angle./2);
    q2 = axis(:,2)*sin(angle./2);
    q3 = axis(:,3)*sin(angle./2); 
    q = [q0 q1 q2 q3];
end

function ab = quaternProd(a, b)
    ab(1) = a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
    ab(2) = a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
    ab(3) = a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
    ab(4) = a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1);
    if ab(1)<0
        ab=-ab;
    end
end


function R = quatern2rotMat(q)

    R(1,1) = q(1)^2+q(2)^2-q(3)^2-q(4)^2;
    R(1,2) = 2*(q(2)*q(3)+q(1)*q(4));
    R(1,3) = 2*(q(2)*q(4)-q(1)*q(3));
    R(2,1) = 2*(q(2)*q(3)-q(1)*q(4));
    R(2,2) = q(1)^2+q(3)^2-q(2)^2-q(4)^2;
    R(2,3) = 2*(q(3)*q(4)+q(1)*q(2));
    R(3,1) = 2*(q(2)*q(4)+q(1)*q(3));
    R(3,2) = 2*(q(3)*q(4)-q(1)*q(2));
    R(3,3) = q(1)^2+q(4)^2-q(2)^2-q(3)^2;
end



function q = accMeg2qRichard(data)


vX=cross(data(1,8:10),data(1,2:4));
vX=vX/norm(vX);
vY=cross(data(1,2:4),vX);
vY=vY/norm(vY);

qX = qUtoV(vX,[1,0,0]);

y= qMultiVec(vY, qX);
qY = qUtoV(y,[0,1,0]);

qx=[-qX(1),qX(2:4)];
qy=[-qY(1),qY(2:4)];

q =qMultiQ(qx,qy);
q=[q(1),-q(2:4)];
if q(1)<0
    q=-q;
end
end


function [qq]=qMultiQ(p,q)   %p*q
qq=[...
        p(1) * q(1) - p(2) * q(2) - p(3) * q(3) - p(4) * q(4)...
       ,p(2) * q(1) + p(1) * q(2) - p(4) * q(3) + p(3) * q(4)...
       ,p(3) * q(1) + p(4) * q(2) + p(1) * q(3) - p(2) * q(4)...
       ,p(4) * q(1) - p(3) * q(2) + p(2) * q(3) + p(1) * q(4)  ];

end

function q = qUtoV(u, v)        %two vetor rotation to quaternions
nu = u/norm(u);
nv = v/norm(v);

if (u*v' == -1)
    q = [0, [1,0,0]];
else
    half = (nu + nv)/norm(nu + nv);
    q = [nu*half',cross(nu, half)];
end
end

function [vector]=qMultiVec(vec,q)  %sensor frame to world frame
x = q(2);
y = q(3);
z = q(4);
w = q(1);

vecx = vec(1);
vecy = vec(2);
vecz = vec(3);

x_ =  w * vecx  +  y * vecz  -  z * vecy;
y_ =  w * vecy  +  z * vecx  -  x * vecz;
z_ =  w * vecz  +  x * vecy  -  y * vecx;
w_ = -x * vecx  -  y * vecy  -  z * vecz;

vector = [x_ * w  +  w_ * -x  +  y_ * -z  -  z_ * -y ...
    , y_ * w  +  w_ * -y  +  z_ * -x  -  x_ * -z ...
    , z_ * w  +  w_ * -z  +  x_ * -y  -  y_ * -x ...
    ];

end



function [R]=accMag2rotMat(data)

VerticalX=cross(data(1,8:10),data(1,2:4));
VerticalX=VerticalX/norm(VerticalX);
VerticalY=cross(data(1,2:4),VerticalX);
VerticalY=VerticalY/norm(VerticalY);
VerticalZ=data(1,2:4)/norm(data(1,2:4));
R=[VerticalX',VerticalY',VerticalZ'];   

end

function [accW]=accWorldframe(R,data)

accW=(R'*data(1,2:4)'-[0;0;9.8])/9.8;


end


function [Qk_plus1]=attitude_update_RK4(Qk,dt,gyro0,gyro1)
% RK4
% conference: A Robust and Easy to implement method for imu
% calibration without External Equipments

q_1=Qk;
k1=(1/2)*omegaMatrix(gyro0)*q_1;
q_2=Qk+dt*(1/2)*k1;
k2=(1/2)*omegaMatrix((1/2)*(gyro0+gyro1))*q_2;
q_3=Qk+dt*(1/2)*k2;
k3=(1/2)*omegaMatrix((1/2)*(gyro0+gyro1))*q_3;
q_4=Qk+dt*k3;
k4=(1/2)*omegaMatrix(gyro1)*q_4;
Qk_plus1=Qk+dt*(k1/6+k2/3+k3/3+k4/6);
Qk_plus1=Qk_plus1/norm(Qk_plus1);

if Qk_plus1(1)<0
    Qk_plus1=-Qk_plus1;
end

end

function [omega]=omegaMatrix(data)

% wx=data(1)*pi/180;
% wy=data(2)*pi/180;
% wz=data(3)*pi/180;
wx=data(1);
wy=data(2);
wz=data(3);

omega=[0  , -wx , -wy , -wz ;...
       wx ,  0  ,  wz , -wy ;...
       wy , -wz ,  0  ,  wx ;...
       wz ,  wy , -wx ,  0   ];

end

function Ang3 = qtpoula(q)
%四元数转欧拉角
    x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
    y = asin(2*(q(1)*q(3) - q(2)*q(4)));
    z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
    Ang3 = [x y z]';
end
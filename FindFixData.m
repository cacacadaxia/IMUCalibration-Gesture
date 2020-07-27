function [PP,fix_point,rotation]=FindFixData(cal,threshold)
%%摘数据的部分我没有按照ICRA2014的算法写，
% 而是自己写了FindFixData这样一个函数。
% ICRA2014的摘数据的算法我也拿matlab写了一遍
% ，算方差实在是太花时间了（可能我的那个算法需要优化），
% FindFixData算法轻便一些，就是要人为的设置参数，ICRA2014
% 不需要人为设置参数。摘出的稳定数据校正角速度参数和磁力计参数，
% 摘出的运动的数据用来校正角速度传感器。
% 
% author  Zhang Xin

n=size(cal,1);
j=1;
for i=1:n
    norm_gyro(i,1)=norm(cal(i,5:7));
    if i==1
        if norm_gyro(i)<threshold
            P(j,1)=i;
        end
    else
        if norm_gyro(i)<=threshold&&norm_gyro(i-1)>threshold
            P(j,1)=i;
        end
        if norm_gyro(i)>threshold&&norm_gyro(i-1)<=threshold
            P(j,2)=i-1;
            j=j+1;
        end
        
    end
end
j=1;
for i=1:size(P,1)-1
    if P(i,2)-P(i,1)>20
        PP(j,1)=P(i,1);
        PP(j,2)=P(i,2);        
        fix_point(j,:)=mean(cal(PP(j,1):PP(j,2),2:end),1);
        if j>=2
            rotation{j-1,1}=cal(PP(j-1,2)-30:PP(j,1)+30,:);
        end
        j=j+1;
    end
end


figure 

plot(1:n,norm_gyro,'b')
for j=1:size(PP,1)
   hold on
   plot(PP(j,1),norm_gyro(PP(j,1)),'ro');
   hold on
   plot(PP(j,2),norm_gyro(PP(j,2)),'ko');
end


end

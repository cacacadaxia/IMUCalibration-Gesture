function [PP,fix_point,rotation]=FindFixData(cal,threshold)
%%ժ���ݵĲ�����û�а���ICRA2014���㷨д��
% �����Լ�д��FindFixData����һ��������
% ICRA2014��ժ���ݵ��㷨��Ҳ��matlabд��һ��
% ���㷽��ʵ����̫��ʱ���ˣ������ҵ��Ǹ��㷨��Ҫ�Ż�����
% FindFixData�㷨���һЩ������Ҫ��Ϊ�����ò�����ICRA2014
% ����Ҫ��Ϊ���ò�����ժ�����ȶ�����У�����ٶȲ����ʹ����Ʋ�����
% ժ�����˶�����������У�����ٶȴ�������
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

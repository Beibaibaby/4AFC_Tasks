Mydata=[result.this_ori;result.this_dur;result.subres;result.RT]';
%����ÿ��duration����ȷ��Ӧ�ʡ�����Ӧ�ʣ�ͬ�����췽��
durs=[20,40,60,80,100,120,140]/1000;
a=[1,2,3,4,
    2,1,4,3,
    3,4,1,2
    4,3,2,1];%ÿ�зֱ����corr,incorr_0,incorr_1,incorr_2
for m=1:length(durs)
    mydata=Mydata(Mydata(:,2)==durs(m),:);%ÿ��duration������
    t=length(mydata);%ÿ��duration�����Դ�
    %��ȷ��Ӧ����cor
    corr=[0,0,0,0]  %corr=0;incorr_0=0;incorr_1=0;incorr_2=0;
    for n=1:t
        c=mydata(n,1);
        r=mydata(n,3);
        w=find(a(:,c)==r);
        corr(w)=corr(w)+1;
    end
    corr_rate(m,:)=corr/t;
end
figure;
col=['r','g','y','m'];
for p=1:4
    scatter([0,durs],[0.25,corr_rate(:,p)'],col(p),'filled');
    hold on
end
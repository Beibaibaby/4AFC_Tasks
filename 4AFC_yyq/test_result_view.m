load('data\20220405004032test_subID-280.mat');
Mydata=[result.this_ori;result.this_dur;result.subres;result.RT]';
%计算每个duration的正确反应率、错误反应率（同方向、异方向）
durs=[20,40,60,80,100,120,140]/1000;
a=[1,2,3,4,
    2,1,4,3,
    3,4,1,2
    4,3,2,1];%每行分别代表：corr,incorr_0,incorr_1,incorr_2
react_num=zeros(length(durs),4);
%反应数据react_num :  corr;incorr_0,incorr_1,incorr_2
for m=1:length(durs)
    mydata=Mydata(Mydata(:,2)==durs(m),:);%每个duration的数据
    t_num(m)=length(mydata);%每个duration的总试次
    for n=1:t_num(m)
        c=mydata(n,1);
        r=mydata(n,3);
        w=find(a(:,c)==r);
        react_num(m,w)=react_num(m,w)+1;
    end
    corr_rate(m,:)=react_num(m,:)/t_num(m);
end

%fit:logistic;mle
PF=@PAL_Logistic;
paramsFree=[1,1,0,0];
searchGrid.alpha=[20:100]/1000;
searchGrid.beta=10.^[-1:.01:2];
searchGrid.gamma=0.25;
searchGrid.lambda=0.02;

[paramsValues LL exitflag]=PAL_PFML_Fit(durs,react_num(:,1)',t_num,...
    searchGrid,paramsFree,PF)
%paramsValues:parameter estimates that define the best-fitting PF.
%LL:log likelihood associated with the fit.
%exitflag:the fit was successful.

%%% graph
StimLevelsFine=[20:0.1:140]/1000;
Fit=PF(paramsValues,StimLevelsFine);
plot(StimLevelsFine,Fit,'b-','linewidth',2);
hold on;
set(gca,'fontsize',12);
axis([0,0.15,0,1]);

style={'b.','y.','g.','m.'};
for p=1:4
    plot(durs,corr_rate(:,p)',style{p},'markersize',10);
    hold on
end

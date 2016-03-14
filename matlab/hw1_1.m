Returns=xlsread('sizeDecileReturns.xls');
DReturns=Returns(2:end,:);
nDReturns=DReturns(:,2:11);
percentage = '%4.2f %% %4.2f%% %4.2f%% %4.2f%% %4.2f%% %4.2f%% %4.2f%% %4.2f%% %4.2f%% %4.2f%%\n';
fid = fopen('results.txt','w');
pnDReturns=100*nDReturns;
% ------------------------- 1.a-----------------------------------------------------
meanr=mean(pnDReturns);
lowestr=min(pnDReturns);
highestr=max(pnDReturns);

fprintf(fid, '1.a\nThe mean returns :\n');
  fprintf(fid, percentage,meanr);
fprintf(fid, 'The lowest returns :\n');
  fprintf(fid, percentage,lowestr);
fprintf(fid, 'The highest returns :\n');
  fprintf(fid, percentage,highestr);
% ------------------------- 1.b-----------------------------------------------------
[RSorted,index] = sort(nDReturns);
indexLarge1=index(end,1);
rLarge1=100*RSorted(end,1);
[rAbsSort,absindex] = sort(abs(nDReturns));
indexSmall10=absindex(1,10);
rSmall10=100*rAbsSort(1,10);

fprintf(fid, '\n1.b\n');
fprintf(fid, 'index of the largest return :%d\n',indexLarge1);
fprintf(fid, 'percentage return :%4.2f%%\n',rLarge1);
fprintf(fid, 'index of the smallest absolute return in decile 10 :%4.2f\n',indexSmall10);
fprintf(fid, 'percentage return :%4.2f%%\n',rSmall10);
% ------------------------- 1.c-----------------------------------------------------
rho=corr(nDReturns(:,1:9),nDReturns(:,10));

fprintf(fid, '\n1.c\n');
fprintf(fid, 'correlations:\n');
  fprintf(fid, '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n',rho);
% ------------------------- 1.d-----------------------------------------------------
NBstd=std(pnDReturns.^2,1);
Bstd=std(pnDReturns,0,2);
 
fprintf(fid, '\n1.d\n');
fprintf(fid, 'The standard deviations of the squared percentage returns :\n');
  fprintf(fid, '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n',NBstd);
fprintf(fid, 'The standard deviation for the first day in the sample : %4.2f' ,Bstd(1));
% ------------------------- 1.e-----------------------------------------------------
mu=0;  
sigma2=1;
x=0;
logLike=normal_log_likelihood(x,mu,sigma2);  
  
fprintf(fid, '\n1.e\n');
fprintf(fid, 'The log likelihood :%4.2f\n',logLike);
% ------------------------- 1.f-----------------------------------------------------
mu5=mean(pnDReturns(:,5));  
sigma25=var(pnDReturns(:,5));
x5=pnDReturns(:,5);
logLike5=modified_normal_log_likelihood(x5,mu5,sigma25);  
  
fprintf(fid, '\n1.f\n');
fprintf(fid, 'The log likelihood for decile 5 :%4.2f\n',logLike5);
% ------------------------- 2.a-----------------------------------------------------
T=10000;
e=randn(T,1);
Y=zeros(T,1);
Y(1)=e(1);
alpha=0.1;
for t=2:T;
  Y(t)=(1-alpha)*Y(t-1)+e(t);
end
plot(Y);
xlabel('Time Trend');
ylabel('Y(t)');
title('Plot of Y(t)-1');
saveas(gcf,'Plot of Y(t)-1.jpeg');
% Compute and save kurtosis
ka=kurtosis(Y);


% ------------------------- 2.b-----------------------------------------------------
e=randn(T,1);
% Initialize a place for Y
Y=zeros(T,1);
% Use initial condition
Y(1)=e(1);
% Loop over t
for t=2:T;
  if(e(t-1) < 0)
    Y(t)=(1-alpha)*Y(t-1)+2*e(t);
  else
    Y(t)=(1-alpha)*Y(t-1)+e(t);
  end
end
% Plot Y
plot(Y);
xlabel('Time Trend');
ylabel('Y(t)');
title('Plot of Y(t)-2');
saveas(gcf,'Plot of Y(t)-2.jpeg');
% Compute and save kurtosis
kb=kurtosis(Y);

% ------------------------- 2.c-----------------------------------------------------
fprintf(fid, '\n2.c\n');
fprintf(fid, 'The kurtosis values for Y(t)--1: %4.2f and the kurtosis values for Y(t)--1: %4.2f\n',ka,kb);
fprintf(fid, 'The y(t) process in part a will have a kurtosis that converge to 3, because y(1) is a normal ditribution and e(1) is a normal distribution, so, y(2) is also a normal distribution, and so on we can prove that y(t) is a normal distribution, so the kurtosis converge to 3.\n');
% ------------------------- 3.a-----------------------------------------------------
n1=sum(nDReturns(:,5)<0);
n2=sum(nDReturns(:,5)==0);
n3=sum(nDReturns(:,5)>0);

fprintf(fid, '\n3.a\n');
fprintf(fid, 'Number of negative :%d , zero:%d,   positive: %d \n',n1,n2,n3);
% ------------------------- 3.b-----------------------------------------------------
n4=sum(abs(nDReturns(:,5))>2*std(nDReturns(:,5)));  

fprintf(fid, '\n3.b\n');
fprintf(fid, 'Number of abs returns greater than 2*std(r):%d \n',n4);
% ------------------------- 3.c-----------------------------------------------------
n5=find(nDReturns(:,5)<0);
mn5=mean(nDReturns(n5,5));
sn5=std(nDReturns(n5,5));

fprintf(fid, '\n3.c\n');
fprintf(fid, 'Mean :%4.2f , standard deviation: %4.2f \n',mn5,sn5);
% ------------------------- 3.d-----------------------------------------------------
n6=sum(all(nDReturns>0,2));
n7=sum(any(nDReturns>0,2));
n8=sum(~any(nDReturns>0,2));

fprintf(fid, '\n3.d\n');
fprintf(fid, 'all positive: %d at least one positive: %d no positive: %d\n',n6,n7,n8);
fclose(fid);
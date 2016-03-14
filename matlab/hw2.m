fid = fopen('results.txt','w');
% ------------------------- 1.a-----------------------------------------------------
%build 6 random array
var5n=var(randn(5,1000));
var20n=var(randn(20,1000));
var25n=var(randn(25,1000));
var50n=var(randn(50,1000));
var100n=var(randn(100,1000));
var500n=var(randn(500,1000));
%plot them
subplot(3,2,1),hist(var5n);
title('sample variance when T=5');
subplot(3,2,2),hist(var20n);
title('sample variance when T=20');
subplot(3,2,3),hist(var25n);
title('sample variance when T=25');
subplot(3,2,4),hist(var50n);
title('sample variance when T=50');
subplot(3,2,5),hist(var100n);
title('sample variance when T=100');
subplot(3,2,6),hist(var500n);
title('sample variance when T=500');
saveas(gcf,'Q1_a.pdf');
% ------------------------- 1.b-----------------------------------------------------
clt5n=sqrt(5)*(var5n-1)/sqrt(2);
clt20n=sqrt(20)*(var20n-1)/sqrt(2);
clt25n=sqrt(25)*(var25n-1)/sqrt(2);
clt50n=sqrt(50)*(var50n-1)/sqrt(2);
clt100n=sqrt(100)*(var100n-1)/sqrt(2);
clt500n=sqrt(500)*(var500n-1)/sqrt(2);
subplot(3,2,1),hist(clt5n);
title('CLT for normal when T=5');
subplot(3,2,2),hist(clt20n);
title('CLT for normal when T=20');
subplot(3,2,3),hist(clt25n);
title('CLT for normal when T=25');
subplot(3,2,4),hist(clt50n);
title('CLT for normal when T=50');
subplot(3,2,5),hist(clt100n);
title('CLT for normal when T=100');
subplot(3,2,6),hist(clt500n);
title('CLT for normal when T=500');
saveas(gcf,'Q1_b.pdf');
% ------------------------- 1.c-----------------------------------------------------
%build 6 random array
var5c=var(chi2rnd(1,5,1000));
var20c=var(chi2rnd(1,20,1000));
var25c=var(chi2rnd(1,25,1000));
var50c=var(chi2rnd(1,50,1000));
var100c=var(chi2rnd(1,100,1000));
var500c=var(chi2rnd(1,500,1000));
%plot them
clt5c=sqrt(5)*(var5c-2)/sqrt(56);
clt20c=sqrt(20)*(var20c-2)/sqrt(56);
clt25c=sqrt(25)*(var25c-2)/sqrt(56);
clt50c=sqrt(50)*(var50c-2)/sqrt(56);
clt100c=sqrt(100)*(var100c-2)/sqrt(56);
clt500c=sqrt(500)*(var500c-2)/sqrt(56);
subplot(3,2,1),hist(clt5c);
title('CLT for chi2 when T=5');
subplot(3,2,2),hist(clt20c);
title('CLT for chi2 when T=20');
subplot(3,2,3),hist(clt25c);
title('CLT for chi2 when T=25');
subplot(3,2,4),hist(clt50c);
title('CLT for chi2 when T=50');
subplot(3,2,5),hist(clt100c);
title('CLT for chi2 when T=100');
subplot(3,2,6),hist(clt500c);
title('CLT for chi2 when T=500');
saveas(gcf,'Q1_c.pdf');
% ------------------------- 2.a-----------------------------------------------------
Portfolios=xlsread('FF_2by3.xls');
Factors=xlsread('FF_Factors.xls');
y=Portfolios(:,2:end);
x=Factors(:,2:4);
r=Factors(:,5);
y=y-[r,r,r,r,r,r];
result1=zeros(4,12);
for i=1:6
    [b,tstat,s2,VCV,VCV_white,R2,Rbar,yhat]=linreg(y(:,i),x,1);
    result1(:,2*i-1:2*i)=[b,tstat];
    fprintf(fid, 'estimations for portfolio %1.0f are %4.2f %4.2f %4.2f %4.2f,\n and the t-statistics are %4.2f %4.2f %4.2f %4.2f\n',i,result1(:,2*i-1),result1(:,2*i));
end
% ------------------------- 2.b-----------------------------------------------------
result2=result1;
for i=1:6
    [b,tstat,s2,VCV,VCV_white,R2,Rbar,yhat]=linreg(y(:,i),x,1);
    tstat=b./diag(sqrt(VCV));
    result2(:,2*i)=tstat;
    fprintf(fid, 'T-statistics using OLS for portfolio %1.0f are %4.2f %4.2f %4.2f %4.2f,\n', i,result2(:,2*i));
end
% ------------------------- 2.c-----------------------------------------------------
result3=zeros(1,6);
for i=1:6
    [b,tstat,s2,VCV,VCV_white,R2,Rbar,yhat]=linreg(y(:,i),x,1);
    result3(:,i)=R2;
    fprintf(fid, 'R-square for portfolio %1.0f is %4.2f ,\n', i,result3(1,i));
end
% ------------------------- 2.d-----------------------------------------------------

fclose(fid);






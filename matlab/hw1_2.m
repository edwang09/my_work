% Question 2.a
% Set T 
T=10000;
% Create some pseudo random normal numbers
e=randn(T,1);
Y=zeros(T,1);
% Use initial condition
Y(1)=e(1);
alpha=0.1;
for t=2:T;
  Y(t)=(1-alpha)*Y(t-1)+e(t);
end
% Plot Y 
plot(Y);
xlabel('Time Trend');
ylabel('Y(t)');
title('Plot of Y(t)-1');
saveas(gcf,'Q2_a.pdf');
% Compute and save kurtosis
ka=kurtosis(Y);


% Question 2.b
%T=10000;
% Create some pseudo random numbers

%rng(0);
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
xlabel('Time Index');
ylabel('Y value');
title('Plot of Y(t)-2');
saveas(gcf,'Q2_b.pdf');
% Compute and save kurtosis
kb=kurtosis(Y);

% Question 2.c
% Print kurtosis values
fprintf(fid, '\r\n\r\nQuestion 2.c\r\n');
fprintf(fid, '\r\nThe kurtosis values for Y(t)--1 is %4.2f and the kurtosis values for Y(t)--1 is %4.2f',ka,kb);
 % fprintf(fid, percentage,ka,kb);

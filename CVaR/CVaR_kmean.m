%CVaR
%The k-means algorithm (kmeans2 function) is provided by Adreas Grothey(2007) 
%and peter richtarik (2013)
clear all;
load ('retm_ascii.mat',"-ASCII");

%Choose different size of original scenarios
retm_ascii = retm_ascii(1:100000,:);

%record beginning time 
t1= clock;

%k-mena function
[J, mu, c] = kmeans2(retm_ascii,100);

%Obtain the good 100 scenarios
ret = mu;

% n is number of data points (scenarios)
[n,m] = size(ret); 
B = 10000;
L = 10100;
a = 0.95;

%solve sub=problem with scenario subset
cvx_begin
  variables x(8);
  variables sfy(n);
  variable t(1);

  minimize (t + sum(sfy)/(n*(1-a)));
  subject to
  B*ret*x +sfy >= L -t;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end


%out-of-sample
oos = t + mean(max(L - B*retm_ascii*x-t, 0))/(1-a);

%calculate total time
t2= clock;
t= t2-t1;
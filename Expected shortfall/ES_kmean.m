%Expected shortfall
%The k-means algorithm (kmeans2 function) is provided by Adreas Grothey(2007) 
%and peter richtarik (2013)
clear all;
load ('retm_ascii.mat',"-ASCII");

%Choose different size of original scenarios
retm_ascii=retm_ascii(1:100000,:);

%record beginning time 
t1=clock;
%k-mena function
[J, mu, c] = kmeans2(retm_ascii,100);

%Obtain the good 100 scenarios
ret = mu;

% n is number of data points (scenarios)
[n,m] = size(ret); 
B = 10000;
L = 9640;

% solve sub problem with sceanrio subset
cvx_begin
  variables x(8);
  variables sfy(n);
 
  minimize (norm(sfy,1)/n);
  subject to
  B*ret*x +sfy >= L;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end

%out-of-sample
oos= mean(max(L - B*retm_ascii*x,0));

%calculate total time
t2=clock;
t=t2-t1;
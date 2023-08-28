%Randomly choose scenarios
clear all;
load ('retm_ascii.mat',"-ASCII");

%record beginning time
t1=clock;
%Choose different size of original scenarios
retm_ascii = retm_ascii(1:10000,:);

%choose scenario subset randomly
n_subset = 100;
[n_total,~] = size(retm_ascii); 
N_AS = randperm(n_total,n_subset);

%Obtain returen matrix for selected subset scenario
ret = retm_ascii(N_AS,:);
% n_total is number of selected scenario subset
[n,~] = size(ret); 


%parameter
B = 10000;
L = 9640;
 
%Solve opt
cvx_begin
  variables x(8);
  variables sfy(n);

  minimize (sum(sfy)/n);
  subject to
  B*ret*x +sfy >= L;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end

%calculate out-of-sample
oos = mean(max(L-B*retm_ascii*x,0));

%calculate total time
t2=clock;
t=t2-t1;

%creat a matrix and Sort it according to all loss
%re_sfy = max(round(L - B*ret*x,5),0);
%data_all_sfy = [re_sfy ret];
%sort_all_sfy = sortrows(data_all_sfy,1);
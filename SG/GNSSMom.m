function [w mse]= gnss_mom(X,d, method)

% X -> feature
% d -> target
% method = 0, LMS regular
%        = 1, LMS com momentum 
%        = 2, Nesterov


max_iter = 100;

[m N] = size(X);

X = [ones(1, N); X];
lr = 0.00000025;

w = ones(m + 1, 1);
w=w/2;
z=zeros(m + 1, 1);

n = 1; %epoch
mse = inf;
max_mse = 0.00003;

beta=.95 ; %.95 %Nesterov 0.97 
if method ==2
    lambda_k=1;
    beta=.97;
end

while mse > max_mse & n <= max_iter

  for i=1:N
    e = d(1,i) - w(:,length(w(1,:)))'*X(:,i);
    if method ==1 
      z=beta*z - lr*e*X(:,i);
      w = w - z;
    elseif method == 2  %Nesterov 
      e = d(1,i) - (w(:,length(w(1,:)))- 1. *beta*z)'*X(:,i); % deu melhor com 5. Para maior que 5 piora mais ou menos da mesma maneira que para menor que 5
      z=beta*z - lr*e*X(:,i);
      w =w - z;
    else % SGD
      w = w + lr*e*X(:,i);
    end
  end

  %.
  % Calculate mean squared error
  %

  e_tot = 0;
  for i = 1:N
    e_tot = e_tot + (d(1,i) - w(:,length(w(1,:)))'*X(:,i))^2;
  end
  mse(n) = e_tot / N;
  n=n+1
      
end









function [f_prime_x0,f_second_x0] = RichardsonDer(f,x0,h_in,N)
% Function RICHARDSONDER implements the Richardson extrapolation algorithm 
% for the approximation of the 1st and 2nd derivative of the single-valued 
% real function F in the point X0, given an initial step size H_IN for the 
% central difference formulas and the order of the extrapolation N.
%
% The outputs are the two tables of extrapolates for the 1st and 2nd
% derivative F_PRIME_X0 and F_SECOND_X0. The first column of each table
% is the array of different step values used, while the other columns contain 
% the derivative extrapolates; each column starting from the 3rd one has a 
% lower order of the error for fixed h than the one on its left (0 entries 
% in the tables are not-computed values). The best approximation is thus
% the last element of each table (least error and smallest step size).
%
%
% Example: Richardson derivative extrapolation using function RICHARDSONDER
% % Function and derivatives
% f=@(x) x*exp(x)
% fp=@(x) exp(x)+x*exp(x);
% fs=@(x) 2*exp(x)+x*exp(x);
% % Point of evaluation
% x0 = 2
% % Exact values of the derivatives
% f_prime_x0 = fp(x0)
% f_second_x0= fs(x0)
% % Step initial value
% h_in = 0.4
% % Order of extrapolation 
% N = 4
% % Calling function RichardsonDer
% [f_prime_x0,f_second_x0] = RichardsonDer(f,x0,h_in,N)

% Step values array
h=ones(N,1);
h(1)=h_in;
for j=2:N
    h(j)=h(j-1)/2;
end

% Preallocating tables for 1st and 2nd derivative
f_prime_x0=zeros(N,N+1);
f_second_x0=f_prime_x0;

% 1st column (step values)
f_prime_x0(:,1)=h;
f_second_x0(:,1)=h;

% 2nd column (central difference formulas for 1st and 2nd derivative)
for j=1:N
    f_prime_x0(j,2)=f(x0+h(j))/(2*h(j))-f(x0-h(j))/(2*h(j));
    f_second_x0(j,2)=f(x0+h(j))/h(j)^2-2*f(x0)/h(j)^2+f(x0-h(j))/h(j)^2;
end

% Other columns (Richardson extrapolates)
for i=3:N+1
    for j=i-1:N
        f_prime_x0(j,i)=(2^(2*(i-2))*f_prime_x0(j,i-1)-f_prime_x0(j-1,i-1))/...
            (2^(2*(i-2))-1);
        f_second_x0(j,i)=(2^(2*(i-2))*f_second_x0(j,i-1)-f_second_x0(j-1,i-1))/...
            (2^(2*(i-2))-1);
    end
end

end
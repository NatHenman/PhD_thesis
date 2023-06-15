clear all

% solution parameters
Re = 1;     % Reynolds number
lambda = 1; % Slip length magnitude

% simulation parameters
dX = 0.005;
dY = 0.01;
tol = 1e-6;
X_end = sqrt(5);

% solution spaces
X = 0:dX:X_end;
Y = 0:dY:1;
u = zeros(length(Y),length(X));
v = zeros(length(Y),length(X));
h = zeros(1,length(X));
J = zeros(length(Y),length(Y)); % matrix of coefficients for Neton algorithm
f = zeros(length(Y),1);         % Newton algorithm: f=0
w = zeros(length(Y),1);         % Newton algorithm correction vector
Jv = zeros(length(Y),length(Y));% constant coefficient matrix for second derivative

% define constant coefficients of Jv
for i = 2:length(Y)-1
    Jv(i,i-1) = 1/dY^2; 
    Jv(i,i) = -2/dY^2; 
    Jv(i,i+1) = 1/dY^2;
end

% x = 0 condition
u(:,1) = 3*(2*lambda+2*Y-Y.^2)/2/(1+3*lambda);
h(1) = 1;

% start marching in X
for i = 2:length(X)
    
    % initial guess form previous X station
    h(i) = h(i-1);
    u(:,i) = u(:,i-1);
    v(:,i) = v(:,i-1);
    
    % start local iteration
    uer = 10; ver = 10; her = 10;
    while uer > tol || ver > tol || her > tol
        
        u_old = u(:,i);
        v_old = v(:,i);
        h_old = h(i);
        
        % calculate u
        wer = 10;
        while wer > tol
            for j = 3:length(Y)-2
                J(j,j-1) = -X(i)*h(i)*v(j,i)/dY - 2*X(i)/Re/dY^2;
                J(j,j) = h(i)^2*(2*u(j,i)-u(j,i-1))/dX + 4*X(i)/Re/dY^2;
                J(j,j+1) = X(i)*h(i)*v(j,i)/dY - 2*X(i)/Re/dY^2;
                f(j) = h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + X(i)*h(i)*v(j,i)*(u(j+1,i)-u(j-1,i))/dY - 2*X(i)*(u(j+1,i)-2*u(j,i)+u(j-1,i))/Re/dY^2;
            end
            for j = length(Y)-1
                J(j,j-1) = -4*X(i)*h(i)*v(j,i)/3/dY - 4*X(i)/3/Re/dY^2;
                J(j,j) = h(i)^2*(2*u(j,i)-u(j,i-1))/dX + 4*X(i)*h(i)*v(j,i)/3/dY + 4*X(i)/3/Re/dY^2;
                f(j) = h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + 4*X(i)*h(i)*v(j,i)*(u(j,i)-u(j-1,i))/3/dY + 4*X(i)*(u(j,i)-u(j-1,i))/3/Re/dY^2;
            end
            for j = 2
                a = lambda*(4*u(j,i)-u(j+1,i))/(2*dY+3*lambda);
                b = 4*lambda/(2*dY+3*lambda);
                J(j,j) = h(i)^2*(2*u(j,i)-u(j,i-1))/dX - X(i)*h(i)*v(j,i)*b/dY - 2*X(i)*(b-2)/Re/dY^2;
                J(j,j+1) = X(i)*h(i)*v(j,i)*(1+b/4)/dY - 2*X(i)*(1-b/4)/Re/dY^2;
                f(j) = h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + X(i)*h(i)*v(j,i)*(u(j+1,i)-a)/dY - 2*X(i)*(u(j+1,i)-2*u(j,i)+a)/Re/dY^2;
            end
            w(2:length(Y)-1) = -J(2:length(Y)-1,2:length(Y)-1) \ f(2:length(Y)-1);
            u(2:length(Y)-1,i) = u(2:length(Y)-1,i) + w(2:length(Y)-1);
            wer = max(abs(w(2:length(Y)-1)))/(sum(abs(u(2:length(Y)-1,i)))/(length(Y)-2));
        end
        u(1,i) = lambda*(4*u(2,i)-u(3,i))/(2*dY+3*lambda);
        u(length(Y),i) = (4*u(length(Y)-1,i)-u(length(Y)-2,i))/3;
        
        %calculate h
        h(i) = (trapz(Y,u(:,i)))^(-1);
        
        %calculate v
        for j = 2:length(Y)-1
            f(j) = -((h(i)-h(i-1))*(u(j+1,i)-u(j-1,i))+h(i)*(u(j+1,i)-u(j+1,i-1)-u(j-1,i)+u(j-1,i-1)))/8/X(i)/dX/dY;
        end
        v(2:length(Y)-1,i) = Jv(2:length(Y)-1,2:length(Y)-1) \ f(2:length(Y)-1);
        
        %errors
        uer = max(abs(u(:,i)-u_old))/(sum(abs(u(:,i)))/length(Y));
        ver = max(abs(v(:,i)-v_old))/(sum(abs(v(:,i)))/length(Y));
        her = abs((h(i)-h_old)/h(i));
        
    end
    
end
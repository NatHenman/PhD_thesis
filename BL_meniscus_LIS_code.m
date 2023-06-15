clear all

% solution parameters
Re = 100;
Ca = 0.01;
epsilon = 0.1;
h0 = 1;
mu = 1;
rho = 1;

% simulation parameters
dX = 0.005;
dY = 0.01;
tol = 1e-6;     % convergence tolerance
Delta = 1e-4;   % small shift parameter in x transform

% solution space
X = 0:dX:1;                     %transformed grid 
Y = 0:dY:1;                     %transformed grid
u = zeros(length(Y),length(X));
v = zeros(length(Y),length(X));
h = zeros(1,length(X));
H = zeros(1,length(X));
P = zeros(1,length(X));
usx = zeros(1,length(X));       % d(u_s)/dx
A = zeros(1,length(X));
B = zeros(1,length(X));
C = zeros(1,length(X));
Ax = zeros(1,length(X));        % dA/dx
Bx = zeros(1,length(X));        % dB/dx
Cx = zeros(1,length(X));        % dC/dx
f0 = zeros(1,length(X));
f0x = zeros(1,length(X));       % d(f_0)/dx
J = zeros(length(Y),length(Y)); % matrix of coefficients for Neton algorithm
f = zeros(length(Y),1);         % Newton algorithm: f=0
w = zeros(length(Y),1);         % Newton algorithm correction vector
Jv = zeros(length(Y),length(Y));% constant coefficient matrix for second derivative
G = zeros(1,length(X));         % H''' = G
gi = zeros(1,length(X));        % int_0^x G
gii = zeros(1,length(X));       % int_0^x gi
giii = zeros(1,length(X));      % int_0^x gii
        

% at X = 0
u(:,1) = 3*Y.*(2-Y)/2;
h(1) = h0;

% define Jv coefficients
for i = 2:length(Y)-1
    Jv(i,i-1) = 1/dY^2; Jv(i,i) = -2/dY^2; Jv(i,i+1) = 1/dY^2;
end

% find x
x = zeros(1,length(X));
myfun = @(xd,c) c-((xd+Delta)^0.5-(1+Delta-xd)^0.5)/2/((1+Delta)^0.5-Delta^0.5)-0.5;  
for i = 1:length(X)
    fun = @(xd) myfun(xd,X(i));
    x(i) = fzero(fun,[0,1]);
end

% Xx = dX/dx
Xx = zeros(1,length(X));
for i = 1:length(X)
    Xx(i) = ((x(i)+Delta)^-0.5+(1+Delta-x(i))^-0.5)/4/((1+Delta)^0.5-Delta^0.5);
end

% f0, A, B, C and derivatives
f0 = (4-erf(2*x/epsilon)+erf(2*(x-1)/epsilon))/3;
f0x = 4*(exp(-4*(x-1).^2/epsilon^2)-exp(-4*x.^2/epsilon^2))/epsilon/pi^0.5;
A = 2*(2-3*f0)./(2*f0-1)./(f0-1);
B = 3*(2*f0.^2-1)./(2*f0-1)./(f0-1);
C = f0.*(3-4*f0)./(2*f0-1)./(f0-1);
Ax = 2*(6*f0.^2-8*f0+3).*f0x./(f0-1).^2./(2*f0-1).^2;
Bx = -3*(6*f0.^2-8*f0+3).*f0x./(f0-1).^2./(2*f0-1).^2;
Cx = (6*f0.^2-8*f0+3).*f0x./(f0-1).^2./(2*f0-1).^2;

% begin global iteration between H and u_s
Her = 10; uer2 = 10; 
while Her > tol || uer2 > tol
    
    H_old = H;
    u_old2 = u(1,:);
    
    % begin marching in X
    for i = 2:length(X)
        
        % guess solutions from previous X station
        u(:,i) = u(:,i-1);
        v(:,i) = v(:,i-1);
        h(i) = h(i-1);
        
        % begin local iteration
        uer = 10; ver = 10; her = 10;
        while uer > tol || ver > tol || her > tol
        
            u_old = u(:,i);
            v_old = v(:,i);
            h_old = h(i);
            
            % find u via Newton algorithm
            wer = 10;
            while wer > tol
                if i ~= length(X) % slip condition when x < 1, no slip when x = 1
                    for j = 3:length(Y)-2
                        J(j,j-1) = -h(i)*v(j,i)/2/dY - 1/epsilon^2/Re/dY^2;
                        J(j,j) = Xx(i)*h(i)^2*(2*u(j,i)-u(j,i-1))/dX + 2/epsilon^2/Re/dY^2;
                        J(j,j+1) = h(i)*v(j,i)/2/dY - 1/epsilon^2/Re/dY^2;
                        f(j) = Xx(i)*h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + h(i)*v(j,i)*(u(j+1,i)-u(j-1,i))/2/dY - (u(j+1,i)-2*u(j,i)+u(j-1,i))/epsilon^2/Re/dY^2;
                    end
                    for j = length(Y)-1
                        J(j,j-1) = -2*h(i)*v(j,i)/3/dY - 2/3/epsilon^2/Re/dY^2;
                        J(j,j) = Xx(i)*h(i)^2*(2*u(j,i)-u(j,i-1))/dX + 2*h(i)*v(j,i)/3/dY + 2/3/epsilon^2/Re/dY^2;
                        f(j) = Xx(i)*h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + 2*h(i)*v(j,i)*(u(j,i)-u(j-1,i))/3/dY + 2*(u(j,i)-u(j-1,i))/3/epsilon^2/Re/dY^2;
                    end
                    for j = 2 % partial slip condition
                        a = (H(i)+1)*(4*u(j,i)-u(j+1,i))/(2*dY*h(i)*(3*A(i)+2*B(i)+C(i))*mu+3*(H(i)+1));
                        b = 4*(H(i)+1)/(2*dY*h(i)*(3*A(i)+2*B(i)+C(i))*mu+3*(H(i)+1));
                        J(j,j) = Xx(i)*h(i)^2*(2*u(j,i)-u(j,i-1))/dX - h(i)*v(j,i)*b/2/dY - (b-2)/epsilon^2/Re/dY^2;
                        J(j,j+1) = h(i)*v(j,i)*(1+b/4)/2/dY - (1-b/4)/epsilon^2/Re/dY^2;
                        f(j) = Xx(i)*h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + h(i)*v(j,i)*(u(j+1,i)-a)/2/dY - (u(j+1,i)-2*u(j,i)+a)/epsilon^2/Re/dY^2;
                    end
                else
                    for j = 2:length(Y)-2
                        J(j,j-1) = -h(i)*v(j,i)/2/dY - 1/epsilon^2/Re/dY^2;
                        J(j,j) = Xx(i)*h(i)^2*(2*u(j,i)-u(j,i-1))/dX + 2/epsilon^2/Re/dY^2;
                        J(j,j+1) = h(i)*v(j,i)/2/dY - 1/epsilon^2/Re/dY^2;
                        f(j) = Xx(i)*h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + h(i)*v(j,i)*(u(j+1,i)-u(j-1,i))/2/dY - (u(j+1,i)-2*u(j,i)+u(j-1,i))/epsilon^2/Re/dY^2;
                    end
                    for j = length(Y)-1
                        J(j,j-1) = -2*h(i)*v(j,i)/3/dY - 2/3/epsilon^2/Re/dY^2;
                        J(j,j) = Xx(i)*h(i)^2*(2*u(j,i)-u(j,i-1))/dX + 2*h(i)*v(j,i)/3/dY + 2/3/epsilon^2/Re/dY^2;
                        f(j) = Xx(i)*h(i)^2*u(j,i)*(u(j,i)-u(j,i-1))/dX + 2*h(i)*v(j,i)*(u(j,i)-u(j-1,i))/3/dY + 2*(u(j,i)-u(j-1,i))/3/epsilon^2/Re/dY^2;
                    end
                end
                w(2:length(Y)-1) = -J(2:length(Y)-1,2:length(Y)-1) \ f(2:length(Y)-1);
                u(2:length(Y)-1,i) = u(2:length(Y)-1,i) + w(2:length(Y)-1);
                wer = max(abs(w))/(sum(abs(u(:,i)))/length(Y));
            end
            u(length(Y),i) = (4*u(length(Y)-1,i)-u(length(Y)-2,i))/3;
            if i ~= length(X)
                u(1,i) = (H(i)+1)*(4*u(2,i)-u(3,i))/(2*dY*h(i)*(3*A(i)+2*B(i)+C(i))*mu+3*(H(i)+1));
            else
                u(1,i) = 0;
            end
            
            % calculate v (recyling f vector)
            for j = 2:length(Y)-1
                f(j) = -Xx(i)*((h(i)-h(i-1))*(u(j+1,i)-u(j-1,i))+h(i)*(u(j+1,i)-u(j+1,i-1)-u(j-1,i)+u(j-1,i-1)))/2/dX/dY;
            end
            v(2:length(Y)-1,i) = Jv(2:length(Y)-1,2:length(Y)-1) \ f(2:length(Y)-1);
            
            % find h
            h(i) = h0/trapz(Y,u(:,i));
        
            % calculate errors of local iteration
            uer = max(abs(u(:,i)-u_old))/(sum(abs(u(:,i)))/length(Y));
            ver = max(abs(v(:,i)-v_old))/(sum(abs(v(:,i)))/length(Y));
            her = abs((h(i)-h_old)/h(i));
                   
        end
        
    end

    % calculate d(u_s)/dx
    for i = 2:length(X)-1
        usx(i) = Xx(i)*(u(1,i+1)-u(1,i-1))/2/dX;
    end
    usx(1) = Xx(1)*(-3*u(1,1)+4*u(1,2)-u(1,3))/2/dX;
    usx(length(X)) = Xx(length(X))*(3*u(1,length(X))-4*u(1,length(X)-1)+u(1,length(X)-2))/2/dX;
    
    % calculate H
    G = -2*mu*Ca*u(1,:).*(3*A+B)/epsilon^3./(H+1).^2+rho*Re*Ca*u(1,:).*(usx.*(A+B+C)+u(1,:).*(Ax+Bx+Cx))/epsilon;
    for i = 3:length(X)-1
        gi(i) = trapz(x(2:i),G(2:i));
    end
    for i = 2:length(X)
        gii(i) = trapz(x(1:i),gi(1:i));
    end
    for i = 2:length(X)
        giii(i) = trapz(x(1:i),gii(1:i));
    end
    H = giii + 6*x.*(x-1)*trapz(x,giii)+x.*(2-3*x)*giii(length(giii));
    
    % calculate global iteration errors
    uer2 = max(abs(u(1,:)-u_old2))/(sum(abs(u(1,:)))/(length(X)));
    Her = max(abs(H-H_old))/(sum(abs(H(:)))/(length(X)));
  
    % terminates code if H deformation is too large
    if min(H)<=-1
        disp("collapse")
        break
    end
   
end
 
% calculate the pressure
P = zeros(1,length(X));
P = -epsilon^3*(gi+12*trapz(x,giii)-6*giii(length(giii)))/Ca;
P(length(x))=NaN;

% calculate angles meniscus makes at x = 0,1
t1 = pi/2-atan(epsilon*(-6*trapz(x,giii)+2*giii(length(x))));
t2 = pi/2+atan(epsilon*(gii(length(x))+6*trapz(x,giii)-4*giii(length(x))));
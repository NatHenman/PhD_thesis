clear all

% solution parameters
Re = 100;

% simulation parametrs
dtheta = pi/150; dxi = dtheta;  % spatial discretisation steps
dt = 0.005;                     % time step
t_start = 0; t_end = 3;         % start and end times
xi_infinity = 5;                % infinite xi boundary
N = 40;                         % number of coeicients in Fourier series

%solution arrays
xi = 0:dxi:xi_infinity;
theta = 0:dtheta:pi;
t = t_start:dt:t_end;
zeta = zeros(length(xi),length(theta),length(t));
psi = zeros(length(xi),length(theta),length(t));
fn = zeros(N,length(xi));                           % psi Fourier coefficients
Fn = zeros(N,length(xi));                           % F_n(xi) = int_0^pi zeta*sin(n*theta) d(theta)
J = zeros(length(xi),length(xi));                   % matrix coefficients to solve fn ode
K = zeros(length(xi),1);                            % J f_n = K

% find initial psi with zeta(1,:,:) = 0
fn(1,:) = exp(xi)-exp(-xi);

% find psi(:,:,1)
psi_sum = zeros(N,length(xi),length(theta));
for n = 1:N
    for i = 1:length(xi)
        for j = 1:length(theta)
            psi_sum(n,i,j) = fn(n,i)*sin(n*theta(j));
        end
    end
end
for i = 1:length(xi)
    for j = 2:length(theta)-1
        psi(i,j,1) = sum(psi_sum(:,i,j));
    end
end

% find zeta(:,:,2)
for j = 2:length(theta)-1
    zeta(1,j,2) = (psi(3,j,1)-8*psi(2,j,1))/2/dxi^2;
end
for i = 2:length(xi)-1
    for j = 2:length(theta)-1
        zeta(i,j,2) = zeta(i,j,1) + (3*dt*exp(-2*xi(i))/2)*(Re^(-1)*((zeta(i+1,j,1)-2*zeta(i,j,1)+zeta(i-1,j,1))/dxi^2+(zeta(i,j+1,1)-2*zeta(i,j,1)+zeta(i,j-1,1))/dtheta^2) - (psi(i,j+1,1)-psi(i,j-1,1))*(zeta(i+1,j,1)-zeta(i-1,j,1))/4/dtheta/dxi + (psi(i+1,j,1)-psi(i-1,j,1))*(zeta(i,j+1,1)-zeta(i,j-1,1))/4/dtheta/dxi);
    end
end


% start time advancement
for m = 2:length(t)
    
    % find Fn
    for n = 1:N
        for i = 1:length(xi)
            Fn(n,i) = 2*trapz(theta,zeta(i,:,m).*sin(n*theta))/pi;
        end
    end
    
    % calculate fn
    fn(1,length(xi)) = 2*sinh(xi(length(xi)));
    for n = 1:N
        for i = 2:length(xi)-1
            J(i,i+1) = 1;
            J(i,i) = -2-n^2*dxi^2;
            J(i,i-1) = 1;
            K(i) = -dxi^2*exp(2*xi(i))*Fn(n,i);
        end
        J(1,1) = 1;
        J(length(xi),length(xi)) = 1;
        if n == 1
            K(length(xi)) = 2*sinh(xi(length(xi)));
        else
            K(length(xi)) = 0;
        end        
        fn(n,:) = J \ K;
    end
    
    % calculate psi
    for n = 1:N
        for i = 1:length(xi)
            for j = 1:length(theta)
                psi_sum(n,i,j) = fn(n,i)*sin(n*theta(j));
            end
        end
    end
    for i = 1:length(xi)
        for j = 2:length(theta)-1
            psi(i,j,m) = sum(psi_sum(:,i,j));
        end
    end
    
    % calculate zeta
    if m ~= length(t)
        for j = 2:length(theta)-1
            zeta(1,j,m+1) = (psi(3,j,m)-8*psi(2,j,m))/2/dxi^2;
        end
        for i = 2:length(xi)-1
            for j = 2:length(theta)-1
                zeta(i,j,m+1) = zeta(i,j,m) + exp(-2*xi(i))*dt*(3*(Re^(-1)*((zeta(i+1,j,m)-2*zeta(i,j,m)+zeta(i-1,j,m))/dxi^2+(zeta(i,j+1,m)-2*zeta(i,j,m)+zeta(i,j-1,m))/dtheta^2) - (psi(i,j+1,m)-psi(i,j-1,m))*(zeta(i+1,j,m)-zeta(i-1,j,m))/4/dtheta/dxi + (psi(i+1,j,m)-psi(i-1,j,m))*(zeta(i,j+1,m)-zeta(i,j-1,m))/4/dtheta/dxi) - (Re^(-1)*((zeta(i+1,j,m-1)-2*zeta(i,j,m-1)+zeta(i-1,j,m-1))/dxi^2+(zeta(i,j+1,m-1)-2*zeta(i,j,m-1)+zeta(i,j-1,m-1))/dtheta^2) - (psi(i,j+1,m-1)-psi(i,j-1,m-1))*(zeta(i+1,j,m-1)-zeta(i-1,j,m-1))/4/dtheta/dxi + (psi(i+1,j,m-1)-psi(i-1,j,m-1))*(zeta(i,j+1,m-1)-zeta(i,j-1,m-1))/4/dtheta/dxi))/2;
            end
        end
    end
   
end

% calculate zeta at r = 1
zeta_b = zeros(length(theta),length(t));
for m = 1:length(t)
    zeta_b(:,m) = zeta(1,:,m);
end

% calculate d(zeta)/dr at r = 1
zeta_r_b = zeros(length(theta),length(t));
for m = 1:length(t)
    zeta_r_b(:,m) = (-3*zeta(1,:,m)+4*zeta(2,:,m)-zeta(3,:,m))/2/dxi;
end

% calculate the pressure at r = 1
pgr1 = zeros(length(theta),length(t));
for m = 1:length(t)
    for j = 2:length(theta)
        pgr1(j,m) = trapz(theta(1:j),zeta_r_b(1:j,m))/Re;
    end
end

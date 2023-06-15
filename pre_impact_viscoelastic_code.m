clear all

% solution parameters
e3 = -1;
e5 = -1;

% simulation parameters
x1 = -20; x2 = 20;          % spatial domain
t1 = -50; t2 = 200;         % start and end times
dt = 0.005;                 % time step
Nx = 801;                   % number of x points
Con_tol = 1e-6;             % global convergence tolerance
tolu = 1e-6;                % Newton algorithm tolerance
touchdown = 0.1;            % termination point
HilbertCoefficients = 10*Nx;% number of Hilbert Coefficients for Hilbert transform

dx = (x2-x1)/(Nx-1);

% solution spaces
x = x1:dx:x2;           
t = t1:dt:t2;
F = zeros(Nx,1);
P = zeros(Nx,1);
G = zeros(Nx,1);
Dmin = zeros(1,1);      % minimum air film thickness
J = zeros(Nx,Nx);       % Newton algorithm matrix
u = zeros(Nx,1);        % Newton algorithm correction vector
HilbPx = zeros(Nx,1);   % Hilbert transform of dP/dX
Px = zeros(Nx,1);       % dP/dX
f = zeros(Nx,1);        % Newton algorithm f = 0

Nt = length(t);

% non-zero initial condition
F(:,1) = x.^2/2-t(1);
Dmin(1) = min(F(:,1)-G(:,1));

% start time loop
for m = 1:Nt-1 
    
    % guess from previous time step
    P(:,m+1) = P(:,m);
    G(:,m+1) = G(:,m);
    F(:,m+1) = F(:,m);
    F(1,m+1) = F(1,m) - dt; F(Nx,m+1) = F(Nx,m)-dt;
   
    % begin iteration
    Gerror = 10; Ferror = 10;   
    while Ferror > Con_tol || Gerror > Con_tol
    
        Fprevious = F(:,m+1);
        Gprevious = G(:,m+1);
        
        % begin Newton algorithm to find G
        uerror = 10;
        while uerror > tolu;
            
            for i = 2:Nx-1
                J(i,i-1) = 3*(F(i,m+1)-G(i,m+1))^2*((e3*(G(i+1,m+1)-G(i-1,m+1))/(2*dx)+e5*(G(i+1,m+1)-G(i-1,m+1)-G(i+1,m)+G(i-1,m))/(2*dx*dt))/(2*dx) + (F(i+1,m+1)-F(i-1,m+1)-G(i+1,m+1)+G(i-1,m+1))*(-e3/(2*dx)-e5/(2*dx*dt))/(2*dx)) + (F(i,m+1)-G(i,m+1))^3*(e3/dx^2+e5/(dx^2*dt));
            end
            for i = 2:Nx-1
                J(i,i) = -6*(F(i,m+1)-G(i,m+1))*(F(i+1,m+1)-F(i-1,m+1)-G(i+1,m+1)+G(i-1,m+1))*(e3*(G(i+1,m+1)-G(i-1,m+1))/(2*dx)+e5*(G(i+1,m+1)-G(i-1,m+1)-G(i+1,m)+G(i-1,m))/(2*dx*dt))/(2*dx) - 3*(F(i,m+1)-G(i,m+1))^2*(e3*(G(i+1,m+1)-2*G(i,m+1)+G(i-1,m+1))/dx^2 + e5*(G(i+1,m+1)-2*G(i,m+1)+G(i-1,m+1)-G(i+1,m)+2*G(i,m)-G(i-1,m))/(dx^2*dt)) + (F(i,m+1)-G(i,m+1))^3*(-2*e3/dx^2-2*e5/(dx^2*dt)) + 12/dt;
            end
            for i = 2:Nx-1
                J(i,i+1) = 3*(F(i,m+1)-G(i,m+1))^2*(-(e3*(G(i+1,m+1)-G(i-1,m+1))/(2*dx)+e5*(G(i+1,m+1)-G(i-1,m+1)-G(i+1,m)+G(i-1,m))/(2*dx*dt))/(2*dx) + (F(i+1,m+1)-F(i-1,m+1)-G(i+1,m+1)+G(i-1,m+1))*(e3/(2*dx)+e5/(2*dx*dt))/(2*dx)) + (F(i,m+1)-G(i,m+1))^3*(e3/dx^2+e5/(dx^2*dt));
            end
            for i = 2:Nx-1
                f(i) = 3*(F(i,m+1)-G(i,m+1))^2*(F(i+1,m+1)-F(i-1,m+1)-G(i+1,m+1)+G(i-1,m+1))*(e3*(G(i+1,m+1)-G(i-1,m+1))/(2*dx)+e5*(G(i+1,m+1)-G(i-1,m+1)-G(i+1,m)+G(i-1,m))/(2*dx*dt))/(2*dx) + (F(i,m+1)-G(i,m+1))^3*(e3*(G(i+1,m+1)-2*G(i,m+1)+G(i-1,m+1))/dx^2+e5*(G(i+1,m+1)-2*G(i,m+1)+G(i-1,m+1)-G(i+1,m)+2*G(i,m)-G(i-1,m))/(dx^2*dt)) - 12*(F(i,m+1)-F(i,m)-G(i,m+1)+G(i,m))/dt;
            end
            
            u(2:Nx-1) = -J(2:Nx-1,2:Nx-1) \ f(2:Nx-1);
            G(2:Nx-1,m+1) = G(2:Nx-1,m+1) + u(2:Nx-1);
            uerror = max(abs(u))/(sum(abs(G(:,m+1)))/Nx);
        
        end
        
        % Find P
        P(:,m+1) = e3*G(:,m+1)+e5*(G(:,m+1)-G(:,m))/dt;
        
        % Find F, using F(:,0) = F(:,1) + dt
        if m == 1 
            w = F(:,m)+dt;
        else
            w = F(:,m-1);
        end
        
        % Calculate dP/dX and Hilbert transorm
        for i = 2:Nx-1
            Px(i) = (P(i+1,m+1)-P(i-1,m+1))/(2*dx);
        end
        Px(1) = (-3*P(1,m+1)+4*P(2,m+1)-P(3,m+1))/2/dx;
        Px(Nx) = (3*P(Nx,m+1)-4*P(Nx-1,m+1)+P(Nx-2,m+1))/2/dx;
        HilbPx = imag(hilbert(Px,HilbertCoefficients));
        for i = 1:Nx          
            F(i,m+1) = dt^2*HilbPx(i) + 2*F(i,m) - w(i);
        end
        
        %calculate errors
        Ferror = rme(F(:,m+1),Fprevious,Nx);
        Gerror = rme(G(:,m+1),Gprevious,Nx);
    
    end

    Dmin(m+1) = min(F(:,m+1)-G(:,m+1));
    
    % terminate code if Dmin<touchdown
    if Dmin(m+1) < touchdown
        tdi = m+1;
        break
    end

end

% relative max error function
 function [ r ] = rme(x,y,N)

    r = max(abs(x-y))/(sum(abs(x))/N);

 end
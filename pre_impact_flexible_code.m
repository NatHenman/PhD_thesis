clear all

% solution parameters
e1 = -1;
e2 = 1;
e3 = 0; 
e4 = -1;
e5 = 0;

% simulation parameters
x1 = -10; x2 = 10;          % spatial domain
t1 = -50; t2 = 30;          % start and end times
dt = 0.005;                 % time step
Nx = 401;                   % number of x points
Con_tol = 1e-6;             % global convergence tolerance              
touchdown = 0.26;           % termination point
HilbertCoefficients = 5*Nx; % number of Hilbert Coefficients for Hilbert transform

dx = (x2-x1)/(Nx-1);

% solution spaces
x = x1:dx:x2;
t = t1:dt:t2;
F = zeros(Nx,1);
P = zeros(Nx,1);
G = zeros(Nx,1);
Dmin = zeros(1,1);      % minimum air film thickness
J = zeros(Nx,Nx);       % lubrication eqn: P coeficients
HilbPx = zeros(Nx,1);   % Hilbert transform of dP/dX
Px = zeros(Nx,1);       % dP/dX
f = zeros(Nx,1);        % lubrication eqn: JP = f

Nt = length(t); 

% Flex equation matrix, constant coefficients
A = zeros(Nx,Nx);
for i = 4:Nx-3
    A(i,i-2) = e1/dx^4;
    A(i,i-1) = -4*e1/dx^4+e2/dx^2;
    A(i,i) = 6*e1/dx^4-2*e2/dx^2+e3+e4/dt^2+e5/dt;
    A(i,i+1) = -4*e1/dx^4+e2/dx^2;
    A(i,i+2) = e1/dx^4;
end
i = 3;
A(i,i) = -7*e2/4/dx^2+e3+e4/dt^2+e5/dt;
A(i,i+1) = -4*e1/dx^4+e2/dx^2;
A(i,i+2) = e1/dx^4;
i = Nx-2;
A(i,i-2) = e1/dx^4;
A(i,i-1) = -4*e1/dx^4+e2/dx^2;
A(i,i) = -7*e2/4/dx^2+e3+e4/dt^2+e5/dt;

% initial conditions
F(:,1) = x.^2/2-t(1);
Dmin(1) = min(F(:,1)-G(:,1));

% start time advancement
for m = 1:Nt-1 
    
    % guess from previous time step
    P(:,m+1) = P(:,m);
    G(:,m+1) = G(:,m);
    F(:,m+1) = F(:,m);
    F(1,m+1) = F(1,m) - dt; F(Nx,m+1) = F(Nx,m)-dt;
    
    % begin iteration
    Perror = 10; Gerror = 10; Ferror = 10;
    while Perror > Con_tol || Ferror > Con_tol || Gerror > Con_tol
    
        Pprevious = P(:,m+1);
        Fprevious = F(:,m+1);
        Gprevious = G(:,m+1);
            
        % calculate P   
        for i = 2:Nx-1
            J(i,i-1) = -3*(F(i,m+1)-G(i,m+1))^2*(F(i+1,m+1)-F(i-1,m+1)-G(i+1,m+1)+G(i-1,m+1))/(4*dx^2)+(F(i,m+1)-G(i,m+1))^3/dx^2;
        end
        for i = 1:Nx
            J(i,i) = -2*(F(i,m+1)-G(i,m+1))^3/dx^2;
        end
        for i = 2:Nx-1
            J(i,i+1) = 3*(F(i,m+1)-G(i,m+1))^2*(F(i+1,m+1)-F(i-1,m+1)-G(i+1,m+1)+G(i-1,m+1))/(4*dx^2)+(F(i,m+1)-G(i,m+1))^3/dx^2;
        end
        for i = 2:Nx-1
            f(i) = 12*(F(i,m+1)-F(i,m)-G(i,m+1)+G(i,m))/dt;
        end        
        P(2:Nx-1,m+1) = J(2:Nx-1,2:Nx-1) \ f(2:Nx-1);
        
        % relaxation with previous time step
        P(2:Nx-1,m+1) = 0.5*(P(2:Nx-1,m+1)+P(2:Nx-1,m));
        
        % calculate F, using F(:,0) = F(:,1) + dt
        if m == 1 
            w = F(:,m)+dt;
        else
            w = F(:,m-1);
        end
        for i = 2:Nx-1
            Px(i) = (P(i+1,m+1)-P(i-1,m+1))/(2*dx);
        end
        HilbPx = imag(hilbert(Px,HilbertCoefficients));
        Px(1) = (-3*P(1,m+1)+4*P(2,m+1)-P(3,m+1))/2/dx;
        Px(Nx) = (3*P(Nx,m+1)-4*P(Nx-1,m+1)+P(Nx-2,m+1))/2/dx;
        for i = 1:Nx         
            F(i,m+1) = dt^2*HilbPx(i) + 2*F(i,m) - w(i);
        end
        
        % calculate G, using G(:,0) = 0
        if m == 1
            q = zeros(Nx,1);
        else
            q = G(:,m-1);
        end
        G(3:Nx-2,m+1) = A(3:Nx-2,3:Nx-2)\(P(3:Nx-2,m+1)+e4*(2*G(3:Nx-2,m)-q(3:Nx-2))/dt^2+e5*G(3:Nx-2,m)/dt);
        G(2,m+1) = G(3,m+1)/4;
        G(Nx-1,m+1) = G(Nx-2,m+1)/4;
        
        % calculate errors
        Perror = rme(P(:,m+1),Pprevious,Nx);
        Ferror = rme(F(:,m+1),Fprevious,Nx);
        Gerror = rme(G(:,m+1),Gprevious,Nx);
    
    end
    
    Dmin(m+1) = min(F(:,m+1)-G(:,m+1));
    
    % terminate code if Dmin<touchdown
    if Dmin(m+1) < touchdown
        break
    end
      
end

% relative max error function
function [ r ] = rme(x,y,N)

    r = max(abs(x-y))/(sum(abs(x))/N);
 
end
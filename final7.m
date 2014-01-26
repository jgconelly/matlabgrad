%Define time stepping scheme
%	Crank-Nicolson : k = 0.5

k = 0.5

% Define dispersion and spatial discretization parameters

alpha = .003;
n = 201;
dx = 2/200;

%basically dx = .01

nsteps = .005/.01
dt = .01

% initialise arrays and vectors 
b = zeros(n,1);
u = zeros(n,1);
B = zeros(n,n);

%initialize u

u = zeros(n,1);
x = 0.0;
for i=1:n
	u(i,1) = sin(pi*x) + 0.5*sin(3*pi*x);
	x = x + dx;
end

% set up tridiagonal entries of A matrix
A = zeros(n,n);
for i = 2:n-1
	A(i,i-1) = -k*((alpha/(dx^2))+ (-.5/(2*dx)));
	A(i,i) = 2*k*((1/dt)+(alpha/(dx^2)));
	A(i,i+1) = k*((.5/(2*dx))-(alpha/(dx^2)));
end

%similarily for B

B = zeros(n,n);
for(i = 2:n-1)
	B(i,i-1) = k*((alpha/(dx^2))+ (.5/(2*dx)));
	B(i,i) = 2*k*((1/dt)-(alpha/(dx^2)));
	B(i,i+1) = k*((-.5/(2*dx))+(alpha/(dx^2)));
end
%set up boundary conditions
A(1,1) = 1.0;
A(1,2) = -1.0;
A(n,n) = 1.0;
A(n,n-1) = -1.0;

%error= log(error of some sort)
for i=1:nsteps
    b = B*u;
    u = A\b;
end

plot (u)
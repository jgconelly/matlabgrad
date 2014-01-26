% Clear old data
clear all
k = 0.5;

% Define dispersion and spatial discretization parameters
alpha = 0.003;
n = 201;
dx = 2/(n-1);

% Loop through all the temporal discretizations required
dt = [0.01]
nsteps = 3.5/.01;
for rum=1:5
  lambda(rum) = alpha*dt(rum)/dx^2;
  lambda = (0.003*0.01)/(dx^2);
  % Initialize arrays and vectors
  b = zeros(n,1);
  u = zeros(n,1);
  B = zeros(n,n);

  % Initialize u
  u = zeros(n,1);
  u(1,1) = 0.0;
  x = dx;
  for i=2:n-1
    u(i,1) = 2*exp(-(x-10)^2);
    x = x + dx;
  end
  u(n,1) = 0.0;

  % Set up tridiagonal entries of A matrix
  A = zeros(n,n);
  for i=2:n-1
    A(i,i-1) = -k*lambda(rum);
    A(i,i) = 1 + 2*k*lambda(rum);
    A(i,i+1) = -k*lambda(rum);
  end

  % Set up tridiagonal entries of B matrix
  B = zeros(n,n);
  for i=2:n-1
    B(i,i-1) = (1-k)*lambda(rum);
    B(i,i) = 1 - 2*(1-k)*lambda(rum);
    B(i,i+1) = (1-k)*lambda(rum);
  end

  % Set up boundary conditions
  A(1,1) = 1.0;
  B(1,1) = 1.0;
  u(1,1) = 0.0;
  A(n,n) = 1.0;
  A(n,n-1) = -1.0;
  B(n,n) = 1.0;
  u(n,1) = 0.0;

  % Solve for each time step
  for i=1:nsteps
    b = B*u;
    u = A\b;
  end
end
  %error(run) = log(norm(u-soln));
%end

% Plot log(error) vs lambda
plot(lambda,'rx-');
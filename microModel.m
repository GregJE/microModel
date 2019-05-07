% microModel
% Greg Evans - 2018

% microModel:
% Calculates and Represents the diaphragm deformation in a microphone, 
% given a acoustic pressure and tension force. This excludes the
% electrostatic interaction. Deformation is imposed as zero on boundaries.

% Created as part of a Bachelors Degree (Hons) in Mechanical Engineering.

% The following script solves the governing second-order differential 
% equation through both analytical and numerical solution methodologies.
% Equation: -T(d^2f/dx^2+d^f/dy^2) = Pa; (x,y) ? Lx x Ly

% INPUTS:
% dx: increment size in x
% dy: increment size in y
% Lx: length of diaphragm in x
% Ly: length of diaphragm in y
% T: tension force
% Pa: acoustic pressure

% OUTPUT:
% 1. 3D Surface Model of the deformed diaphragm
% 2. x-y Comparison of numerical and analytical solution methods*
% 3. x-y Averaged cumulative error between solution methods
% 4. 2D Contour Plot
%     *Only accurate for a sqaure mesh


Pa = 0.1; % Acoustic Pressure
T = 100; % Tension Force

Lx = 2; % Size of domain in x-direction
dx = 0.05; % Space between grid points in x-direction
Ly = 1; % Size of domain in y-direction
dy = (Lx/Ly)*dx; % Space between grid points in y-direction
Ny = Ly/dy+1; % Number of grid points in y-direction
Nx = Ny;
LxD = 2e-3;
f0 = LxD^2*Pa/T;
x = linspace(0,Lx,Ny); % x-coordinate of grid points
y = linspace(0,Ly,Ny); % y-coordinate of grid points
[X,Y] = meshgrid(x,y); % Matrix of x and y coordinates (useful for plotting)

% Initialising Matrices
f = zeros(Nx,Ny);
fm = zeros(Nx,Ny);
fd = zeros(Nx,Ny);
fr = zeros(Nx,Ny);
fx = zeros(Nx,Ny);
rel_error = zeros(Nx,Ny);
rel_errord = zeros(Nx,Ny);
rel_errorm = zeros(Nx,Ny);
rel_errorr = zeros(Nx,Ny);
rel_errorx = zeros(Nx,Ny);

% Set the deformation at the boundary
for i=1:Nx
    f(end,i) = 0;
    f(1,i) = 0;
    f(i,1) = 0;
    f(i,end) = 0;
    fm(end,i) = 0;
    fm(1,i) = 0;
    fm(i,1) = 0;
    fm(i,end) = 0;
    fd(end,i) = 0;
    fd(1,i) = 0;
    fd(i,1) = 0;
    fd(i,end) = 0;
    fr(end,i) = 0;
    fr(1,i) = 0;
    fr(i,1) = 0;
    fr(i,end) = 0;
    fx(end,i) = 0;
    fx(1,i) = 0;
    fx(i,1) = 0;
    fx(i,end) = 0;
end
% Set the deformation at the corners
% Dimensionless Solution
f(1,1) = 0;
f(end,end) = 0;
f(end,1) = 0;
f(1,end) = 0;
fold = f;
% Manufactured Solution
fm(1,1) = 0;
fm(end,end) = 0;
fm(end,1) = 0;
fm(1,end) = 0;
fmold = fm;
% Dimensional Solution
fd(1,1) = 0;
fd(end,end) = 0;
fd(end,1) = 0;
fd(1,end) = 0;
fdold = fd;
% Over-Relaxation Solution
fr(1,1) = 0;
fr(end,end) = 0;
fr(end,1) = 0;
fr(1,end) = 0;
frold = fr;
% Non-Square Solution
fx(1,1) = 0;
fx(end,end) = 0;
fx(end,1) = 0;
fx(1,end) = 0;
fxold = fx;

loop = 1;
iter = 1;
epsi = 1e-6; % Stopping criterion for Liebmann Method
lambda = 1.8;
tic
while loop
    for j=2:(Ny-1)
        for i=2:(Nx-1)
            Q(j,i) = 1-50*pi^2*sin(pi*x(i))*sin(2*pi*y(j)); % Extra Term (Q5)
            f(j,i) = 0.25*(f(j,i+1)+f(j,i-1)+f(j+1,i)+f(j-1,i)+dx^2); % Dimensionless (Q4)
            fm(j,i) = 0.25*((-dx^2*Q(j,i))+dx^2+fm(j,i+1)+fm(j,i-1)+fm(j+1,i)+fm(j-1,i)); % Manufactured (Q5)
            fr(j,i) = lambda*f(j,i)+(1-lambda)*frold(j,i); % Relaxation Term (Q6)
            fd(j,i) = 0.25*(((Pa/T)*dx^2)+fd(j,i+1)+fd(j,i-1)+fd(j+1,i)+fd(j-1,i)); % Dimensional Solution (Q7)
            fx(j,i) = (((fx(j,i+1)+fx(j,i-1))/dx^2)+((fx(j+1,i)+fx(j-1,i))/dy^2)+1)/(2*((1/(dx^2))+(1/(dy^2)))); % Non-Sqaure Solution (Q8)
            % Relative Error for each system
            rel_error(j,i) = 100*(f(j,i)-fold(j,i))/f(j,i);
            rel_errord(j,i) = 100*(fd(j,i)-fdold(j,i))/fd(j,i);
            rel_errorm(j,i) = 100*(fm(j,i)-fmold(j,i))/fm(j,i);
            rel_errorr(j,i) = 100*(fr(j,i)-frold(j,i))/fr(j,i);
            rel_errorx(j,i) = 100*(fx(j,i)-fxold(j,i))/fx(j,i);
        end
    end
    % Norm Error for each system to determine convergence
    norm_error(iter) = norm(rel_error,2);
    norm_errord(iter) = norm(rel_errord,2);
    norm_errorm(iter) = norm(rel_errorm,2);
    norm_errorr(iter) = norm(rel_errorr,2);
    norm_errorx(iter) = norm(rel_errorx,2);
    if (norm_errorx(iter)<epsi)||(iter>5000)
        loop = 0;
    else
        iter = iter+1;
        fold = f;
        fmold = fm;
        fdold = fd;
        frold = fr;
        fxold = fx;
    end
end
% Making fx dimensional
fxx = fx * f0;
toc
% Plotting the 3D deformation surface
disp('Number of iteration to convergence = ');disp(iter)
figure(1)
surf(X,Y,f)
xlabel('x')
ylabel('y')
zlabel('f(x,y)')

% Question 5: Determining Values for comparison
disp(y(floor(Nx/2)+1)); % Checking Correct Y (0.5) was used
man = (fm(:,floor(Nx/2)+1)); % Manufactured values
fman = @(x,y) (10.*sin(pi.*x).*sin(2.*pi.*y))+((x-1).*x); % Defining original eqn
xman = linspace(0,1,1000); % Making sample x array
yman = linspace(0,1,1000); % Making sample y array
xmanerr = linspace(0,1,Nx);
ymanerr = linspace(0,1,Ny);
manarray = fman(0.5,yman); % Finding analytical value to compare
manerr = fman(0.5,ymanerr);


% Plotting the comparison between manufactured and original fm
figure(2)
hold on
plot(x,man);
plot(xman,manarray);
hold off
grid minor
xlabel('x [mm]')
ylabel('Magnitude')
ylim([-15 15])
legend('Numerically Solved','Analytically Solved')

% Plotting error function between analytical and manufactured
toterror = 0;
for q = 2:Nx
    error(q) = abs(man(q)-manerr(q));
    toterror = toterror + error(q);
    errorArray(q) = toterror;
end
figure(3)
hold on
plot(xmanerr,errorArray/Nx);
grid minor
hold off
xlabel('x [mm]')
ylim([0 0.4])
ylabel('Averaged Cumulative Error')

% Plotting contours
figure(4)
[C,h] = contour(X,Y,fxx,10);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
xlabel('x [mm]')
ylabel('y [mm]')

function numerical_methods_assignment
% NUMERICAL_METHODS_ASSIGNMENT
% Numerical Methods Assignment using MATLAB functions.
% Methods: Newton-Raphson, Secant, and Runge-Kutta 4th Order.
%
% The main function calls local functions. Each method is implemented
% in both FOR-loop and WHILE-loop versions.

clc;
clear;

fprintf('=============================================\n');
fprintf('       NUMERICAL METHODS ASSIGNMENT\n');
fprintf('=============================================\n\n');

fprintf('1. NEWTON-RAPHSON METHOD - FOR LOOP\n');
root1 = newtonRaphsonFor(1.5, 0.0001, 20);
fprintf('Root = %.6f\n\n', root1);

fprintf('2. NEWTON-RAPHSON METHOD - WHILE LOOP\n');
root2 = newtonRaphsonWhile(1.5, 0.0001, 20);
fprintf('Root = %.6f\n\n', root2);

fprintf('3. SECANT METHOD - FOR LOOP\n');
root3 = secantFor(1, 2, 0.0001, 20);
fprintf('Root = %.6f\n\n', root3);

fprintf('4. SECANT METHOD - WHILE LOOP\n');
root4 = secantWhile(1, 2, 0.0001, 20);
fprintf('Root = %.6f\n\n', root4);

fprintf('5. RUNGE-KUTTA 4TH ORDER - FOR LOOP\n');
[yFor, xFor] = RK4For(0, 1, 0.1, 0.2);
fprintf('y(%.1f) = %.6f\n\n', xFor, yFor);

fprintf('6. RUNGE-KUTTA 4TH ORDER - WHILE LOOP\n');
[yWhile, xWhile] = RK4While(0, 1, 0.1, 0.2);
fprintf('y(%.1f) = %.6f\n', xWhile, yWhile);

end


% =========================================================
% NEWTON-RAPHSON: FOR LOOP
% =========================================================
function root = newtonRaphsonFor(x, tol, N)

f = @(x) x^3 - x - 2;
df = @(x) 3*x^2 - 1;

for i = 1:N
    xnew = x - f(x)/df(x);
    error = abs(xnew - x);

    fprintf('Iteration %d: x = %.6f, Error = %.6f\n', ...
        i, xnew, error);

    if error < tol
        x = xnew;
        break;
    end

    x = xnew;
end

root = x;

end


% =========================================================
% NEWTON-RAPHSON: WHILE LOOP
% =========================================================
function root = newtonRaphsonWhile(x, tol, maxIter)

f = @(x) x^3 - x - 2;
df = @(x) 3*x^2 - 1;

i = 1;

while i <= maxIter
    xnew = x - f(x)/df(x);
    error = abs(xnew - x);

    fprintf('Iteration %d: x = %.6f, Error = %.6f\n', ...
        i, xnew, error);

    if error < tol
        x = xnew;
        break;
    end

    x = xnew;
    i = i + 1;
end

root = x;

end


% =========================================================
% SECANT METHOD: FOR LOOP
% =========================================================
function root = secantFor(x0, x1, tol, N)

f = @(x) x^3 - x - 2;

for i = 1:N
    denominator = f(x1) - f(x0);

    if denominator == 0
        error('Secant method failed: denominator became zero.');
    end

    x2 = x1 - (f(x1)*(x1 - x0))/denominator;
    errorValue = abs(x2 - x1);

    fprintf('Iteration %d: x = %.6f, Error = %.6f\n', ...
        i, x2, errorValue);

    if errorValue < tol
        break;
    end

    x0 = x1;
    x1 = x2;
end

root = x2;

end


% =========================================================
% SECANT METHOD: WHILE LOOP
% =========================================================
function root = secantWhile(x0, x1, tol, maxIter)

f = @(x) x^3 - x - 2;

i = 1;

while i <= maxIter
    denominator = f(x1) - f(x0);

    if denominator == 0
        error('Secant method failed: denominator became zero.');
    end

    x2 = x1 - (f(x1)*(x1 - x0))/denominator;
    errorValue = abs(x2 - x1);

    fprintf('Iteration %d: x = %.6f, Error = %.6f\n', ...
        i, x2, errorValue);

    if errorValue < tol
        break;
    end

    x0 = x1;
    x1 = x2;
    i = i + 1;
end

root = x2;

end


% =========================================================
% RUNGE-KUTTA 4TH ORDER: FOR LOOP
% =========================================================
function [y, x] = RK4For(x, y, h, xend)

f = @(x,y) x + y;
N = round((xend - x)/h);

for i = 1:N
    k1 = h*f(x,y);
    k2 = h*f(x + h/2, y + k1/2);
    k3 = h*f(x + h/2, y + k2/2);
    k4 = h*f(x + h, y + k3);

    y = y + (k1 + 2*k2 + 2*k3 + k4)/6;
    x = x + h;

    fprintf(['Step %d: x = %.2f, k1 = %.6f, k2 = %.6f, ' ...
             'k3 = %.6f, k4 = %.6f, y = %.6f\n'], ...
             i, x, k1, k2, k3, k4, y);
end

end


% =========================================================
% RUNGE-KUTTA 4TH ORDER: WHILE LOOP
% =========================================================
function [y, x] = RK4While(x, y, h, xend)

f = @(x,y) x + y;
i = 1;

while x < xend - 1e-12
    k1 = h*f(x,y);
    k2 = h*f(x + h/2, y + k1/2);
    k3 = h*f(x + h/2, y + k2/2);
    k4 = h*f(x + h, y + k3);

    y = y + (k1 + 2*k2 + 2*k3 + k4)/6;
    x = x + h;

    fprintf(['Step %d: x = %.2f, k1 = %.6f, k2 = %.6f, ' ...
             'k3 = %.6f, k4 = %.6f, y = %.6f\n'], ...
             i, x, k1, k2, k3, k4, y);

    i = i + 1;
end

end

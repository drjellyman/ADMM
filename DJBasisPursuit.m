clear all;
close all; 

m = 10;
n = 30;
A = randn(m,n);
xOrig = randn(n,1);
b = A*xOrig;
rho = 1;
toleranceAbs = 1e-5;
toleranceRel = 5e-3;

x = zeros(n,1);
z = zeros(n,1); 
u = zeros(n,1);

% Run the loop
for k = 1:1000
    x = (eye(n)-A'*inv(A*A')*A)*(z-u)-A'*inv(A*A')*b;
    zold = z;
    z = max(0,x+u-(1/rho)) - max(0,-(x+u)-(1/rho));    
    u = u + x - z; 
    history.rNorm(k) = norm(x-z);
    history.sNorm(k) = norm(-rho*(z-zold));
    history.epsPri(k) = sqrt(n)*toleranceAbs+toleranceRel*max(norm(x), norm(-z));
    history.epsDual(k) = sqrt(n)*toleranceAbs+toleranceRel*norm(rho*u);
    
    if (history.rNorm(k) < history.epsPri(k) && history.sNorm(k) ...
            < history.epsDual(k))
        break        
    end
end
IterDim = [1:length(history.rNorm)];
figure; semilogy(IterDim, history.rNorm, 'b', IterDim, history.sNorm, 'r', ...
    IterDim, history.epsPri, '--b', IterDim, history.epsDual, '--r'); 
legend('rNorm', 'sNorm', 'epsPri', 'epsDual'); grid on; 
xlabel('No of iterations'); 


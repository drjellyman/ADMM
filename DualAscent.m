clear all; 
close all; 

x = [-3:0.01:3];
f_x = x.^2;
A = -2; 
b = -3; 
c_x = A*x-b;
lambda = 0.1; 
rho = 0.1;

figure; hold on; plot(x, f_x, x, c_x);
for k = 1:50
    % Calculate x_k+1 using previous lambda
    L = f_x+lambda*(A*x-b);
    [g, i] = min(L);
    % Calculate lambda update using x_k+1
    lambda = lambda + rho*(A*x(i)-b);
    % Store progress
    % plot(x(i), f_x(i), 'go');
    history.objval(:,k) = [x(i), f_x(i)]'; 
    history.lambda(k) = lambda;
    history.g(k) = g;
    history.L(:,k) = L; 
end

plot_txt_size = 18; 
plot(history.objval(1,:), history.objval(2,:), 'go'); 
legend('f(x)','Ax=b','min(L(x,lambda))'); set(gca,'fontsize',plot_txt_size); 
xlim([-3, 3]); xlabel('x'); grid on;
figure; plot(history.objval(2,:)); legend('Objective value'); set(gca,'fontsize',plot_txt_size); xlabel('No. of iterations'); grid on;
figure; plot(history.lambda); legend('Lambda'); set(gca,'fontsize',plot_txt_size); xlabel('No. of iterations'); grid on;
figure; plot(history.g); legend('g'); set(gca,'fontsize',plot_txt_size); xlabel('No. of iterations'); grid on;
figure; plot(history.L); legend('L'); set(gca,'fontsize',plot_txt_size); xlabel('x'); grid on; xlim([0, 600]);

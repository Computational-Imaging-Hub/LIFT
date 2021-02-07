function [x,convergence] = Solver_PlugPlay_FISTA(A,At,b,x0,opt)
% The objective function: F(x) = 1/2 ||y - hx||^2 + lambda |x|
% Input: 
%   A: the forward operator function, input is a 3d cube, output a 2d
%   matrix
%   At: backward operator function, input is a 2D image, output a 3d matrix
%   b: degraded image (2D matrix)
%   x0: initialization (3D matrix)
%   opt.lambda: weight constant for the regularization term
% Output:
%   x: output image stack (3D)
%
% Author: Xiaohua Feng

fprintf(' - Running FISTA with fixed step size method\n');

lambda = opt.lambda;
maxiter = opt.maxiter;
tol = opt.tol;
vis = opt.vis;
xk = x0;
yk = xk;
zk = yk;
tk = 1;
L0 = opt.step;
POSCOND= opt.POScond;
% k-th (k=0) function, gradient, hessian
objk  = func(x0,b,A,lambda);
fprintf('%6s %9s %9s\n','iter','f','sparsity');
fprintf('%6i %9.2e %9.2e\n',0,objk,nnz(xk)/numel(xk));
convergence = zeros(maxiter,1,'gpuArray');
for i = 1:maxiter
%     i
    x_old = xk;
    y_old = yk;
    t_old = tk;

    yg = y_old - 1/L0*(At(A(y_old)-b));
    switch (opt.denoiser)        
        case 'ProxTV'
            denoise = @denoise_ProxTV;
        case 'Prox_l1'
            denoise = @denoise_l1;
        case 'Prox_l1_GPU'
            denoise = @denoise_l1_GPU;
        otherwise
            denoise = @denoise_l1;
    end
    zk = denoise(yg,lambda/L0);
%     fx = func(x_old,b,A,lambda);
    xk = zk;
    
    if(POSCOND)
        xk(xk<0)=0; % positiveness constraint
    end
    tk = (1/20+sqrt(1/2+4*t_old*t_old))/2;
    yk = xk + (t_old-1)/tk*(xk-x_old)+(t_old)/tk*(zk-xk);
    if(POSCOND)
        yk(yk<0)=0; % positiveness constraint
    end
    if vis > 0
        fprintf('%6i %9.2e %9.2e\n',i,func(xk,b,A,lambda),nnz(xk)/numel(xk));
    end
    
%     convergence(i) = fx;

end
x = xk;
end

function norm_val = normNdMatrix(x,n)
    norm_val_temp = x.^n;
    norm_val = sum(norm_val_temp(:));
end

function Fx = func(xk,b,A,lambda)
    e = b - A(xk);
    % Fx = 0.5*normNdMatrix(e,2) + lambda*norm_tv3d(xk);
    Fx = 0.5*normNdMatrix(e,2);
end


function img = denoise_ProxTV(noisy,sigma_hat)
	img = prox_tv(noisy,sigma_hat);
end

function img_estimated = denoise_l1(noisy,sigma_hat)
    img_estimated = prox_l1(noisy, sigma_hat);
end

function img_estimated = denoise_l1_GPU(noisy,sigma_hat)
    img_estimated = prox_l1_GPU(noisy, sigma_hat);
end
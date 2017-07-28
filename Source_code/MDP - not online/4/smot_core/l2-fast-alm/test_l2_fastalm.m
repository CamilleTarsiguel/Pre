clear;
u = [1:1:8;2:3:23];

[D N] = size(u);

up = u;
Omega = ones(1,N);
up(:,4:5) = 0;
Omega(4:5) = 0;

lambda = 10.0;

u_hat_alm = l2_fastalm_mo(up,lambda,'omega',Omega);
% u_hat_alm = srpca_alm(u,lambda);
% % u_hat_alm = srpca_alm_not_working(u,lambda);

P = sparse(N*D,N*D);
P(1:N*D+1:N^2*D^2) = reshape(repmat(Omega,[D 1]),[1 N*D]);

cvx_begin sdp
    variable u_hat_cvx(size(up))    
    minimize norm_nuc(cvx_hankel_mo(u_hat_cvx)) + lambda/2 * pow_pos(norm(P*u_hat_cvx(:) - P*up(:),'fro'),2)
cvx_end

[u; u_hat_alm; u_hat_cvx]



%%

u = [1:1:8;2:3:23];

[D N] = size(u);

up = u;
Omega = ones(1,N);
up(:,4:5) = 0;
Omega(4:5) = 0;

tic
for i=1:100
    u_hat1 = hstln_mo(up,3,'omega',Omega);
end
toc;

profile on;
tic
for i=1:100
    u_hat2 = fast_hstln_mo(up,3,'omega',Omega);
end
toc;
profile off;
profile viewer;
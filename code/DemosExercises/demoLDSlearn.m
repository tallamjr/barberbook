function demoLDSlearn
% DEMOLDSLEARN learning a linear dynamical system using EM
close all
% generate a latent sine wave:
T=100;
lat=sin(0.2*(1:T)); B=randn(2,1);
v=B*lat; v=v+0.1*randn(size(v));

% learn a two-dimensional latent LDS to represent the data
H=2;
opts.init.CovV=ones(2,1); % diagonal covV

[A B CovH CovV meanP CovP meanV meanH loglik]=LDSlearn(v,H,opts);
[f,F,g,G,Gp,loglik]=LDSsmooth(v,A,B,CovH,CovV,CovP,meanP,meanH,meanV);

vv=B*g; % mean reconstruction of output

subplot(1,2,1); plot(v'); title('original data');
subplot(1,2,2); plot(vv'); title('mean LDS reconstruction')
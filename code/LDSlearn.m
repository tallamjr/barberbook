function [A B CovH CovV meanP CovP meanV meanH loglik]=LDSlearn(v,H,varargin)
%LDSLEARN Learn parameters for a Latent Linear Dynamical System (Kalman Filter)
% [A B CovH CovV priormean CovP meanV meanH loglik]=LDSlearn(v,H,<opts>)
%
% The method is based on maximum likelihood and the Expectation Maximisation algorithm
%
% Inputs:
% v : observation matrix of dimension V x T. Each column contains a datapoint, v(:,1:T)
% H : hidden dimension of the latent LDS
% opts.maxits : maximum number of iterations of the EM algorithm
% opts.plotprogress : set to 1 to plot the log likelihood evolution
% opts.tol : log likelihood change termination criterion
%
% opts.init.A
% opts.init.B
% opts.initv.CovV
% opts.init.CovH
% opts.init.CovP
% opts.init.meanP
% opts.init.meanH
% opts.init.meanV
%
% If either A, CovV, CovH, CovP, CovV are initialised to a (non-zero) vector, the the solution learned is diagonal.
% Eg if we use opts.init.A=ones(2,1); the learned A will be diagonal
% If either of meanP, meanH, or MeanV is initialised to the zero
% vector, then they will remain the zero vector and not be updated.
%
% Outputs:
% A : H x H transition matrix
% B : V x H emission matrix
% CovH : H x H transition covariance
% CovV : V x V emission covariance
% meanH : H x 1 transition mean
% meanV : V x 1 emission mean
% covP : H x H initial prior
% meanP : H x 1 initial mean
% loglik: log likelihod of the sequence log p(v)
%
% example see demoLDSlearn.m

[V T]=size(v); 
init.A=eye(H);
[usvd ssvd vsvd]=svd(v);
init.B=usvd(:,1:H);
init.CovH=eye(H);
init.CovP=eye(H);
vv=usvd(:,1:H)*ssvd(1:H,1:H)*vsvd(1:H,:);
init.CovV=cov((v-vv)');
init.meanP=zeros(H,1);
init.meanH=zeros(H,1);
init.meanV=zeros(V,1);

opts=[];if nargin==3; opts=varargin{1}; end
opts=setfields(opts,'tol',1e-5,'maxits',100,'plotprogress',1,'init',init); %    default options

loglikold=-realmax;

A=opts.init.A;
B=opts.init.B;
CovV=opts.init.CovV;
CovH=opts.init.CovH;
CovP=opts.init.CovP;
meanH=opts.init.meanH;
meanV=opts.init.meanV;
meanP=opts.init.meanP;

diagCovH=0; diagCovV=0; diagCovP=0; diagA=0;
if isvector(CovH); diagCovH=1; end
if isvector(CovV); diagCovV=1; end
if isvector(CovP); diagCovP=1; end
if isvector(A); diagA=1; end

for loop=1:opts.maxits
    
    [f,F,g,G,Gp,lik]=LDSsmooth(v,A,B,CovH,CovV,CovP,meanP,meanH,meanV,1);
    loglik(loop)=lik;
    if opts.plotprogress; plot(loglik,'-o');  title('log likelihood'); drawnow; end
    if loop>1 & loglik(loop)<loglik(loop-1); warning('log likelihood decreased'); end
    if loglik(loop) - loglikold < opts.tol; break; end
    
    if ~all(meanP==0); meanP=g{1}; end
    CovP=G{1}+g{1}*g{1}'-g{1}*meanP'-meanP*g{1}'+meanP*meanP';
    CovP=(CovP+CovP')/2;
    if diagCovP; CovP=diag(diag(CovP)); end
    
    HH=zeros(H);
    HHp=zeros(H);
    VH=zeros(V,H);
    for t=1:T-1
        HH=HH+G{t}+g{t}*g{t}';
        HHp=HHp+Gp{t}'-meanH*g{t}';
        VH=VH+(v(:,t)-meanV)*g{t}';
    end
    A=HHp/HH;
    if diagA; A=diag(diag(A)); end
    HH=HH+G{T}+g{T}*g{T}';
    VH=VH+(v(:,T)-meanV)*g{T}';
    B=VH/HH;
    CovVtmp=zeros(V);
    meanVtmp=zeros(V,1);
    for t=1:T
        dv=v(:,t)-meanV;
        CovVtmp=CovVtmp+dv*dv'-dv*g{t}'*B'-B*g{t}*dv'+B*(G{t}+g{t}*g{t}')*B';
        meanVtmp=meanVtmp+v(:,t)-B*g{t};
    end
    if ~all(meanV==0); meanV=meanVtmp/T; end
    CovV=CovVtmp/T;CovV=(CovV+CovV')/2;
    if diagCovV; CovV=diag(diag(CovV)); end
    
    CovHtmp=zeros(H); meanHtmp=zeros(H,1);
    for t=1:T-1
        CovHtmp=CovHtmp+G{t+1}+g{t+1}*g{t+1}'-A*Gp{t}-Gp{t}'*A'+A*(G{t}+g{t}*g{t}')*A';
        CovHtmp=CovHtmp-g{t+1}*meanH'-meanH*g{t+1}'+meanH*meanH'+meanH*g{t}'*A'+A*g{t}*meanH';
        meanHtmp=meanHtmp+g{t+1}-A*g{t};
    end
    CovH=CovHtmp/(T-1); CovH=(CovH+CovH')/2;
    if diagCovH; CovH=diag(CovH); end
    
    if ~all(meanH==0); meanH=meanHtmp/(T-1); end
    loglikold=loglik;
end
loglik=lik;
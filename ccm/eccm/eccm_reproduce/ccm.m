function [SugiCorr, SugiX, origX] = ccm(X, Y, tau, E, tp)
%%
switch nargin
    case 5
    case 4
        tp = 0;
    otherwise
        error('Bad input')
end
%%
X=reshape(X,[],1);
Y=reshape(Y,[],1);
if tp >= 0
    X=X(1:end-tp);
    Y=Y(1+tp:end);
else
    X=X(1-tp:end);
    Y=Y(1:end+tp);
end
%%
L=length(X);
T=1+(E-1)*tau;
Xm=zeros((L-T+1),E);
Ym=zeros((L-T+1),E);
SugiN=E+1;
N = L-T+1;
%% RECONTRUCTIONS OF ORIGINAL SYSTEMS

for t=1:(L-T+1)
    Xm(t,:)=X((T+t-1):-tau:(T+t-1-(E-1)*tau));
    Ym(t,:)=Y((T+t-1):-tau:(T+t-1-(E-1)*tau));
end
%%
SugiX=zeros(N,1);
SugiY=zeros(N,1);

%% CMM

dat=floor((L-T+1)/2);

for ii=(dat+1):(L-T+1)
    [n1s,d1s]=knnsearch(Xm((ii-dat):(ii-1),:),Xm(ii,:),'k',SugiN);
    [n2s,d2s]=knnsearch(Ym((ii-dat):(ii-1),:),Ym(ii,:),'k',SugiN);
    
    u1s=exp(-d1s/d1s(1));
    w1s=u1s/sum(u1s);
    SugiY(ii)= w1s*Y(n1s+T-1+ii-(dat+1));
    
    u2s=exp(-d2s/d2s(1));
    w2s=u2s/sum(u2s);
    SugiX(ii)= w2s*X(n2s+T-1+ii-(dat+1));
end

origY=Y(T:end);
origY=origY((dat+1):(L-T+1));
SugiY=SugiY((dat+1):(L-T+1));
origX=X(T:end);
origX=origX((dat+1):(L-T+1));
SugiX=SugiX((dat+1):(L-T+1));

%%

SugiCorr1=corrcoef(origY,SugiY,'Rows','complete');
SugiCorr(2,1)=SugiCorr1(1,2);

SugiCorr2=corrcoef(origX,SugiX,'Rows','complete');
SugiCorr(1,1)=SugiCorr2(1,2);

% SugiR(2,1)=sqrt((sum((origY-SugiY).^2)/numel(origY)))/std(origY);
% SugiR(1,1)=sqrt((sum((origX-SugiX).^2)/numel(origX)))/std(origX);

end
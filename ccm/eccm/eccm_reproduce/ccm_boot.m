function [SugiCorr, SugiX, origX, origY, SugiY] = ccm_boot(X, Y, tau, E, tp, num_samples)

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
SugiX=zeros(N,num_samples);
SugiY=zeros(N,num_samples);
%%
origY=Y(T:end);
origX=X(T:end);
%% boot
parfor ii=1:N
    for jj = 1:num_samples
        r = randsample(N, N, true);

        [n1s,d1s]=knnsearch(Xm(r,:),Xm(ii,:),'k',SugiN+1);
        [n2s,d2s]=knnsearch(Ym(r,:),Ym(ii,:),'k',SugiN+1);
        
        n1s = r(n1s);
        n2s = r(n2s);
        
        d1s = d1s(:,2:SugiN+1);
        d2s = d2s(:,2:SugiN+1);
        
        u1s=exp(-d1s/d1s(1));
        w1s=u1s/sum(u1s);
        SugiY(ii,jj)= w1s*origY(n1s(2:end));
    
        u2s=exp(-d2s/d2s(1));
        w2s=u2s/sum(u2s);
        SugiX(ii,jj)= w2s*origX(n2s(2:end));
    end
end

%%
SugiCorr = zeros(2, num_samples);
for zz = 1:num_samples
    SugiCorr1=corrcoef(origY,SugiY(:,zz),'Rows','complete');
    SugiCorr(2,zz)=SugiCorr1(1,2);

    SugiCorr2=corrcoef(origX,SugiX(:,zz),'Rows','complete');
    SugiCorr(1,zz)=SugiCorr2(1,2);
end
% SugiR(2,1)=sqrt((sum((origY-SugiY).^2)/numel(origY)))/std(origY);
% SugiR(1,1)=sqrt((sum((origX-SugiX).^2)/numel(origX)))/std(origX);

end



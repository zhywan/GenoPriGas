%Pearson's chi-squared test for a matrix
%input a matrix with each column as an attribute, each row as an
%observation. Values are 0,1,2,...,k-1.

function [r2,p]=chi2p(X,k)

C=zeros(k,k);
E=C;
C_rs=zeros(k,1);
C_cs=zeros(1,k);
C_sum=0;
[n,m]=size(X);
XX=zeros(n,m,k);
S=zeros(k,k,n);%contingency table
r2=eye(m);
p=zeros(m,m);
for i=1:k
    XX(:,:,i)=X==(i-1);
end
for i=1:m
    i
    for j=1:i
        for ii=1:k
            for jj=1:k
                S(ii,jj,:)=XX(:,i,ii)&XX(:,j,jj);
                C(ii,jj)=sum(S(ii,jj,:));
            end
        end
        C_rs=sum(C,2);
        C_cs=sum(C);
        C_sum=sum(C_rs);
        E=C_rs*C_cs/C_sum;
        Z_rs=sum(C_rs==0);
        Z_cs=sum(C_cs==0);
        if C_cs(2)+C_cs(3)+C_rs(2)+C_rs(3)==0 % to prevent NaN (compare all zero)
            r2(i,j)=n;
            p(i,j)=0;
        elseif C_cs(2)+C_cs(3)==0 ||C_rs(2)+C_rs(3)==0% prevent NaN (nonzero & zero)
            r2(i,j)=0;
            p(i,j)=1;
        else
            if (Z_rs+Z_cs)==0
                r2(i,j)=sum(sum((C-E).^2./E));
                p(i,j)=1-chi2cdf(r2(i,j),(k-1)^2);
            else
                E_valid=E~=0;
                r2(i,j)=sum(sum((C(E_valid)-E(E_valid)).^2./E(E_valid)));
                freedom=(k-Z_rs-1)*(k-Z_cs-1);
                p(i,j)=1-chi2cdf(r2(i,j),freedom);
            end
        end
    end
end
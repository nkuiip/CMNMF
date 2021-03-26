function [L,W,H1,H2] = CMNMF_LF( MaxIter,normalization,A1,A2,W,H1,H2,M,alpha,beta,lambda1,lambda2)
%L为目标函数值
D1 = diag(sum(M,2));
D2 = diag(sum(M,1));
%L=zeros(1+2*Inner_MaxIter,MaxIter);
L=zeros(MaxIter*3,16);
[L(1,1),L(1,5),L(1,9),L(1,13)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
t=1;
for i=1:MaxIter
    disp(['Iter= ' num2str(i)]);
    %update W
    W=W.*(A1*H1'+alpha*A2*H2')./(W*(H1*H1')+alpha*W*(H2*H2')+lambda1*W);
    W(isnan(W))=0;
    W(isinf(W))=0;
    t = t+1;
    [L(t,1),L(t,5),L(t,9),L(t,13)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
    L(t,2)=L(t-1,1)-L(t,1);
    L(t,3)=L(t,2)/L(t-1,1);
    L(t,4)=L(t,1)/L(2,1);
    L(t,6)=L(t-1,5)-L(t,5);
    L(t,7)=L(t,6)-L(t-1,5);
    L(t,8)=L(t,5)/L(2,5);
    L(t,10)=L(t-1,9)-L(t,9);
    L(t,11)=L(t,10)/L(t-1,9);
    L(t,12)=L(t,9)/L(2,9);
    L(t,14)=L(t-1,13)-L(t,13);
    L(t,15)=L(t,14)/L(t-1,13);
    L(t,16)=L(t,13)/L(2,13);
    %update H1
    H1=H1.*((W'*A1 + beta*H2*M')./(W'*W*H1 + beta*H1*D1));
    H1(isnan(H1))=0;
    H1(isinf(H1))=0;
    t = t+1;
    [L(t,1),L(t,5),L(t,9),L(t,13)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
    L(t,2)=L(t-1,1)-L(t,1);
    L(t,3)=L(t,2)/L(t-1,1);
    L(t,4)=L(t,1)/L(2,1);
    L(t,6)=L(t-1,5)-L(t,5);
    L(t,7)=L(t,6)-L(t-1,5);
    L(t,8)=L(t,5)/L(2,5);
    L(t,10)=L(t-1,9)-L(t,9);
    L(t,11)=L(t,10)/L(t-1,9);
    L(t,12)=L(t,9)/L(2,9);
    L(t,14)=L(t-1,13)-L(t,13);
    L(t,15)=L(t,14)/L(t-1,13);
    L(t,16)=L(t,13)/L(2,13);
    %update H2
    H2=H2.*((alpha*W'*A2 + beta*H1*M)./(alpha*(W'*W*H2) + beta*H2*D2));
    H2(isnan(H2))=0;
    H2(isinf(H2))=0;
    t = t+1;
    [L(t,1),L(t,5),L(t,9),L(t,13)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
    L(t,2)=L(t-1,1)-L(t,1);
    L(t,3)=L(t,2)/L(t-1,1);
    L(t,4)=L(t,1)/L(2,1);
    L(t,6)=L(t-1,5)-L(t,5);
    L(t,7)=L(t,6)-L(t-1,5);
    L(t,8)=L(t,5)/L(2,5);
    L(t,10)=L(t-1,9)-L(t,9);
    L(t,11)=L(t,10)/L(t-1,9);
    L(t,12)=L(t,9)/L(2,9);
    L(t,14)=L(t-1,13)-L(t,13);
    L(t,15)=L(t,14)/L(t-1,13);
    L(t,16)=L(t,13)/L(2,13);
    disp(['L=' num2str(L(t,1))]);
    
    if normalization==5
        D = diag(1./sqrt(sum(W.^2,1)));
        D_inverse = diag(sqrt(sum(W.^2,1)));
        W = W*D;
        H1 = D_inverse * H1;
        H2 = D_inverse * H2;
    end
    if normalization==6
        D = diag(1./sqrt(sum(W.^2,2)));
        D_inverse=pinv(D);
        W = D*W;
        H1 = H1*D_inverse;
        H2 = H2*D_inverse;
    end
    if normalization==7
        D = diag(1./sum(W,1));
        D_inverse = diag(sum(W,1));
        W = W*D;
        H1 = D_inverse * H1;
        H2 = D_inverse * H2;
    end
    if normalization==8
        D = diag(1./sum(W,2));
        D_inverse=pinv(D);
        W = D*W;
        H1 = H1*D_inverse;
        H2 = H2*D_inverse;
    end
    
    if normalization>=5&&normalization<=8
        t=t+1;
        [L(t,1),L(t,5),L(t,9),L(t,13)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
        L(t,2)=L(t-1,1)-L(t,1);
        L(t,3)=L(t,2)/L(t-1,1);
        L(t,4)=L(t,1)/L(2,1);
        L(t,6)=L(t-1,5)-L(t,5);
        L(t,7)=L(t,6)-L(t-1,5);
        L(t,8)=L(t,5)/L(2,5);
        L(t,10)=L(t-1,9)-L(t,9);
        L(t,11)=L(t,10)/L(t-1,9);
        L(t,12)=L(t,9)/L(2,9);
        L(t,14)=L(t-1,13)-L(t,13);
        L(t,15)=L(t,14)/L(t-1,13);
        L(t,16)=L(t,13)/L(2,13);
        disp(['After normalization ,L=' num2str(L(t,1))]);
    end
end

if normalization==1
    D = diag(1./sqrt(sum(W.^2,1)));
    D_inverse = diag(sqrt(sum(W.^2,1)));
    W = W*D;
    H1 = D_inverse * H1;
    H2 = D_inverse * H2;
end
if normalization==2
    D = diag(1./sqrt(sum(W.^2,2)));
    W = D*W;
    D_inverse=pinv(D);
    H1 = H1*D_inverse;
    H2 = H2*D_inverse;
end
if normalization==3
    D = diag(1./sum(W,1));
    D_inverse = diag(sum(W,1));
    W = W*D;
    H1 = D_inverse * H1;
    H2 = D_inverse * H2;
end
if normalization==4
    D = diag(1./sum(W,2));
    D_inverse=pinv(D);
    W = D*W;
    H1 = H1*D_inverse;
    H2 = H2*D_inverse;
end

if normalization>=1&&normalization<=4
    t=t+1;
    [L(t,1),L(t,5),L(t,9),L(t,13)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
    L(t,2)=L(t-1,1)-L(t,1);
    L(t,3)=L(t,2)/L(t-1,1);
    L(t,4)=L(t,1)/L(2,1);
    L(t,6)=L(t-1,5)-L(t,5);
    L(t,7)=L(t,6)-L(t-1,5);
    L(t,8)=L(t,5)/L(2,5);
    L(t,10)=L(t-1,9)-L(t,9);
    L(t,11)=L(t,10)/L(t-1,9);
    L(t,12)=L(t,9)/L(2,9);
    L(t,14)=L(t-1,13)-L(t,13);
    L(t,15)=L(t,14)/L(t-1,13);
    L(t,16)=L(t,13)/L(2,13);
    disp(['After normalization ,L=' num2str(L(t,1))]);
end

%     W=diag(1./sum(W,2))*W;
%     W(isnan(W))=0;
%     H1=H1*diag(1./sum(H1));
%     H2=H2*diag(1./sum(H2));
end

function [O,O1,O2,O3]= objective(alpha,beta,lambda1,lambda2,V1,V2,W,H1,H2,M)

O1=sum(sum((V1-W*H1).^2));
O2 = sum(sum((V2-W*H2).^2));
%O2=-gama*trace(H1*M*H2')+lamat1*sum(sum((W).^2))+lamat2*(sum(sum(H1.^2))+sum(sum(H2.^2)));
D1 = diag(sum(M,2));
D2 = diag(sum(M,1));
O3=(trace(H1*D1*H1') + trace(H2*D2*H2')-2*trace(H1*M*H2'));
O=O1+alpha*O2+beta*O3;

end



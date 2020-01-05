function TSVD()
    %do quad to get wk, tk
    %do logarithmically distributed to get sj
    N = 40;
    S = log_dis(N);
    Y = getY(S);
    [T W] = legpts(N,[0,5],'FAST');
    %use s t to get A
    A = getA(W,S,T);
    [U,S,V] = svd(A);
    %singuNum = 40;
    tol = 1e-6
    
    %Utrun = U(:,1:singuNum);
    Strun = zeros(N,N);
    for i=1:N 
        if(S(i,i)<tol)
            Strun(i,i) = 0;
        else 
            Strun(i,i) = 1./S(i,i);
        end
    end
    Strun
    Xcal = V*Strun*U'*Y;
    %Xcal = pinv(A)*Y
    Xtrue = getTrueX(T);
    diag(Strun)
    
    norm(Xcal - Xtrue,2)
end

function Xtrue = getTrueX(T)
    N = size(T);
    Xtrue = zeros(N);
    for i = 1:N 
        t = T(i);
        if(t<=1)
            Xtrue(i) = t;
        elseif(1<=t && t<3)
            Xtrue(i) = 3/2-t/2;
        elseif(3<=t)
            Xtrue(i) = 0;
        end
    end
end
function A = getA(W,S,T)
    J = size(S,1);
    K = size(T,1);
    A = zeros(J,K);
    for j = 1:J
        for k = 1:K
            A(j,k) = W(k)*exp((-1)*S(j)*T(k));
        end
    end
end

function Y = getY(S)
    N = size(S);
    Y = zeros(N);
    for i = 1:N 
        Y(i) = getLf(S(i));
    end
end

function S = log_dis(N)
    S = zeros(N,1);
    for j = 1:N 
        temp = (-1 + (j-1)/20)*log(10);
        S(j) = exp(temp);
    end
end

function Lf = getLf(s)
    Lf = (2-3*exp((-1)*s)+exp((-3)*s))/(2*(s^2));
end
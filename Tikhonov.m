function Tik()
    N = 40;
    S = log_dis(N);
    Y = getY(S);
    %add random noise
    rad = 1e-4
    noise = (-1)*rad + 2*rad*rand(N,1)
    Ynoise = Y+noise;

    [T W] = legpts(N,[0,5],'GW');
    %use s t to get A
    A = getA(W,S,T);
    %use normal equation to get Xcal
    %(A^TA+alpha*I)xa = A^Tb
    delta = 1e-3;
    I = eye(N);
    Aplus = (A'*A+delta*delta*I);
    B = A'*Ynoise;
    %Aplus*X = B
    Xcal = Aplus \ B;
    Xtrue = getTrueX(T);
    norm(Xcal-Xtrue,2)
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
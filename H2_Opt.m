function [Jopt,Zopt] = H2_Opt(A,B,E,C,D,G,H,nc,options,Zstruct,Jfeas,Zfeas,P)

%Options
epsilon = options(1);
itmax = options(2);

%Matrices Dimensions
nx = size(A,2);
nu = size(B,2);
ny = size(G,1);
nw = size(E,2);
nz = size(C,1); 
Nsub = size(Zstruct,1);

%Satisfying Orthogonality Condition
An = A - B/(D'*D)*D'*C;
Cn = C - D/(D'*D)*D'*C;

%Interations
ito = 0;
Jopt(ito+1) = Jfeas;

while true 
    ito = ito + 1;

    Ap = An - B/(D'*D)*B'*P;
    Cp = Cn - D/(D'*D)*B'*P;

    %LMI description
    setlmis([])
    
    S = lmivar(1,[nx 1]);
    W = lmivar(1,[nw 1]);

    if nc == 0
        Z = lmivar(2,[nu ny]);
    elseif nc > 0
        for i = 1:Nsub
            [R{i},n,sR{i}] = lmivar(2,Zstruct(i,:));
        end
        Z = lmivar(3,[sR{1},sR{2};sR{3},zeros(Zstruct(3,1),Zstruct(2,2))]);
    end

    ct = 1;
    lmiterm([-ct,1,1,S],1,1);
    
    ct = ct+1;
    lmiterm([ct,1,1,S],Ap',1,'s');
    lmiterm([ct,1,1,0],Cp'*Cp);
    lmiterm([ct,2,1,Z],1,G);
    lmiterm([ct,2,1,S],(D'*D)\B',1);
    lmiterm([ct,2,1,0],(D'*D)\D'*C);
    lmiterm([ct,2,2,0],-inv(D'*D));
    
    ct = ct+1;
    lmiterm([-ct,1,1,W],1,1);
    lmiterm([-ct,2,1,0],E);
    lmiterm([-ct,2,1,Z],B,H);
    lmiterm([-ct,2,2,0],2*inv(P));
    lmiterm([-ct,2,2,S],-inv(P),inv(P));

    lmisys = getlmis;

    %Objective Function
    options = [1e-4,200,1e9,100,1];
    np = decnbr(lmisys);
    c = zeros(np,1);
    for i=1:np
        Wi = defcx(lmisys,i,W);
        c(i) = trace(Wi);
    end

    [copt,xopt] = mincx(lmisys,c,options);

    if (isempty(copt))
        Ss = [];
        Zs = [];
        disp('not able to perfom de minimization on optimality');
        return
    end

    Sopt = dec2mat(lmisys,xopt,S);
    Zopt = dec2mat(lmisys,xopt,Z);

    L = Zopt*G;
    F = B*Zopt;
    AL = A + B*L;
    CL = C + D*L;
    EF = E + F*H;
    P = lyap(AL',CL'*CL);
    Jopt(ito+1) = trace(EF'*P*EF);
    if Jopt(ito+1)/Jopt(ito) >= 1 - epsilon
        break;  
    end
    
end
return;
end

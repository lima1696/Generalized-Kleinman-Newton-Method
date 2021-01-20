function [Jfeas,Zfeas,P,f] = H2_Feas(A,B,E,Cfeas,Dfeas,G,H,nc,options,Zstruct,C,D)

%Options
epsilon = options(1);
itmax = options(2);
tol = options(3);

%Matrices Dimensions
nx = size(A,2);
nu = size(B,2);
ny = size(G,1);

Nsub = size(Zstruct,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Feasibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solve Riccati Equation
[~,~,Yo] = care(A,B,Cfeas'*Cfeas,Dfeas'*Dfeas,Cfeas'*Dfeas);
Yo = -Yo;

%Matrices AL and CL
AL = A + B*Yo;
CL = Cfeas + Dfeas*Yo;
P = lyap(AL',CL'*CL);

%Initial Penalty
penalty = 2*sum(svd(Yo));

%Satisfying Orthogonality Condition
An = A - B/(Dfeas'*Dfeas)*Dfeas'*Cfeas;
Cn = Cfeas - Dfeas/(Dfeas'*Dfeas)*Dfeas'*Cfeas;

%Interations Number
it = 0;


f(it+1) = penalty;
while true 
    
    it = it + 1;

    An = A - B/(Dfeas'*Dfeas)*Dfeas'*Cfeas;
    Cn = Cfeas - Dfeas/(Dfeas'*Dfeas)*Dfeas'*Cfeas;

    Ap = An - B/(Dfeas'*Dfeas)*B'*P;
    Cp = Cn - Dfeas/(Dfeas'*Dfeas)*B'*P;

    %LMI descriptions
    setlmis([])

    
    S = lmivar(1,[nx 1]);
    W = lmivar(1,[nx 1]);
    X = lmivar(1,[nu 1]);
    Y = lmivar(2,[nu nx]);
    
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
    lmiterm([ct,2,1,Y],1,1);
    lmiterm([ct,2,1,Z],1,G);
    lmiterm([ct,2,1,S],(Dfeas'*Dfeas)\B',1);
    lmiterm([ct,2,1,0],(Dfeas'*Dfeas)\Dfeas'*Cfeas);
    lmiterm([ct,2,2,0],-inv(Dfeas'*Dfeas));
    
    if nc > 0
         ct = ct+1;
         lmiterm([-ct,1,1,0],tol);
         lmiterm([-ct,2,1,Z],1,1);
         lmiterm([-ct,2,2,0],tol);
    end

    ct = ct+1;
    lmiterm([-ct,1,1,W],1,1);
    lmiterm([-ct,2,1,Y],1,1);
    lmiterm([-ct,2,2,X],1,1);

    lmisys = getlmis;

    %Objective function
    options = [1e-4,200,1e9,100,1];

    np = decnbr(lmisys);
    c = zeros(np,1);
    for i=1:np
        Wi = defcx(lmisys,i,W);
        Xi = defcx(lmisys,i,X);
        c(i) = trace(Wi) + trace(Xi);
    end

    [copt,xopt] = mincx(lmisys,c,options);

    if (isempty(copt))
        Ss = [];
        Xs = [];
        Zs = [];
        Ys = [];
        Ws = [];
        disp('not able to perfom de minimization on feasiability');
        return
    end

    Sfeas = dec2mat(lmisys,xopt,S);
    Yfeas = dec2mat(lmisys,xopt,Y);
    Zfeas = dec2mat(lmisys,xopt,Z);


    L = Zfeas*G + Yfeas ;
    AL = A+B*L;
    CL = Cfeas+Dfeas*L;
    P = lyap(AL',CL'*CL);
    eig_max = max(real(eig(A+B*Zfeas*G)));
    
    %First Stop Condition
    if eig_max <  0
        L = Zfeas*G;
        F = B*Zfeas;
        AL = A + B*L;
        CL = C + D*L;
        EF = E + F*H;
        P = lyap(AL',CL'*CL);
        Jfeas = trace(EF'*P*EF);
        flag = 1;
        break;
    end


    penalty_i = 2*sum(svd(Yfeas));
    penalty = penalty_i;
    f(it+1) = penalty;

    %Second Stop Condition
    if it > itmax
        flag = 3;
        break;
    end
end

if flag == 1
    return;  

elseif flag == 2
    disp('Max number of iterations');
    return;
end


end

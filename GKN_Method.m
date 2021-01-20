% This function can run any static or dynamic feedback problem design
% H2 objective and convex gain constraints
% The system under consideration is
%                   xdot = Ax+Bu+Ew
%                      y = Gx+Hw
%                      z = Cx+Du
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Jopt,Jfeas,Zopt,Zfeas] = GKN_Method(A,B,E,C,D,G,H,nc,options)

epsilon = options(1);

%Matrix Dimensions
nx = size(A,2);
nu = size(B,2);
ny = size(G,1);
nw = size(E,2);
nz = size(C,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Static Output Feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nc ==  0
    
    %Z structure
    Zstruct = [nu ny];
    
    %Step 3 - Feasibility Phase
    [Jfeas,Zfeas,P,f] = H2_Feas(A,B,E,C,D,G,H,nc,options,Zstruct,C,D);
    
    %Step 4 - Optimality Phase
    [Jopt,Zopt] = H2_Opt(A,B,E,C,D,G,H,nc,options,Zstruct,Jfeas,Zfeas,P);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dynamic Output Feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nc > 0
    
    %Z structure
    Zstruct = [nc nc;nc ny;nu nc];
    
    %Matrices with compatible dimensions
    Aa = [zeros(nc) zeros(nc,nx);zeros(nx,nc) A]; 
    Ba = [eye(nc) zeros(nc,nu);zeros(nx,nc) B]; 
    Ga = [eye(nc) zeros(nc,nx);zeros(ny,nc) G]; 
    Ea = [zeros(nc,nw) ; E]; 
    Ha = [zeros(nc,nw) ; H]; 
    Ca = [epsilon*eye(nc,nc) zeros(nc,nx);
         zeros(nc,nc) zeros(nc,nx);
         zeros(nz,nc) C];
    Da = [zeros(nc,nc) zeros(nc,nu); 
         epsilon*eye(nc,nc) zeros(nc,nu);
         zeros(nz,nc) D];
    Cafeas = [eye(nc,nc) ones(nc,nx);
             ones(nx,nc) eye(nx,nx);
             zeros(nc+nu,nc+nx)];
    Dafeas = [zeros(nc+nx,nc+nu); eye(nc+nu,nc+nu)];

    %Step 3 - Feasibility Phase
    [Jfeas,Zfeas,P,f] = H2_Feas(Aa,Ba,Ea,Cafeas,Dafeas,Ga,Ha,nc,options,Zstruct,Ca,Da);
    
    %Step 4 - Optimality Phase
    [Jopt,Zopt] = H2_Opt(Aa,Ba,Ea,Ca,Da,Ga,Ha,nc,options,Zstruct,Jfeas,Zfeas,P);
   
end

end

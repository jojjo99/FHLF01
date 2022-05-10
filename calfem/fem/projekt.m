load e.mat
load p.mat
load t.mat

coord=0.01*p';
enod=t(1:3,:)'; % nodes of elements
nelm=size(enod,1); % number of elements
nnod=size(coord,1); % number of nodes
dof=(1:nnod)'; % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number
for ie=1:nelm
edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
edof(ie,:)=[ie,enod(ie,:)];
end

er = e([1 2 5],:); % Reduced e

conv_segments = [10 11 12];   % Choosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
if ismember(er(3,i),conv_segments)
edges_conv = [edges_conv er(1:2,i)];
end
end

ndof= nelm +1;
Kt=zeros(nnod);
f=zeros(nnod,1);
kpmma=2.8;
kglas=0.8;
q=0;
alfa=100;
T0= 20+273;
Q=3*10^6;
nen=3;
ep=[0.5*0.01];
C=zeros(nnod);

D=eye(2);

[ex, ey]=coordxtr(edof, coord, dof, nen);
for elnr= 1:nelm
    
    if t(4,elnr)==1
        Q=0;
        kparam=kpmma;
        c_rho=1466*1185;
    else
        Q=3*10^6;
        kparam=kglas;
        c_rho=670*3860;
    end
    [k, fe]=flw2te(ex(elnr,:), ey(elnr,:), ep, kparam.*D, Q+(alfa*T0)/ep(1));
    kc=plantml(ex(elnr,:), ey(elnr,:), alfa);
    Ce=plantml(ex(elnr,:), ey(elnr,:), c_rho*ep);
    Kte=k+kc;
    indx = edof(elnr,2:end);
    Kt(indx,indx) = Kt(indx,indx)+Kte;
    C(indx,indx) = C(indx,indx)+Ce;
    f(indx) = f(indx) + fe;
    
end

%%
a = solveq(Kt,f);
ed=extract(edof,a);
patch(ex', ey', ed');

%%
a0=293*ones(nnod,1);
deltat=1;
a_prev=a0;
f_next=f;


for i=1:5
a_next= (C+deltat.*Kt)\(C*a_prev+deltat.*f_next);
a_prev=a_next;
end

eT=extract(edof,a_next);

figure()
patch(ex',ey',eT','EdgeColor','none')
title('Temperature distribution [C]')
colormap();
colorbar;
hold on
xlabel('x-position [m]')
caxis([290 300]);









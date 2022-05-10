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


for i=1:20
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

%% Del 3

ndof= nelm +1;
Kt=zeros(nnod*2);
f=zeros(nnod*2,1);
kpmma=2.8;
kglas=0.8;
q=0;
alfa=100;
T0= 20+273;
Q=3*10^6;
nen=3;
ep=[1 0.5*0.01];
C=zeros(nnod);
Eglas=2.8*10^9;
Nyglas=0.35;
Elins=67*10^9;
Nylins=0.2;


[ex, ey]=coordxtr(edof, coord, dof, nen);
for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t(4,elnr)==1
        D=hooke(1, Eglas, Nyglas);
        es=alfa*Eglas*deltaT/(1-Nyglas)*[1 1 0];
    else
         D=hooke(1, Elins, Nylins);
         es=alfa*Elins*deltaT/(1-Nylins)*[1 1 0];
    end
    [ke]=plante(ex(elnr,:), ey(elnr,:), ep, D);
    fe= plantf(ex(elnr,:), ey(elnr,:), ep, es);
    indx = edof_S(elnr,2:end);
    Kt(indx,indx) = Kt(indx,indx)+ke;
    f(indx) = f(indx) + fe;
   
end

vec1=[e(1,:) e(2,:)];
u=unique(vec1)';
u=sort(u);
bc(:,1)=u;
bc(:,2)=zeros(length(u),1);
bc(:,3)=zeros(length(u),1);
a=solveq(Kt,f,bc);
ed=extract(edof_S,a);

for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t(4,elnr)==1
        D=hooke(1, Eglas, Nyglas);
        es1=alfa*Eglas*deltaT/(1-Nyglas)*[1 1 0];
    else
         D=hooke(1, Elins, Nylins);
         es1=alfa*Elins*deltaT/(1-Nylins)*[1 1 0];
    end
    
    [es,et]=plants(ex(elnr,:),ey(elnr,:),ep,D, ed(elnr,:));
    Sigma(elnr,:)= es-es1;
    Tao(elnr,:)=et;
end
sigx=Sigma(:,1);
sigy=Sigma(:,2);
tauxy=Sigma(:,3);

Sigmatot=sqrt(sigx.^2 + sigy.^2 - sigx.*sigy + 3*tauxy);

ed=extract(edof,Sigmatot);
patch(ex', ey', ed');

 









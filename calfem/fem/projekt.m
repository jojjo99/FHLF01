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
%caxis([290 300]);

%% Del 3 Plot av displacement

Kt2=zeros(nnod*2);
f2=zeros(nnod*2,1);
kpmma=2.8;
kglas=0.8;
q=0;
T0= 20+273;
Q=3*10^6;
nen=3;
ep=[1 0.5*0.01];
C=zeros(nnod);
EPmma=2.8*10^9;
NyPmma=0.35;
Elins=67*10^9;
Nylins=0.2;
alphaPMMA=70*10^-6;
alphaGlass=7*10^-6;


[ex, ey]=coordxtr(edof, coord, dof, nen);
for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t(4,elnr)==1
        D=hooke(1, EPmma, NyPmma);
        es=alphaPMMA*EPmma*deltaT/(1-NyPmma)*[1 1 0];
    else
         D=hooke(1, Elins, Nylins);
         es=alphaGlass*Elins*deltaT/(1-Nylins)*[1 1 0];
    end
    ke=plante(ex(elnr,:), ey(elnr,:), ep, D);
    fe= plantf(ex(elnr,:), ey(elnr,:), ep, es);
    indx = edof_S(elnr,2:end);
    Kt2(indx,indx) = Kt2(indx,indx)+ke;
    f2(indx) = f2(indx) + fe;
   
end

er = e([1 2 5],:); % Reduced e
conv_segments = [5 4 6 1 8 2 7 3]; % Choosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
if ismember(er(3,i),conv_segments)
edges_conv = [edges_conv er(1:2,i)];
end
end
u=unique(edges_conv);
bc=[];
for i=1:length(u)
    bc= [bc; dof_S(u(i),1) 0; dof_S(u(i), 2) 0];
end

a2=solveq(Kt2,f2,bc);
%%a=zeros(length(a),1);
ed2=extract(edof_S,a2);

mag = 5000; % Magnification (due to small deformations)
exd = ex + mag*ed2(:,1:2:end);
eyd = ey + mag*ed2(:,2:2:end);
figure()
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3);
axis equal
title('Displacement field [Magnitude enhancement 100]')
%% Plot av von Mises stress

for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t(4,elnr)==1
        D=hooke(1, EPmma, NyPmma);
        es1=alphaPMMA*EPmma*deltaT/(1-NyPmma)*[1 1 0];
    else
         D=hooke(1, Elins, Nylins);
         es1=alphaGlass*Elins*deltaT/(1-Nylins)*[1 1 0];
    end
    
    [es,et]=plants(ex(elnr,:),ey(elnr,:),ep,D, ed2(elnr,:));
    Sigma(elnr,:)= es-es1;
end
sigx=Sigma(:,1);
sigy=Sigma(:,2);
tauxy=Sigma(:,3);


Sigmatot=sqrt(sigx.^2 + sigy.^2 - sigx.*sigy + 3*tauxy.^2);

for i=1:size(coord,1)
[c0,c1]=find(edof(:,2:4)==i);
Seff_nod(i,1)=sum(Sigmatot(c0))/size(c0,1);
end


ed3=extract(edof,Seff_nod);
figure()
patch(ex', ey', ed3');
colorbar;
%% Del 4.

T0= 20+273;
Q=3*10^6;
nen=3;
ep=[1 0.5*0.01];
C=zeros(nnod);
EPmma=2.8*10^9;
NyPmma=0.35;
Elins=67*10^9;
Nylins=0.2;
alphaPMMA=70*10^-6;
alphaGlass=7*10^-6;
Kt3=zeros(nnod*2);
[ex, ey]=coordxtr(edof, coord, dof, nen);
Kc = diag(ones(2,1),2) + diag(ones(2,1),-2)+ 2*diag(ones(4,1),0);

for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t(4,elnr)==1
        D=hooke(1, EPmma, NyPmma);
        es=alphaPMMA*EPmma*deltaT/(1-NyPmma)*[1 1 0];
    else
         D=hooke(1, Elins, Nylins);
         es=alphaGlass*Elins*deltaT/(1-Nylins)*[1 1 0];
    end
    
    ke=plante(ex(elnr,:), ey(elnr,:), ep, D);
    fe= plantf(ex(elnr,:), ey(elnr,:), ep, es);
    indx = edof_S(elnr,2:end);
    Kt2(indx,indx) = Kt2(indx,indx)+ke;
    f2(indx) = f2(indx) + fe;
   
    for i=1:length(er)
        if sum(ismember(er(1:2,i), edges_conv))==2
           ind1=er(1,i);
           ind2=er(2,i);
           L = sqrt((p(1,ind1)-p(1,ind2))^2 +(p(2,ind1)-p(2,ind2))^2);
           Kce = Kc*L/6;
           indx1=dof_S(ind1,:);
           indx2=dof_S(ind2,:);
           indtot=[indx1 indx2];
           Kt3(indtot,indtot)=Kt3(indtot,indtot)+Kce;
        end
    end
    
    
end

er = e([1 2 5],:); % Reduced e
conv_segments = [5 4 6 1 8 2 7 3]; % Choosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
if ismember(er(3,i),conv_segments)
edges_conv = [edges_conv er(1:2,i)];
end
end
u=unique(edges_conv);
bc=[];
for i=1:length(u)
    bc= [bc; dof_S(u(i),1) 0; dof_S(u(i), 2) 0];
end
Ktot=Kt2-Kt3;
a3=solveq(Ktot,f2);

ed3=extract(edof_S,a3);

mag = 100; % Magnification (due to small deformations)
exd = ex + mag*ed3(:,1:2:end);
eyd = ey + mag*ed3(:,2:2:end);
figure()
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3);
axis equal
title('Displacement field [Magnitude enhancement 100]')

ed2=ed3;

for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t(4,elnr)==1
        D=hooke(1, EPmma, NyPmma);
        es1=alphaPMMA*EPmma*deltaT/(1-NyPmma)*[1 1 0];
    else
         D=hooke(1, Elins, Nylins);
         es1=alphaGlass*Elins*deltaT/(1-Nylins)*[1 1 0];
    end
    
    [es,et]=plants(ex(elnr,:),ey(elnr,:),ep,D, ed2(elnr,:));
    Sigma(elnr,:)= es-es1;
end
sigx=Sigma(:,1);
sigy=Sigma(:,2);
tauxy=Sigma(:,3);


Sigmatot=sqrt(sigx.^2 + sigy.^2 - sigx.*sigy + 3*tauxy.^2);

for i=1:size(coord,1)
[c0,c1]=find(edof(:,2:4)==i);
Seff_nod(i,1)=sum(Sigmatot(c0))/size(c0,1);
end


ed3=extract(edof,Seff_nod);
figure()
patch(ex', ey', ed3');
colorbar;
%Loopa i elementen,
% Kolla om tv? noder p? rand
% Calc length, 
% Calc Kce
%Assembla
%sum(isember(enod(i,:), edges_conv(
% dof_s

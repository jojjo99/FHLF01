%% Project FEM %%
% Code for using the symmetry boundary conditions.
% Start reading the files t_sym.mat, p_sym.mat, e_sym.mat

%% Part A. Stationary temperature distribution
load 't_sym.mat'
load 'p_sym.mat'
load 'e_sym.mat'

coord=0.01*p_sym'; % Correct coordinates to meters
enod=t_sym(1:3,:)'; % Nodes of elements
nelm=size(enod,1); % Number of elements
nnod=size(coord,1); % Number of nodes
dof=(1:nnod)'; % Dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % Give each dof a number

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end
er = e_sym([1 2 5],:); % Reduced e

conv_segments = [10 11 12];   % Choosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
end

% Initializing constants and matrices
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
    
    if t_sym(4,elnr)==1
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
patch(ex', -ey'-0.01, ed');
colorbar;

%% Part B. Transient temperature distribution

%%
a0=293*ones(nnod,1); % Set temperature in all nodes to 20 degrees
deltat=1; % 1s timestep
a_prev=a0;
f_next=f;


for i=1:5
    a_next= (C+deltat.*Kt)\(C*a_prev+deltat.*f_next);
    a_prev=a_next;
end

eT=extract(edof,a_next);

figure()
patch(ex',ey',eT','EdgeColor','none')
patch(ex',-0.01-ey',eT','EdgeColor','none')
title('Temperature distribution [C]')
colormap();
colorbar;
hold on
xlabel('x-position [m]')
caxis([290 300]);

%% Part C. Effective von Mises stress field

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
    
    if t_sym(4,elnr)==1
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

er = e_sym([1 2 5],:); % Reduced e
conv_segments = [7 2 8 3 1]; % Choosen boundary segments
edges_conv1 = [];
for i = 1:size(er,2)
if ismember(er(3,i),conv_segments)
edges_conv1 = [edges_conv1 er(1:2,i)];
end
end

er = e_sym([1 2 5],:); % Reduced e
conv_segments = [4 5 6]; % Choosen boundary segments
edges_conv2 = [];
for i = 1:size(er,2)
if ismember(er(3,i),conv_segments)
edges_conv2 = [edges_conv2 er(1:2,i)];
end
end

u1=unique(edges_conv1);
u2=unique(edges_conv2);
bc=[];
for i=1:length(u1)
    bc= [bc; dof_S(u1(i),1) 0; dof_S(u1(i), 2) 0];
end

for i=1:length(u2)
    bc= [bc;dof_S(u2(i),2) 0];
end

bc1=unique(bc(:,1));
bc1=[bc1, zeros(length(bc1),1)];

a2=solveq(Kt2,f2,bc1);
%%a=zeros(length(a),1);
ed2=extract(edof_S,a2);

mag = 5000; % Magnification (due to small deformations)
exd = ex + mag*ed2(:,1:2:end);
eyd = ey + mag*ed2(:,2:2:end);
figure()
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3);

patch(exd',-0.01-eyd',[0 0 0],'EdgeColor', 'none','FaceAlpha',0.3);
hold on
patch(exd',-0.01-eyd',[0 0 0],'FaceAlpha',0.3);

axis equal
title('Displacement field [Magnitude enhancement 100]')
%% Part C. Plot of von Mises stress field.
for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t_sym(4,elnr)==1
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
patch(ex', -0.01-ey', ed3');
colorbar;
%% Part D. Spring boundaries.
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
    
    if t_sym(4,elnr)==1
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
        if sum(ismember(er(1:2,i), edges_conv1))==2
           ind1=er(1,i);
           ind2=er(2,i);
           L = sqrt((p_sym(1,ind1)-p_sym(1,ind2))^2 +(p_sym(2,ind1)-p_sym(2,ind2))^2);
           Kce = Kc*L/6;
           indx1=dof_S(ind1,:);
           indx2=dof_S(ind2,:);
           indtot=[indx1 indx2];
           Kt3(indtot,indtot)=Kt3(indtot,indtot)+Kce;
        end
    end
    
    
end

er = e_sym([1 2 5],:); % Reduced e
conv_segments = [4 5 6]; % Choosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
if ismember(er(3,i),conv_segments)
edges_conv = [edges_conv er(1:2,i)];
end
end
u=unique(edges_conv);
bc=[];
for i=1:length(u)
    bc= [bc;dof_S(u2(i),2) 0];
end
Ktot=Kt2-Kt3;
a3=solveq(Ktot,f2,bc);

ed3=extract(edof_S,a3);

mag = 100; % Magnification (due to small deformations)
exd = ex + mag*ed3(:,1:2:end);
eyd = ey + mag*ed3(:,2:2:end);
figure()
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3);
patch(ex',-0.01-ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',-0.01-eyd',[0 0 0],'FaceAlpha',0.3);
axis equal
title('Displacement field [Magnitude enhancement 100]')

ed2=ed3;

for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t_sym(4,elnr)==1
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
patch(ex', -0.01-ey', ed3');
colorbar;


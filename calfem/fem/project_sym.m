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
nen=3;
Kt=zeros(nnod);
f=zeros(nnod,1);
kPMMA=2.8;
kGLASS=0.8;
alpha=100;
T0= 293;
Q=3*10^6;
ep=[0.5*0.01];
C=zeros(nnod);
D=eye(2);
[ex, ey]=coordxtr(edof, coord, dof, nen);

% Calculating and assembling K and f matrices
for elnr= 1:nelm
    
    if t_sym(4,elnr)==1
        Q=0;
        kparam=kPMMA;
        c_rho=1466*1185;
    else
        Q=3*10^6;
        kparam=kGLASS;
        c_rho=670*3860;
    end
    [k, fe]=flw2te(ex(elnr,:), ey(elnr,:), ep, kparam.*D, Q+(alpha*T0)/ep(1));
    kc=plantml(ex(elnr,:), ey(elnr,:), alpha); %% element convection k matrix
    Ce=plantml(ex(elnr,:), ey(elnr,:), c_rho*ep); %%for time integration
    Kte=k+kc;
    indx = edof(elnr,2:end);
    Kt(indx,indx) = Kt(indx,indx)+Kte;
    C(indx,indx) = C(indx,indx)+Ce;
    f(indx) = f(indx) + fe;
    
end

% Solving the diffenrential equation
a = solveq(Kt,f);
ed=extract(edof,a);
patch(ex', ey', ed', 'Edgecolor', 'none');
patch(ex', -ey'-0.01, ed', 'Edgecolor', 'none');
colormap();
cbar=colorbar;
xlabel('x-position [m]');
ylabel('y-position [m]');
cbar.Title.String = 'Temperature [K]';
str = 'Label2';
set(get(cbar, 'xlabel'), 'string', str, 'rotation', 0);
%% Part B. Transient temperature distribution

%%
a0=293*ones(nnod,1); % Set temperature in all nodes to 293 K
delta_t=1; % 1s timestep
a_prev=a0;
f_next=f;

% Numerical integration
for i=1:20
    a_next= (C+delta_t.*Kt)\(C*a_prev+delta_t.*f_next);
    a_prev=a_next;
end

eT=extract(edof,a_next);

figure()
patch(ex',ey',eT','EdgeColor','none')
patch(ex',-0.01-ey',eT','EdgeColor','none')
colormap();
cbar=colorbar;
hold on
xlabel('x-position [m]')
ylabel('y-position [m]')
caxis([290 300]);
cbar.Title.String = 'Temperature [K]';
str = 'Label2';
set(get(cbar, 'xlabel'), 'string', str, 'rotation', 0);

%% Part C. Displacement field, fixed boundaries.

Kt2=zeros(nnod*2);
f2=zeros(nnod*2,1);
kPMMA=2.8;
kGLASS=0.8;
T0= 293;
Q=3*10^6;
nen=3;
ep=[1 0.5*0.01];
C=zeros(nnod);
EPMMA=2.8*10^9;
NyPMMA=0.35;
EGLASS=67*10^9;
NyGLASS=0.2;
alphaPMMA=70*10^-6;
alphaGLASS=7*10^-6;


[ex, ey]=coordxtr(edof, coord, dof, nen);
for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t_sym(4,elnr)==1
        D=hooke(1, EPMMA, NyPMMA);
        es=alphaPMMA*EPMMA*deltaT/(1-NyPMMA)*[1 1 0];
    else
         D=hooke(1, EGLASS, NyGLASS);
         es=alphaGLASS*EGLASS*deltaT/(1-NyGLASS)*[1 1 0];
    end
    ke=plante(ex(elnr,:), ey(elnr,:), ep, D);
    fe= plantf(ex(elnr,:), ey(elnr,:), ep, es);
    indx = edof_S(elnr,2:end);
    Kt2(indx,indx) = Kt2(indx,indx)+ke;
    f2(indx) = f2(indx) + fe;
   
end

% Finding nodes on the boundary
er = e_sym([1 2 5],:); % Reduced e
conv_segments = [7 2 8 3 1]; % Choosen boundary segments
edges_conv1 = [];
for i = 1:size(er,2)
if ismember(er(3,i),conv_segments)
edges_conv1 = [edges_conv1 er(1:2,i)];
end
end

% Finding nodes on the symmetry line
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

% Setting the degree of freedom in y-direction to zero
for i=1:length(u1)
    bc= [bc; dof_S(u1(i),1) 0; dof_S(u1(i), 2) 0];
end

for i=1:length(u2)
    bc= [bc;dof_S(u2(i),2) 0];
end

bc1=unique(bc(:,1));
bc1=[bc1, zeros(length(bc1),1)];

a2=solveq(Kt2,f2,bc);
ed2=extract(edof_S,a2);

mag = 2000; % Magnification (due to small deformations)
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
xlabel('x-position [m]')
ylabel('y-position [m]')

%% Part C. Effective von Mises stress field.
ed2=extract(edof_S,a2);
for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t_sym(4,elnr)==1
        D=hooke(1, EPMMA, NyPMMA);
        es1=alphaPMMA*EPMMA*deltaT/(1-NyPMMA)*[1 1 0];
    else
         D=hooke(1, EGLASS, NyGLASS);
         es1=alphaGLASS*EGLASS*deltaT/(1-NyGLASS)*[1 1 0];
    end
    
    [es,et]=plants(ex(elnr,:),ey(elnr,:),ep,D, ed2(elnr,:));
    Sigma(elnr,:)= es-es1;
end
sigx=Sigma(:,1);
sigy=Sigma(:,2);
tauxy=Sigma(:,3);


Sigmatot=sqrt(sigx.^2 + sigy.^2 - sigx.*sigy + 3*tauxy.^2);

% Calculating stress of the nodal points
for i=1:size(coord,1)
[c0,c1]=find(edof(:,2:4)==i);
Seff_nod(i,1)=sum(Sigmatot(c0))/size(c0,1);
end

ed3=extract(edof,Seff_nod);
ed3=ed3./10^6;

figure()
patch(ex', ey', ed3', 'Edgecolor', 'none');
patch(ex', -0.01-ey', ed3', 'Edgecolor', 'none');
cbar=colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
cbar.Title.String = 'Stress [MPa]';

%% Part D. Spring boundaries.
T0= 20+273;
Q=3*10^6;
nen=3;
ep=[1 0.5*0.01];
C=zeros(nnod);
EPMMA=2.8*10^9;
NyPMMA=0.35;
EGLASS=67*10^9;
NyGLASS=0.2;
alphaPMMA=70*10^-6;
alphaGLASS=7*10^-6;
Kt3=zeros(nnod*2);
[ex, ey]=coordxtr(edof, coord, dof, nen);
Kc = diag(ones(2,1),2) + diag(ones(2,1),-2)+ 2*diag(ones(4,1),0);
kspring=-1000;


for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t_sym(4,elnr)==1
        D=hooke(1, EPMMA, NyPMMA);
        es=alphaPMMA*EPMMA*deltaT/(1-NyPMMA)*[1 1 0];
    else
         D=hooke(1, EGLASS, NyGLASS);
         es=alphaGLASS*EGLASS*deltaT/(1-NyGLASS)*[1 1 0];
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
           Kce = Kc*L/6*kspring*0.005;
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

ed4=extract(edof_S,a3);

mag = 100; % Magnification (due to small deformations)
exd = ex + mag*ed4(:,1:2:end);
eyd = ey + mag*ed4(:,2:2:end);
figure()
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3);
patch(ex',-0.01-ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exd',-0.01-eyd',[0 0 0],'FaceAlpha',0.3);
axis equal
xlabel('x-position [m]')
ylabel('y-position [m]')

ed2=ed4;

for elnr= 1:nelm
    
    deltaT=sum(ed(elnr,:))/3 - T0;
    
    if t_sym(4,elnr)==1
        D=hooke(1, EPMMA, NyPMMA);
        es1=alphaPMMA*EPMMA*deltaT/(1-NyPMMA)*[1 1 0];
    else
         D=hooke(1, EGLASS, NyGLASS);
         es1=alphaGLASS*EGLASS*deltaT/(1-NyGLASS)*[1 1 0];
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

ed5=extract(edof,Seff_nod);
ed5=ed5./10^6;

figure()
patch(ex', ey', ed5', 'Edgecolor', 'none');
patch(ex', -0.01-ey', ed5', 'Edgecolor', 'none');
cbar=colorbar;
title('Effective von Mises stress field [Pa]')
xlabel('x-position [m]')
ylabel('y-position [m]')
cbar.Title.String = 'Stress [MPa]';


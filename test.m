A=10;
k=5;
L=6;
Q=100;
nen=2; %%antal noder/element
nelm=3; %%antal element
ndof= nelm +1; %%antal frihetsgrader ges av totalt antal noder i detta fall

coord=linspace(2,8,ndof)';
vec=(1:nelm)';
edof=[vec, vec, vec+1];  %%frihetsgrader per element
dof=(1:ndof)';  %%vilken frihetsgrad varje nod har

[Ex]=coordxtr(edof, coord, dof, nen); %%koordinat för varje nod

K=zeros(ndof);
f=zeros(ndof, 1);

eleL=L/nelm;
for elnr= 1:nelm
    ke=spring1e(k*A/eleL);
    fa=eleL*Q/2;
    fe=[fa;fa];
    [K,f]=assem(edof(elnr,:), K, ke, f,fe);
end

f(end)=f(end)+15*A;
bc=[1 0];
a=solveq(K,f,bc);
plot(coord,a);

%%
nen=3;

K=zeros(ndof);
f=zeros(ndof,1);
D=eye(2);
ep=[1];

for elnr= 1:nelm
    ke=flw2te(ex(elnr,:), ey(elnr,:), ep, D);
    K=assem(edof(elnr,:), K, ke);
    
end
a=solveq(K,f,bc);
ed=extract(edof,a);
patch(ex', ey', ed');



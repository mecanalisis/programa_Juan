% Programa de calculo de cojinetes.
% primero cargar la base de datos Bearing.N
% N es la velocidad de rotacion en rps, es un vectos de 100:2000 cada 100
% T1 es la tempoeratura de entrada del aceite
% T2 representa la temperatura de salida del aceite
% T3 intenta estimar la temperatura de Metal 
% Aceite tipo ISO46 siacosidad [46 6]cSt a [40 100]C



D=.445;
L=.75*D;
r=D/2;
Cb=0.0015*r; %huelgo radial del cojinete
W=30000*9.81; %N
N=(100:100:2000)'/60; %rps
densidad1=850; %Kg/m^3
Cesp=1670; %J/Kg/K % valor tipico para aceite mineral
alpha=0.2;

%---------------Calculo
LoverD=L/D;
Cp=Cb*2;
phi=Cp/r;
m=1-Cb/Cp;
omega=N*2*pi;

T1=40; %temperatura de entrada de aceite
[myu,densidad,Cv,v2]=VisCal([46 6],densidad1,[40 100],T1);
OutVar=[];
OutKC=[];
for ii=1:length(omega)
    T2=1;
    T2new=3*T2;
    k=0;
    temp=[];
    while abs(T2-T2new)>1
        T2=T2new;
        S=myu*N(ii)*D*L/phi^2/W;
        
        Y=interp1(Bearing.Som,Bearing.Static,S,'spline');
        e=Y(1);theta=Y(2);Qs1=Y(3);Qe1=Y(4);fjPhi=Y(5);
        
        fj=phi*fjPhi;
        Qs=r*omega(ii)*Cp*L*Qs1;
        Qe=r*omega(ii)*Cp*L*Qe1;
        
        T2new=T1+fj*r*W*omega(ii)/(Cesp*Qs*densidad);
       
        T=alpha*T1+(1-alpha)*T2new; %temp de trabajo
        temp=[temp T];
        [myu,densidad,Cv,v2]=VisCal([46 6],850,[40 100],T);
        k=k+1;
        if k==200; break; end
    end
    T3=T1+fj*r*W*omega(ii)/(Cesp*Qe*densidad);
    Y=interp1(Bearing.Som,Bearing.KC,S,'spline');
    K=Y(1:4)*W/Cp;
    C=Y(5:8)*W/Cp/omega(ii);    
    OutVar(ii,:)=[N(ii)*60 e theta Qs Qe T1 T2 T3 T myu densidad];
    OutKC(ii,:)=[K C];
end

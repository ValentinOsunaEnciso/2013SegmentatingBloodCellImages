% Differential Evolution + Mixture of Gaussian Functions (DE_MGF) %%%%%%%%%
% Valentin Osuna-Enciso, CIC-IPN, Abril, 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper: Segmentation of blood cell images using evolutive Methods. %%%%%%%
% En esta version, utilizo distancia Hellinger en vez de distancia %%%%%%%%
% Euclideana. Mejores resultados, mejor convergencia, no es necesario el %%
% factor de penalizacion. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [temp t EvaluacioneFO tiempo]=DE_MGF
clear all
format long
DB=imread('Im019_1.jpg'); %images from ALL_IDB
numClases=3;    %Numero de clases que deseo obtener(a-priori)
DB=rgb2gray(DB);%Convierto a escala de grises
H=imhist(DB);   %Calculo del Histograma
H=H/sum(H);     %Se normaliza Histograma experimental(suma de Hi=1)
Amax=max(H);
L=size(H,1);
x_high=[Amax;Amax;Amax;L-1;L-1;L-1;(L-1)/12;(L-1)/12;(L-1)/12];
x_low=[0;0;0;1;1;1;0;0;0];x=[]; xp=0:1:255;
%INICIALIZACION DE1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np=90;          %Tamano de la poblacion.90
F=0.25;         %Factor de escalamiento
Cr=0.8;         %Probabilidad de cruza
D=numClases*3;  %Para 3 clases, 9 dimensiones.
gmax=400;       %Numero maximo de iteraciones
t=0;            %Contador de iteraciones
EvaluacioneFO=0;
for ind1=1:Np
    for ind2=1:D
       x(ind2,ind1)=x_low(ind2,1)+rand()*(x_high(ind2,1)-x_low(ind2,1));  
    end
    error(ind1)=MGF(x(:,ind1),H);
    EvaluacioneFO=EvaluacioneFO+1;
end
[f_x ind]=sort(error);%Ordeno deMenor aMayor la evaluacion de la funcion
x_best=x(:,ind(1));%Guardo el mejor;Van Np Evaluaciones De La Funcion.
% rng('shuffle');
tic
while f_x(1,1)>0.1186 && t<gmax 
   for ind1=1:Np
       r1=randi(Np);r2=randi(Np);
       while(r1==r2==ind1)
          r1=randi(Np); r2=randi(Np);       %Generados sean diferentes.
       end
       v(:,1)=x_best+F*(x(:,r1)-x(:,r2));   %Mutacion
       u(:,1)=x(:,ind1);
       j_rand=randi(D);
       for ind2=1:D                         %Genero vector de prueba
          if(rand()<Cr || ind2==j_rand) && v(j_rand,1)>0     %Cruza
             u(j_rand,1)=v(j_rand,1);
          end
       end
        temp=MGF(u(:,1),H);
        EvaluacioneFO=EvaluacioneFO+1;
        if(temp<error(ind1))
          x(:,ind1)=u(:,1);
          error(ind1)=temp;
          if (temp<f_x(1,1))
             x_best=u(:,1);
             f_x(1,1)=temp;
          end
      end
   end
   t=t+1;
   %disp(sprintf('  Iterac=%d, Fitnes=%f, Evaluaciones:%d\n',t,f_x(1,1),EvaluacioneFO));
end
disp(sprintf('  Iterac=%d, Fitnes=%f, Evaluaciones:%d\n',t,f_x(1,1),EvaluacioneFO));
tiempo=toc;
%x_best=x_best';
Haus=grafica(x_best,H,DB);
temp=f_x(1,1);
end
%% FUNCION GRAFICA GAUSSIANAS E IMAGEN SEGMENTADA: %%%%%%%%%%%%%%%%%%%%%%%%
% Recibe x_best: La mejor particula, D: dimensiones de cada particula
% H: histograma de la imagen, DB: imagen en escala de gris
function mhd=grafica(x_best,H,DB)
    xp=0:1:255;
    valM1=round(x_best(4,1));valM2=round(x_best(5,1));valM3=round(x_best(6,1));
    valA1=x_best(1,1);valA2=x_best(2,1);valA3=x_best(3,1);
    valDE1=x_best(7,1);valDE2=x_best(8,1);valDE3=x_best(9,1);
    [M ind1]=sort([valM1,valM2,valM3]);
    if(ind1(1)==1)
        DE1=valDE1;A1=valA1;M1=valM1;
    elseif(ind1(1)==2)
        DE1=valDE2;A1=valA2;M1=valM2;
    else
        DE1=valDE3;A1=valA3;M1=valM3;
    end
    if(ind1(2)==1)
        DE2=valDE1;A2=valA1;M2=valM1;
    elseif(ind1(2)==2)
        DE2=valDE2;A2=valA2;M2=valM2;
    else
        DE2=valDE3;A2=valA3;M2=valM3;
    end
    if(ind1(3)==1)
        DE3=valDE1;A3=valA1;M3=valM1;
    elseif(ind1(3)==2)
        DE3=valDE2;A3=valA2;M3=valM2;
    else
        DE3=valDE3;A3=valA3;M3=valM3;
    end
    Resultado=(A1*exp(-((xp-M1).^2)/(2*(DE1^2))))+...
        (A2*exp(-((xp-M2).^2)/(2*(DE2^2))))+...
        (A3*exp(-((xp-M3).^2)/(2*(DE3^2))));
    plot(Resultado,'k--','LineWidth',2)%,Hold on,plot((ampli1*exp(-((x-media1).^2)/(2*(DE1^2))))),plot((ampli2*exp(-((x-media2).^2)/(2*(DE2^2))))),plot((ampli3*exp(-((x-media3).^2)/(2*(DE3^2)))))
    hold on
    plot(H,'r','LineWidth',2),figure
    plot((A3*exp(-((xp-M3).^2)/(2*(DE3^2)))),'k--','LineWidth',2),hold on
    plot((A2*exp(-((xp-M2).^2)/(2*(DE2^2)))),'k-.','LineWidth',2)
    plot((A1*exp(-((xp-M1).^2)/(2*(DE1^2)))),'k','LineWidth',2)
    %plot(Resultado,'k'),title('Resultado')
%%%%%Realizo umbralizacion imagen escala de grises:%%%%%%%%%%%%%%%%%%%%%%%
    a1=(DE1^2)-(DE2^2);
    a2=(DE2^2)-(DE3^2);
    b1=2*((M1*(DE2^2))-(M2*(DE1^2)));
    b2=2*((M2*(DE3^2))-(M3*(DE2^2)));
    c1=((DE1*M2)^2)-((DE2*M1)^2)+(2*((DE1*DE2)^2)*log((DE2*A1)/(DE1*A2)));
    c2=((DE2*M3)^2)-((DE3*M2)^2)+(2*((DE3*DE2)^2)*log((DE3*A2)/(DE2*A3)));
    T1a=(-b1+sqrt((b1^2)-(4*a1*c1)))/(2*a1);
    T1b=(-b1-sqrt((b1^2)-(4*a1*c1)))/(2*a1);
    T2a=(-b2+sqrt((b2^2)-(4*a2*c2)))/(2*a2);
    T2b=(-b2-sqrt((b2^2)-(4*a2*c2)))/(2*a2);
    [fila columna]=size(DB);
    for ind1=1:fila
        for ind2=1:columna
            if (DB(ind1,ind2)<=T1b)&&(DB(ind1,ind2)>=0)
                DBsegmented(ind1,ind2)=0;
            elseif (DB(ind1,ind2,1)<=T2b)&&(DB(ind1,ind2)>T1b)
                DBsegmented(ind1,ind2)=.5;
            elseif(DB(ind1,ind2,1)>T2b)
                DBsegmented(ind1,ind2)=1;
            end
        end
    end
    figure,DBsegmented=mat2gray(DBsegmented);
    imshow(DBsegmented)
    %% Modified Hausdorff distance:
    AI=DBsegmented;
    BI=imread('Im248_0_GT.tif');
    BI=rgb2gray(BI(:,:,1:3));
    BI=BI>=254;
    [AIbordes t]=edge(AI,'canny',.1);
    [BIbordes t]=edge(BI,'canny',.1);
    [A(:,1) A(:,2)]=find(AIbordes);
    [B(:,1) B(:,2)]=find(BIbordes);
    [ mhd ] = ModHausdorffDist( A, B );
end
%% Mixture of Gaussian Functions: it works with 3 Gaussian functions. %%%%%
function error=MGF(x,H)
    xp=0:1:255;
    mix(:,1)=(x(1,1)*exp(-((xp-round(x(4,1))).^2)...
   /(2*(x(7,1)^2))))+(x(2,1)*exp(-((xp-round(x(5,1))).^2)...
   /(2*(x(8,1)^2))))+(x(3,1)*exp(-((xp-round(x(6,1))).^2)...
   /(2*(x(9,1)^2))));
%     medias=[round(x(4,1)),round(x(5,1)),round(x(4,1))];
%     medias=sort(medias);
%     error=HellingerDistance(H,mix)+(3/((medias(3)-medias(2))+(medias(2)-medias(1))+eps));
    error=HellingerDistance(H,mix);
end
%% 
    %[M temp]=sort(round(x(4:6,1)));
    %penalizacion=w*abs(sum(mix)-1);            % Original
    %penalizacion=1/((M(3)-M(1))+(M(2)-M(1)));  % Modificada
    % Distancia Euclideana:
%     error=((H-mix)).^2;
%     error=sum(error)/size(H,1);
%     error=error+penalizacion;
%     error1=MahalanobisDistance(H,mix);
%     error2=BhattacharyyaDistance(H,mix);
    % Distancia Minkowski:
%   error4=abs(H-mix);[error4 temp]=max(error4);error4=error4+penalizacion;
%     error4=abs(H-mix);error4=sum(error4)+penalizacion;
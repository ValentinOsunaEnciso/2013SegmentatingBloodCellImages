% Artificial Bee Colony + Mixture of Gaussian Functions (ABC_MGF)  %%%%%%%%
% Valentin Osuna-Enciso, CIC-IPN, Abril, 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper: Segmentation of blood cell images using evolutive Methods. %%%%%%%
% En esta version, utilizo distancia Hellinger en vez de distancia %%%%%%%%
% Euclideana. Mejores resultados, mejor convergencia, no es necesario el %%
% factor de penalizacion. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [temp k EvaluacioneFO tiempo]=ABC_MGF
%/* ABC algorithm coded using MATLAB language */
%/* Artificial Bee Colony (ABC) is one of the most recently defined algorit
%hms by Dervis Karaboga in 2005, motivated by the intelligent behavior of 
%honey bees. Referance Papers:
%1)D. Karaboga,AN IDEA BASED ON HONEY BEE SWARM FOR NUMERICAL OPTIMIZATION,
%TECHNICAL REPORT-TR06, Erciyes University, Engineering Faculty, Computer 
%Engineering Department 2005. 2)D. Karaboga, B. Basturk, A powerful and 
%Efficient Algorithm for Numerical Function Optimization: Artificial Bee 
%Colony (ABC) Algorithm, Journal of Global Optimization, Volume:39,Issue:3,
%pp:459-171, November 2007,ISSN:0925-5001 , doi: 10.1007/s10898-007-9149-x 
%3)D. Karaboga, B. Basturk, On The Performance Of Artificial Bee Colony 
%(ABC) Algorithm, Applied Soft Computing,Volume 8, Issue 1, January 2008, 
%Pages 687-697. 4)D. Karaboga, B. Akay, A Comparative Study of Artificial 
%Bee Colony Algorithm,  Applied Mathematics and Computation, 214, 108-132, 
%2009.
%Copyright,2009,ErciyesUniversity,IntelligentSystemsResearchGroup,TheDept. 
%of Computer Engineering. Contact:Dervis Karaboga (karaboga@erciyes.edu.tr)
%Bahriye Basturk Akay (bahriye@erciyes.edu.tr)
clear all
DB=imread('Im019_1.jpg'); 
numClases=3;    %Numero de clases que deseo obtener(a-priori)
DB=rgb2gray(DB);%Convierto a escala de grises
H=imhist(DB);   %Calculo del Histograma
H=H/sum(H);     %Se normaliza Histograma experimental(suma de Hi=1)
Amax=max(H);    %Alturas maximas
L=size(H,1);    % Numero de niveles de gris 
% CONTROL PARAMETERS OF ABC ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NP=90;          %The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2;%TheNumberOfFoodSourcesEqualsTheHalf of the colony size
%AFoodSourceWhichCouldn'tBeImprovedThrough"limit"trialsIsAbandonedByIts
limit=70; %employed bee.
Nmax=600; %/*The number of cycles for foraging {a stopping criteria}*/
%/* Problem specific variables*/
%objfun='Sphere'; %cost function to be optimized
D=numClases*3;  % Num de parametros de funciones Gaussianas(alt,med,desv)
x_high=[Amax;Amax;Amax;L-1;L-1;L-1;(L-1)/12;(L-1)/12;(L-1)/12];
x_low=[0;0;0;1;1;1;0;0;0];x=[]; xp=0:1:255;
EvaluacioneFO=0;
for ind1=1:FoodNumber
    for ind2=1:D 
        Foods(ind1,ind2)=(x_low(ind2,1)+rand()*...
           (x_high(ind2,1)-x_low(ind2,1))); 
    end
    F_x_(ind1)=MGF(Foods(ind1,:)',H);
    EvaluacioneFO=EvaluacioneFO+1;
end
%Foods [FoodNumber][D]; /*Foods is the population of food sources. 
%Each row of Foods matrix is a vector holding D parameters to be optimized. 
%The number of rows of Foods matrix equals to the FoodNumber
%F_x_[FoodNumber];  /*f is a vector holding objective function values 
%associated with food sources. Fitness[FoodNumber]; fitness is a vector 
%holding fitness (quality) values associated with food sources*/
%trial[FoodNumber]; /*trial is a vector holding trial numbers through which
%solutions can not be improved. prob[FoodNumber]; /*prob is a vector 
%holding probabilities of food sources (solutions) to be chosen*/
%solution [D]; /*New solution (neighbour) produced by 
%v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter 
%and k is a randomlu chosen solution different from i.ObjValSol; Objective 
%function value of new solution. FitnessSol; Fitness value of new solution
%neighbour, param2change; /*param2change corrresponds to j, neighbour 
%corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
%GlobalMin; /*Optimum solution obtained by ABC algorithm*/
%GlobalParams[D]; /*Parameters of the optimum solution*/
%GlobalMins[runtime];GlobalMins holdsThe GlobalMin ofEach runInMultipleRuns
%GlobalMins=zeros(1,runtime); 
% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has 
%different range, use arrays lb[j], ub[j] instead of lb and ub 
% Range = repmat((ub-lb),[FoodNumber 1]);
% Lower = repmat(lb, [FoodNumber 1]);
% Foods = rand(FoodNumber,D) .* Range + Lower;
%F_x_=feval(objfun,Foods);
Fitness=calculateFitness(F_x_*10000);
%reset trial counters
trial=zeros(1,FoodNumber);
%/*The best food source is memorized*/
BestInd=find(F_x_==min(F_x_));
BestInd=BestInd(end);
GlobalMin=F_x_(BestInd);
GlobalParams=Foods(BestInd,:);
k=1;
% rng('shuffle');
tic
while (k <= Nmax && GlobalMin>0.1186)
    %%%%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for ind1=1:(FoodNumber)        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;        
       %ArandomlyChosenSolutionIsUsedInProducingaMutantSolutionOfsolution i
       neighbour=fix(rand*(FoodNumber))+1;       
        %/*Randomly selected solution must be different from the solution i       
            while(neighbour==ind1)
                neighbour=fix(rand*(FoodNumber))+1;
            end;        
       sol=Foods(ind1,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(ind1,Param2Change)+...
           (Foods(ind1,Param2Change)-Foods(neighbour,Param2Change))*...
           (rand-0.5)*2;        
       %IfGeneratedParameterValueIsOutOfBoundaries,ItIsShiftedOntoTheBounda
       ind=find(sol<x_low(:,1)');
       sol(ind)=x_low(ind,1)';
       ind=find(sol>x_high(:,1)');
       sol(ind)=x_high(ind,1)';        
       %Evaluate new solution
       ObjValSol=MGF(sol',H);
       EvaluacioneFO=EvaluacioneFO+1;
       FitnessSol=calculateFitness(ObjValSol*10000);        
       %GreedySelectionIsAppliedBetweenTheCurrentSolution i and its mutant
       %IfTheMutantSolutionIsBetterThanTheCurrent solution i, replace the 
       %solution with the mutant and reset the trial counter of solution i
       if (FitnessSol>Fitness(ind1)) %¿>, o <? Originalmente >
          Foods(ind1,:)=sol;
          Fitness(ind1)=FitnessSol;
          F_x_(ind1)=ObjValSol;
          trial(ind1)=0;
       else
          %Ifsolution i can'tBeImproved,IncreaseTrialCounter  
          trial(ind1)=trial(ind1)+1;
       end;                  
    end;
    % CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A food source is chosen with the probability which is proportioal to 
    % its quality Different schemes can be used to calculate the 
    % probability values. For example prob(i)=fitness(i)/sum(fitness) or
    % in a way used in the metod below prob(i)=a*fitness(i)/max(fitness)+b
    % probability values are calculated by using fitness values and 
    % normalized by dividing maximum fitness value*/
    prob=(0.9.*Fitness./max(Fitness))+0.1;  
    % ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind1=1;
    t=0;
    while(t<FoodNumber)
        if(rand<prob(ind1))
            t=t+1;
            %/*The parameter to be changed is determined randomly*/
            Param2Change=fix(rand*D)+1;        
    %ArandomlyChosenSolutionIsUsedInProducingaMutantSolutionOfTheSolution i
            neighbour=fix(rand*(FoodNumber))+1;       
            %RandomlySelectedSolutionMustBeDifferentFromTheSolution i        
            while(neighbour==ind1)
                neighbour=fix(rand*(FoodNumber))+1;
            end;        
            sol=Foods(ind1,:);
            %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
            sol(Param2Change)=Foods(ind1,Param2Change)+...
                (Foods(ind1,Param2Change)-Foods(neighbour,Param2Change))*...
                (rand-0.5)*2;        
       %IfGeneratedParameterValueIsOutOfBoundaries,ItIsShiftedOntoTheBounda
        ind=find(sol<x_low(:,1)');
        sol(ind)=x_low(ind,1)';
        ind=find(sol>x_high(:,1)');
        sol(ind)=x_high(ind,1)';        
        %Evaluate new solution
        ObjValSol=MGF(sol',H); 
        EvaluacioneFO=EvaluacioneFO+1;
        FitnessSol=calculateFitness(ObjValSol*10000);        
        %GreedySelectionIsAppliedBetweenTheCurrentSolution i and its mutant
        %If the mutant solution is better than the current solution i, 
        %replaceTheSolutionWithTheMutantAndResetTheTrialCounterOfSolution i
            if (FitnessSol>Fitness(ind1))
                Foods(ind1,:)=sol;
                Fitness(ind1)=FitnessSol;
                F_x_(ind1)=ObjValSol;
                trial(ind1)=0;
            else
                %IfSolution ind1 can'tBeImproved,IncreaseIts trial counter
                trial(ind1)=trial(ind1)+1; 
            end;
        end;    
        ind1=ind1+1;
        if (ind1==(FoodNumber)+1) 
            ind1=1;
        end;   
    end; 
    %/*The best food source is memorized*/
    ind=find(F_x_==min(F_x_));
    ind=ind(end);
    if (F_x_(ind)<GlobalMin)
         GlobalMin=F_x_(ind);
         GlobalParams=Foods(ind,:);
    end;       
    %%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DetermineThe food sources whose trial counter exceeds the "limit"value 
    %In Basic ABC, only one scout is allowed to occur in each cycle
    ind=find(trial==max(trial));
    ind=ind(end);
    if (trial(ind)>limit)
        trial(ind)=0;
        %sol=(ub-lb).*rand(1,D)+lb;
        for ind2=1:D 
            sol(1,ind2)=(x_low(ind2,1)+rand()*...
           (x_high(ind2,1)-x_low(ind2,1))); 
        end
        ObjValSol=MGF(sol',H); 
        EvaluacioneFO=EvaluacioneFO+1;
        FitnessSol=calculateFitness(ObjValSol*10000); 
        Foods(ind,:)=sol;
        Fitness(ind)=FitnessSol;
        F_x_(ind)=ObjValSol;
    end;
    %fprintf('Iter=%d F_x_=%g, Evaluaciones=%d\n',k,GlobalMin,EvaluacioneFO);
    k=k+1;
end % End of ABC
fprintf('Iter=%d F_x_=%g, Evaluaciones=%d\n',k,GlobalMin,EvaluacioneFO);
tiempo=toc;
%x_best=GlobalParams';
%GlobalMins(r)=GlobalMin;
grafica(GlobalParams',H,DB);
% temp=GlobalMin;
end

%% FUNCION GRAFICA GAUSSIANAS E IMAGEN SEGMENTADA: %%%%%%%%%%%%%%%%%%%%%%%%
% Recibe x_best: La mejor particula, D: dimensiones de cada particula
% H: histograma de la imagen, DB: imagen en escala de gris
function grafica(x_best,H,DB)
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
%     error3=HellingerDistance(H,mix)+(3/((medias(3)-medias(2))+(medias(2)-medias(1))+eps));
    error=HellingerDistance(H,mix);
end
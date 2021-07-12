% close all
clear all

% ############## Parameters ################
CommunicationTime = 5; %timestep in minute for realtime communication between VPP platform and batteries
ForecastTimeStepInit = 2; %Time Step in hour (between 1/2 and 12)
Days_Considered = 30;
Forecast_Inaccuracy = 1;  %1 if we ake into account a 25% max error in the forecasts, 0 if we consider perfect forecast
Percentage_Power_for_Critical_Commitment = 1;  %this value corresponds to the percetnage of total batteries power of an export commitment that requires to have all batteries with energy.
%An example is : if we have 10 batteries at max power of 10kW, then if we commit for energy that corresponds to 100kW during 30 minutes (50kWh), then we need all batteries to have enough energy so they all contribute to the required power.
% It is important to consider that power as if we only consider energy, 8 batteries might have enough energy, but they cannot export such power --> the commitment will not be met. 
Safety_Coefficient = 1; %this represents the percentage of the big export that we can commit to without doubt.
Type_of_Selling_Price = 1;  %1 for Market Day-ahead price, 2 for no selling tariff
Type_of_Buying_Price = 1; %1 for Agile tarif, 2 for ToU with 2 prices (Peak and low price), 3 for flat tariff

Flat_Buying_Price = 17/100/1000; %17p£/kWh
Peak_Buying_Price = 19/100/1000;
Low_Buying_Price = 12/100/1000;



load('All_Households_Data.mat');
% We remove house 44 as the owner has been removed from the project trial
Porigin(:,44)=[];
Dorigin(:,44)=[];

%We add the VPP (PV and Wind)
%Scaling of Wind corresponding to the approximation given by the owner from Findhorn
for k=1:size(Wind,1)
   if Wind(k,1)>5
       Wind(k,1)=5;
   end
end
PV_VPP =  Porigin(:,size(Porigin,2)-14);

Porigin = [Porigin, PV_VPP*20, Wind/max(Wind)*max(PV_VPP)*22];
Dorigin = [Dorigin, zeros(size(Dorigin,1),2)];

Households_Demand = Dorigin;
Households_Production = Porigin;
Import_Cost = zeros(size(Households_Demand+1,2),1);
Export_Revenue = zeros(size(Households_Demand+1,2),1);
Total_Bill_RT = zeros(size(Households_Demand+1,2),1);
Total_Bill_No_Battery = zeros(size(Households_Demand+1,2),1); 
Total_Bill_Heuristic = zeros(size(Households_Demand+1,2),1);
Total_Bill_HeuristicNoSP =  zeros(size(Households_Demand+1,2),1);
Total_Bill_HeuristicNoBatt = zeros(size(Households_Demand+1,2),1);
Expected_Bill_After_Optim = zeros(size(Households_Demand+1,2),1);  % This is the bill the agent is 
% expecting after he run his optimization, based on the forcast power
% demand & production, it does not include Real Time decisions
% Beware it accurate only if there is no error forecast. To get real Value run it without uncertainty and 1/2 hour forecast)

epsilon = 0.1;
N_Bat = size(Dorigin,2);%length(Names_Households);
Efficiency = 0.9; % Ec; %Battery Efficiency
EfficiencyC = Efficiency;
EfficiencyD = EfficiencyC;
ForecastTimeStep = 60*ForecastTimeStepInit; %30; %forecast for every 30 minutes
MarketTimeStep = 30; %Time Step in minutes %60*Timestep; %30; %market bid for every 30 minutes
Timestep = MarketTimeStep/60;
Window = 3*24;% Optimization for one day
Length = Window/Timestep; %number of data points
    CostCycle =0*7.4074e-05*(9000*1.2)/1000 ;
    CostCyce0=8.4074e-05*(9000*1.2)/1000;
    CostCycle =1*CostCyce0;
    NdataConsidered = 24/Timestep*Days_Considered; %round((size(Timestamp,1)-1)/60/24); %computation for the whole week
Start =0*24/Timestep;
kWhCost = 200; %Battery cost /kWh
% k=1;
Capacity_Bat = 10300;
PowerRating_Bat = 5300;
MBCRange =   N_Bat*Capacity_Bat; % 6300 13500; %Range of different battery capacity
Pbatmax = N_Bat*PowerRating_Bat*Timestep;% 2*3300*Timestep max battery Power = 3.3kW% 3*max(P); %W  We consider there is no limit of battery power
WHCRange = MBCRange'; % Range of Battery capacity in Watt-Hour (WH)
% Ec = 0.85;% Battery  charging efficiency (Ec)

    MaxBcap = MBCRange ;% Maximum battery capcity in WHr [Variable Parameter]
    MinBcap = 0.001*MaxBcap;% Minimum battery capacity at 80% DoD
%     IBC = MaxBcap*0.61; % Initial Battery Capacity (IBC)of 80%MacBcap
    IBC = MaxBcap*0.51; % Initial Battery Capacity (IBC)of 80%MacBcap
RealTimeDataTimeStep = 1; %1 min data time step for RT


t = (Timestamp(1):minutes(RealTimeDataTimeStep):Timestamp(size(Timestamp,1)))';
t2 = (Timestamp(1):minutes(ForecastTimeStep):Timestamp(size(Timestamp,1)))';
t3 = (Timestamp(1):minutes(MarketTimeStep):Timestamp(size(Timestamp,1)))';  

SumDorigin = sum(Dorigin,2);
SumPorigin = sum(Porigin,2);
Ndata =  size(Porigin,1);
j = 1;
Pforecast = zeros(size(t2));
Dforecast = zeros(size(t2));


for i =1:ForecastTimeStep/RealTimeDataTimeStep:size(Porigin,1)-1
    Pforecast(j)= sum(SumPorigin(i:i+ForecastTimeStep/RealTimeDataTimeStep-1))*(1+Forecast_Inaccuracy*(0.25-rand(1)/2));
    Dforecast(j)= sum(SumDorigin(i:i+ForecastTimeStep/RealTimeDataTimeStep-1))*(1+Forecast_Inaccuracy*(0.25-rand(1)/2));
    j = j+1;
end

Pbid = zeros(size(t3));
Dbid = zeros(size(t3));
j = 1;
Pbid(1) = Pforecast(j);
Dbid(1) = Dforecast(j);
for i = 2:size(Pbid,1)
    if i== (j)*ForecastTimeStep/MarketTimeStep+1
        Pbid(i) = Pforecast(j+1)*MarketTimeStep/ForecastTimeStep;
        Dbid(i) = Dforecast(j+1)*MarketTimeStep/ForecastTimeStep;
        j = j+1;
    else
     Pbid(i) =  Pbid(i-1);
     Dbid(i) =  Dbid(i-1);
    end     
end

Pref = Pbid;
Dref = Dbid;



%% Price Definition           
BPorigin = ones(size(Timestamp));
SPorigin = zeros(size(Timestamp));

j =1;

switch Type_of_Selling_Price 
    case 1 %Selling Price = Day ahead Auction Price from Nordpool
         for i = 1:size(Timestamp,1)
            SPorigin(i) = NordPoolSellingPrice(j);
           if  Timestamp(i)== t3(j+1)
               j=j+1;
           end

        end
%                     SP = interp1(t,SPorigin,t3);
        SP3 = interp1(t,SPorigin,t3);             

    case 2 %No selling Price
%                      SP = interp1(t,SPorigin,t3);
        SP3 = interp1(t,SPorigin,t3);                                 
end
j = 1;
switch Type_of_Buying_Price 
    case 1 %Buying Price = Agile ToU (Octopus)
         for i = 1:size(Timestamp,1)
            BPorigin(i) = BuyingPriceOctopus30min(j);
           if  Timestamp(i)== t3(j+1)
               j=j+1;
           end
        end
%                     BP = interp1(t,BPorigin,t3);
        BP3 = interp1(t,BPorigin,t3); 
    case 2
        HHP1 = 7.5;
        HHP2 = 22.5;
        T1 = t(1:HHP1*60+1);% Hours of  LHBP tariff 
        T2 = t(HHP1*60+2:HHP2*60+1);% Hours of HHBP tariff
        T3 = t(HHP2*60+2:1440);% Hours of LHBP tariff
        LHBP1 = Low_Buying_Price*ones(length(T1),1);% Low buying price for all the hours
        HHBP1 = Peak_Buying_Price*ones(length(T2),1);% High buying price for all the hours
        LHBP2 = Low_Buying_Price*ones(length(T3),1);% Low buying price for all the hours
        Daybuyprice = [LHBP1;HHBP1;LHBP2];
        BPorigin = repmat(Daybuyprice,30+1,1);
        BPorigin = [BPorigin; BPorigin(size(BPorigin,1))];
%                     BP = interp1(t,BPorigin,t3);
        BP3 = interp1(t,BPorigin,t3); 

    case 3 %Flat Buying Price
        BPorigin = Flat_Buying_Price*BPorigin;
%                     BP = interp1(t,BPorigin,t3);
        BP3 = interp1(t,BPorigin,t3); 
end






%% we cut all vectors to the size we need for the simulation          
t2 = t2(Start*MarketTimeStep/ForecastTimeStep+1:Start*MarketTimeStep/ForecastTimeStep+NdataConsidered*MarketTimeStep/ForecastTimeStep+1);
t3 = t3(Start+1:Start+NdataConsidered+1);
BP3 = BP3(Start+1:Start+NdataConsidered+1);
SP3 = SP3(Start+1:Start+NdataConsidered+1);

SystemBuyPriceMWh = SystemBuyPriceMWh(Start+1:Start+NdataConsidered+1);
NetImbalanceVolumeMWh = NetImbalanceVolumeMWh(Start+1:Start+NdataConsidered+1);
BP = BP3;
SP = SP3;
Pforecast = Pforecast(Start*MarketTimeStep/ForecastTimeStep+1:Start*MarketTimeStep/ForecastTimeStep+NdataConsidered*MarketTimeStep/ForecastTimeStep+1);
Dforecast = Dforecast(Start*MarketTimeStep/ForecastTimeStep+1:Start*MarketTimeStep/ForecastTimeStep+NdataConsidered*MarketTimeStep/ForecastTimeStep+1);


Pbid = Pbid(Start+1:Start+NdataConsidered+1);
Dbid= Dbid(Start+1:Start+NdataConsidered+1);
P = Pbid;
D = Dbid;           



SumPorigin = SumPorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
Porigin = Porigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1,:);
Dorigin = Dorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1,:);
% Dorigin= Dorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
SumDorigin= SumDorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
t = t(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
BPorigin = BPorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
SPorigin = SPorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
Porigininit = Porigin;
Dorigininit = Dorigin;

%% We start the optimization
% optimization constraints
Pmax = 3*max(max(P),max(D));
BPrice = BP(1:min(Length,size(BP,1))); %0.1* [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1.9 1.9 1.9 1.9 1 1 1 1 1]'; %0.9*ones(size(Demand));
SPrice = SP(1:min(Length,size(BP,1))); %0.1*[1 1 1 1 1 1 1 1 1 0.1 0.1 0.1 1 1 1 1 1 1 1 1 1 10000.9 10000.9 1.0]'; %0.01*ones(size(Demand));
    SoC = zeros(length(P),1);
Energysold = zeros(length(P),1);
Energybought = zeros(length(P),1);
    FromGrid = zeros(length(P),1);
    ToGrid = zeros(length(P),1);
%     SoC(1) = IBC;
    SoCmax = MBCRange ; %kWh 
    SoCmin = MinBcap;
    SoCInit = IBC; %kWh
    Powerbat = zeros(length(P),1);
    ellapsedtime = zeros(100,1);
    indice = 0;
    for i = 1:Length:length(P)
            indice = indice +1;
%         for j = 0:1:Length-1
            tic
            LengthOptim = min(i+Length,length(D))-i+1;  % Because the size of P might not be a multiple of the Length used for the optimization
            Demand = D(i:i+LengthOptim-1);
            Production = P(i:i+LengthOptim-1);
            BPrice = BP(i:i+LengthOptim-1);
            SPrice = SP(i:i+LengthOptim-1);
            f = [BPrice;-SPrice;zeros(LengthOptim,1); zeros(LengthOptim,1);zeros(LengthOptim,1);zeros(LengthOptim,1);zeros(LengthOptim,1); zeros(LengthOptim,1) ];

            Index_binary = size(f,1)-2*LengthOptim+1:size(f,1); %[121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144];
            Aeq = [eye(LengthOptim) zeros(LengthOptim) EfficiencyD*eye(LengthOptim) -eye(LengthOptim)/EfficiencyC -eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %Pgrin in = Demand - Production + Pbat_recharge/efficiency -Pbat_Discharge*Efficiency +X1   || X1 is here to make t possible to have Pgrid IN !=0 when Pgrid out =0
                   zeros(LengthOptim) eye(LengthOptim) -EfficiencyD*eye(LengthOptim) eye(LengthOptim)/EfficiencyC zeros(LengthOptim) -eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim)]; %Pgrin out = - Demand + Production - Pbat_recharge/efficiency + Pbat_Discharge*Efficiency +X2 || X2 is here to make t possible to have Pgrid OUT !=0 when Pgrid in =0
            beq = [ Demand - Production ; Production - Demand ];

            A = [eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) -Pmax*eye(LengthOptim) zeros(LengthOptim); %makes sure Pgrid in = 0 when Pgrid out != 0
                 zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) Pmax*eye(LengthOptim) zeros(LengthOptim);  %makes sure Pgrid out = 0 when Pgrid in != 0
                 zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) -Pmax*eye(LengthOptim) zeros(LengthOptim); %makes sure X2  = 0 when X1 != 0
                 zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) Pmax*eye(LengthOptim) zeros(LengthOptim); %makes sure X1  = 0 when X2 != 0
                 zeros(LengthOptim) zeros(LengthOptim) -tril(ones(LengthOptim)) tril(ones(LengthOptim)) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure SoC <= SoC max
                 zeros(LengthOptim) zeros(LengthOptim) tril(ones(LengthOptim)) -tril(ones(LengthOptim)) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure SoC >= SoCmin
                 zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) -Pbatmax*eye(LengthOptim) ; %makes sure Pbatdischarge  = 0 when Pcharge != 0
                 zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) Pbatmax*eye(LengthOptim)]; %makes sure Pbat_charge  = 0 when Pdischarge != 0

            b = [zeros(LengthOptim,1); Pmax*ones(LengthOptim,1); zeros(LengthOptim,1); Pmax*ones(LengthOptim,1); (SoCmax-SoCInit)*ones(LengthOptim,1); (SoCInit - SoCmin)*ones(LengthOptim,1); zeros(LengthOptim,1); Pbatmax*ones(LengthOptim,1)]; 
            lb = [zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1);  zeros(LengthOptim,1);  zeros(LengthOptim,1); zeros(LengthOptim,1)]  ;
            ub = [Pmax*ones(LengthOptim,1); Pmax*ones(LengthOptim,1); Pbatmax*ones(LengthOptim,1);  Pbatmax*ones(LengthOptim,1); Pmax*ones(LengthOptim,1);  Pmax*ones(LengthOptim,1);  ones(LengthOptim,1);  ones(LengthOptim,1) ]  ;

            %if we wanted to change the battery power for discharge for the night (to
            %follow a simulation where price was very high for selling and see if the simulation is able to find a good behaviour)
%             ub(LengthOptim*2+21:LengthOptim*2+24)=3;
            % ub(LengthOptim*3+21:LengthOptim*3+24)=3;
%             x0=[];
%              options=optimoptions('intlinprog','Display','off');
%             options.Display = 'off';
            [x,fval1,exitflag1,output1]  = intlinprog(f, Index_binary, A,b, Aeq, beq, lb, ub);
            elapsedTime = toc ;
%             ellapsedtime(indice) = cputime - t;
            ellapsedtime(indice) = elapsedTime;

            SoC(i:i+LengthOptim-1) = -tril(ones(LengthOptim))*x(LengthOptim*2+1:LengthOptim*3)+tril(ones(LengthOptim))*x(LengthOptim*3+1:LengthOptim*4)+SoCInit*ones(LengthOptim,1);
            PgridIN = x(1:LengthOptim);
            PgridOUT = x(LengthOptim+1:LengthOptim*2);
            Pbat = x(LengthOptim*2+1:LengthOptim*3)-x(LengthOptim*3+1:LengthOptim*4);
            X1 = x(LengthOptim*4+1:LengthOptim*5);
            X2 = x(LengthOptim*5+1:LengthOptim*6);
            Alpha = x(LengthOptim*6+1:LengthOptim*7);
            Beta = x(LengthOptim*7+1:LengthOptim*8);

            SoCInit = SoC(i+LengthOptim-1);
            for j = 0:LengthOptim-1
                Energybought(i+j) = PgridIN(1+j)*BPrice(j+1);
                Energysold(i+j) =PgridOUT(1+j)*SPrice(j+1);          
                FromGrid(i+j) = max(0,(D(i+j)-P(i+j)-max(0,Pbat(j+1))*EfficiencyD+max(0,-Pbat(j+1)/EfficiencyC)));
                ToGrid(i+j) = max(0, (P(i+j)-D(i+j)+max(0,Pbat(j+1)*EfficiencyD)-max(0,-Pbat(j+1)/EfficiencyC)));
            end
            Powerbat (i:i+LengthOptim-1) = Pbat;

    end



    %
% Market Bid
    SoCmax = MBCRange ; %kWh 
    SoCmin = MinBcap;
Pbatmax = N_Bat*PowerRating_Bat*Timestep;% 2*3300*Timestep max battery Power = 3.3kW% 3*max(P); %W  We consider there is no limit of battery power


    %we will now determine which mode the battery is in.
    % We list 4 modes:
    %1. Battery follows P-D
    %    1.a: SoC increases (mode 11)
    %    1.b: SoC decreases (mode 12)
    % 2. Battery must ensure a SoC at the end of the time interval
    %    2.a: SoC must be equal or higher than the one obtained during the (mode 21)
    %    optimization at the end of the considered time interval --> use
    %    the algo to reach the SoC
    %    2.b: SoC must be equal or lower than the one obtained during the  (mode 22)
    %    optimization at the end of the considered time interval --> use
    %    the algo to reach the SoC


    % Note that sometimes, the optimization does not follow a logic choice
    % and stops recharging to recharge later as it gives the same result in
    % terms of revenue. The algorithm does identify this


   %Determine the time intervals (everytime there is a new price or that the battery does not  

  % first, we determine the SoC targets (every time a SoC evolution changes)  
   k=0;
InitialMode =0;
while InitialMode==0
    k=k+1;
    if SoC(k) - SoC(k+1) <-epsilon
        InitialMode =1; %1 for Charging Mode 

    elseif SoC(k) - SoC(k+1) > epsilon
       InitialMode =-1; %-1 for discharging Mode 
    end
end
    j=1;
SoCTargets =  zeros(size(t3));
SoCTargetsType =  zeros(size(t3)); %0 if target does not correspond to a change of price, 1 otherwise
PreviousMode = InitialMode;
last_index = 1;
ChargeMode = zeros(size(P));
ChargeMode(1) = InitialMode;
CurrentMode = InitialMode;
for k = 1:1:size(P,1)-1
    if SoC(k) - SoC(k+1) <-epsilon
        CurrentMode =1; %1 for Charging Mode 
    elseif SoC(k) - SoC(k+1) > epsilon
       CurrentMode =-1; %-1 for discharging Mode 
    end
    if CurrentMode*PreviousMode <0  || abs(BP(k)-BP(k+1))>0 ||  abs(SP(k+1)-SP(k))>0 %5|| (t3(k).Hour== t3(startfrr).Hour) || t3(k).Hour == t3(endfrr).Hour
       SoCTargets(last_index:k)=SoC(k);
        last_index = k+1;
    end
    PreviousMode = CurrentMode;
    ChargeMode(k+1)=PreviousMode;

end
% Determine the times of change of time interval:
Counter = 0; %counter that indicates the number of time interval 
for k = 1:1:size(P,1)-1

    if abs(SoCTargets(k)- SoCTargets(k+1))>epsilon || abs(BP(k)-BP(k+1))>0 ||  abs(SP(k+1)-SP(k))>0 %|| (t3(k).Hour== t3(startfrr).Hour) || t3(k).Hour == t3(endfrr).Hour
       Counter = Counter+1;
    end
end

TimeInterval = datetime(zeros(Counter+2,1),0,0,0,0,0);
IndexTimeInterval = ones(Counter+2,1);
j = 1;
TimeInterval(j)=t3(1);
IndexTimeInterval(j)=1;
j = j+1;
for k = 2:1:size(P,1)-1
    if abs(SoCTargets(k)- SoCTargets(k-1))>epsilon || abs(BP(k-1)-BP(k))>0 ||  abs(SP(k-1)-SP(k))>0 %|| (t3(k).Hour== t3(startfrr).Hour) || t3(k).Hour == t3(endfrr).Hour
       TimeInterval(j)=t3(k);
       IndexTimeInterval(j)=k;
       j = j+1;
    elseif abs(SoCTargets(k)- SoCTargets(k-1))<epsilon/10  && k>3%if there was no change in the target during this timeslot, we assume that the target is to be met for the next timeslot
        TimeInterval(j-1)=t3(k);
        IndexTimeInterval(j-1)=k;
    end
end

TimeInterval(j)=t3(k+1);
IndexTimeInterval(j)=k+1;
% Determine the Mode within each time interval
j = 1;
CurrentModeVector = zeros(size(P));
for k = 1:size(P,1)-1
    if t3(k) > TimeInterval(j) || k==1
        %Notice that P and D are energy vectors. ot power... Hence we do
        %not need to include the time step. However, Efficiency must be
        %added.
if k>=3*48+14*2
    tempre=0;
end
        if (ChargeMode(k) ==1  && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1))))*Efficiency - ...
                   abs(SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j))))<= epsilon) || (ChargeMode(k) ==1 && abs(SoCTargets(IndexTimeInterval(j+1)) - SoCmax)<epsilon...
                   &&  sum(FromGrid(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))<epsilon    )% || ( ChargeMode(k) ==1 && t3(k).Hour>= t3(startfrr).Hour && t3(k).Hour <= t3(endfrr).Hour)%abs(SoCTargets(IndexTimeInterval(j+1)-1) - SoCTargets(IndexTimeInterval(j+1)) )<epsilon) %charging mode && charge by P -D (mode 11) or if the battery was charged by P until it was full
        CurrentModeVector(k) = 11;
        elseif (ChargeMode(k) ==-1 && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))/Efficiency - ...
                   (SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j)))))<= epsilon) || (ChargeMode(k) ==-1 &&  abs(SoCTargets(IndexTimeInterval(j+1)) - SoCmin)<epsilon...
                          &&  sum(ToGrid(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))<epsilon ) %|| (ChargeMode(k) ==-1 && t3(k).Hour>= t3(startfrr).Hour && t3(k).Hour <= t3(endfrr).Hour) %abs(SoCTargets(IndexTimeInterval(j+1)-1) - SoCTargets(IndexTimeInterval(j+1)) )<epsilon) %discharge mode
        CurrentModeVector(k) = 12;    
        elseif ChargeMode(k) ==1  && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1))))*Efficiency - ...
                   abs(SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j))))> epsilon  %charging mode && charge by P -D (mode 11)
        CurrentModeVector(k) = 21;
        elseif ChargeMode(k) ==-1 && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))/Efficiency - ...
                   (SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j)))))> epsilon  %discharge mode
        CurrentModeVector(k) = 22;    
        end
        j = j+1;        
    else
       CurrentModeVector(k) = CurrentModeVector(k-1); 
    end
end

% We now compute the bid curves
% t3, Pbid, Dbid are the vectors to be considered

%We determine SoCTargetBids, same as SoCTarget, but with the time scale of
%the Market
SoCTargetBid = zeros(size(t3));
ChargeModeBid = zeros(size(t3));
CurrentModeVectorbid = zeros(size(t3));

SoCTargetBid(1) = SoCTargets(1);
ChargeModeBid(1) = ChargeMode(1);
CurrentModeVectorbid(1) = CurrentModeVector(1);
LengthTimeIntervalBid = zeros(size(t3));
j=1;
previousindex = 1;
for k = 1:size(t3,1)-1
if k>=3*48+11.5*2
    tempre=0;
end
    if t3(k) == t3(IndexTimeInterval(j))
        SoCTargetBid(k) = SoCTargets(IndexTimeInterval(j+1));
        ChargeModeBid(k) = ChargeMode(IndexTimeInterval(j+1)); 
        LengthTimeIntervalBid(previousindex:k-1)=k-previousindex;
        CurrentModeVectorbid(k) = CurrentModeVector(IndexTimeInterval(j+1));
        previousindex = k;
        j = j+1;
    elseif k==1
        SoCTargetBid(1) = SoCTargets(1);
        ChargeModeBid(1) = ChargeMode(1);
        CurrentModeVectorbid(1) = CurrentModeVector(1);
    else
        SoCTargetBid(k) = SoCTargetBid(k-1);
        ChargeModeBid(k) = ChargeModeBid(k-1);
        CurrentModeVectorbid(k) = CurrentModeVectorbid(k-1);
    end
end
SoCTargetBid(k+1)= SoCTargets(size(SoCTargets,1));
SoCTargetBid = [SoCTargetBid(1);SoCTargetBid];

LengthTimeIntervalBid(previousindex:k+1) = k-previousindex;
%We now compute the forecast of the SoC for the battery 
SoCbid = zeros(size(t3));
Ebatbid = zeros(size(t3));
j=1;

SoCbid(1) = SoC(1);
%Initialization
       RealisedOutputprevious = 0;   
SoCprevious = SoC(1);
timeinterval = LengthTimeIntervalBid(1);
LastSoCtarget = SoCprevious;
NextSoCtarget = SoCTargetBid(1);
        if CurrentModeVectorbid(1)==21 || CurrentModeVectorbid(1)==22
            soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
        else 
            soctargettype = 0;
        end
          StartTime = 0; % time (indice) at which we start this market timeslot 
        Et=0;  %energy that is exported from the battery at the time t (j)
        slopeprevious = 0; %slope of the battery power at the previous time step k
       time1previous = 0;  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
       RealisedOutput = 0;  
       NeededEnergy = SoCTargetBid(1)  - SoCprevious; 
CurrentModeVectorbid = [CurrentModeVectorbid(1);CurrentModeVectorbid];
CurrentModeVector = [CurrentModeVector;CurrentModeVector(size(CurrentModeVector,1))];
IndexTimeInterval = [IndexTimeInterval; IndexTimeInterval(size(IndexTimeInterval,1))];
for k = 1:size(t3,1)

if k>=3*24*2+3.5*2+1
    tempre=0;
end

    if t3(k) == t3(IndexTimeInterval(j)+1) %if we are in a new timestep
        StartTime = k-1; % time (indice) at which we start this market timeslot 
        Et=0;  %energy that is exported from the battery at the time t (j)
        slopeprevious = 0; %slope of the battery power at the previous time step k
       time1previous = 0;  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
       RealisedOutput = 0;
       RealisedOutputprevious = 0;        
        timeinterval = LengthTimeIntervalBid(k);
        LastSoCtarget = SoCprevious;
        NextSoCtarget = SoCTargetBid(k);
         NeededEnergy = NextSoCtarget - LastSoCtarget; 

        if CurrentModeVector(IndexTimeInterval(j)+1)==21 || CurrentModeVector(IndexTimeInterval(j)+1)==22% ||  (t3(k).Hour>= t3(startfrr).Hour && t3(k).Hour <= t3(endfrr).Hour)
           soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
        else 
            soctargettype = 0;
        end
        j = j+1;

    end

    CurrentTime = k-1-StartTime; %the time from which we recompute the energy that we still need
     time1 = 2*(sign(NeededEnergy)*Pbatmax*Efficiency*(timeinterval-CurrentTime/2+1/2)-Et*(CurrentTime/2+1/2)+RealisedOutput-sign(NeededEnergy)*min(abs(NeededEnergy),Pbatmax*Efficiency*timeinterval))/(sign(NeededEnergy)*Pbatmax*Efficiency-Et+0.000000001); %determines the the time at which we need to use the battery at its maximum power
    if time1 > timeinterval %if we do not need to be at full power within this time interval
        slope = 2*(sign(NeededEnergy)*min(abs(NeededEnergy),Pbatmax*Efficiency*timeinterval)-RealisedOutput - Et*(timeinterval-CurrentTime))/((timeinterval-CurrentTime)*(timeinterval-CurrentTime+1));
    else
        slope = (sign(NeededEnergy)*Pbatmax*Efficiency-Et)^2/(sign(NeededEnergy)*Pbatmax*Efficiency*(2*timeinterval+1-2*CurrentTime)-Et+2*(RealisedOutput-sign(NeededEnergy)*min(abs(NeededEnergy),Pbatmax*Efficiency*timeinterval))+0.000001);
    end 
    if abs(sign(NeededEnergy)*Pbatmax*Efficiency-Et ) < epsilon  %if we already reached the maximum power at previous step while we were still ramping (happens when time1 is an integer
        time1 = time1previous; %we keep the same value as previous iteration, as it will still be the maximum power for the battery
        slope=slopeprevious;
    end

    Ebatbid(k) =slope+Et;% Ebatbid(k)+slope; %The battery power is equal to what should be produced to reach the bid minus what is already produced by the PV
    Et = Ebatbid(k);
    s = sign(Pbid(k)-Dbid(k));
    if sign(NeededEnergy)*s>0
       Ebatbid(k) = sign(Et)*max(abs(Ebatbid(k)),abs( Pbid(k)-Dbid(k))*((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency) );
    end
    Ebatbid(k)= min(Pbatmax, max(-Pbatmax,Ebatbid(k)));
%     if CurrentModeVectorbid(k)==11 || CurrentModeVectorbid(k)==12
    if soctargettype ==0
        Ebatbid(k) = max(-Pbatmax,min(Pbatmax,(Pbid(k)-Dbid(k))*((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency))) ;

    end

    if k-StartTime - time1 >= -epsilon  && time1 >0 % if we have already reached the time at which the battery should provide all its energy, then we must set its power to Pmax
       Ebatbid(k)=sign(NeededEnergy)*Pbatmax;
    end   
    SoCbid(k) = SoCprevious+Ebatbid(k);
    [Ebatbid(k),SoCbid(k)] = coercsoc(SoCbid(k),SoCmax,SoCmin,Ebatbid(k),1,LastSoCtarget, NextSoCtarget,soctargettype,Pbatmax);  %1 because we consider energy and not power.
    SoCbid(k) = SoCprevious+Ebatbid(k);
    RealisedOutput = RealisedOutputprevious + Ebatbid(k);

%     if RealisedOutput > NeededEnergy %if we already exported too much energy compared to our bid, we need to reduce the amount of eergy exported
    if sign( NeededEnergy - RealisedOutput)*sign(NeededEnergy) <0  %if there is more energy provided than required %if we already exported too much energy compared to our bid, we need to reduce the amount of eergy exported
       Ebatbid(k) = Ebatbid(k) + NeededEnergy - RealisedOutput ;%max(-Pbatmax, min(Pbatmax,NeededEnergy-RealisedOutput)); %so we charge the battery with the surplus of energy. making sure the battery's max power is below Pmax
    end    
    SoCbid(k) = SoCprevious+Ebatbid(k);
    RealisedOutput = RealisedOutputprevious + Ebatbid(k);
    RealisedOutputprevious = RealisedOutput;
    time1previous = time1;
    slopeprevious = slope;    
    SoCprevious = SoCbid(k);
%     Et = Ebatbid(k);


end

Extra_Power_Needed = zeros(size(t3));
To_Grid_bid = zeros(size(t3));
From_Grid_bid = zeros(size(t3));

for k = 1:size(t3,1)

    s = -sign(Ebatbid(k));
    a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);


    From_Grid_bid(k) = max(0,Dbid(k)-Pbid(k)+ Ebatbid(k)*a);
    To_Grid_bid(k) = max(0,Pbid(k)-Dbid(k)- Ebatbid(k)*a);
%Now we detect when there is a need for big power
%(= when the requested power = MaxPower*Percentage_Power_for_Critical_Commitment)

    if -Ebatbid(k)>= Pbatmax*Percentage_Power_for_Critical_Commitment
        Extra_Power_Needed(k)=1;
        To_Grid_bid(k) = To_Grid_bid(k)*Safety_Coefficient;
    end       

end


% We now compute the real Battery operation curves
% t, Porigin, Dorigin are the vectors to be considered

Pbatmax = PowerRating_Bat*RealTimeDataTimeStep/60*1;% max battery Power for single battery
SoCmax= SoCmax/N_Bat;
SoCmin = SoCmin/N_Bat;

%We determine SoCTargetreals, same as SoCTarget, but with the time scale of
%the Market
SoCTargetBid = SoCTargetBid/N_Bat;  %We determine the SoC target for all batteries (assumed identical here)
SoCbidscaled = SoCbid/N_Bat;
SoCTargetreal = zeros(size(t));
ChargeModereal = zeros(size(t));
CurrentModeVectorreal = zeros(size(t));
SoCTargets = SoCTargets/N_Bat;
SoCTargetreal(1) = SoCTargetBid(1);
ChargeModereal(1) = ChargeMode(1);
CurrentModeVectorreal(1) = CurrentModeVector(1);
LengthTimeIntervalreal = zeros(size(t));
j=1;
previousindex = 1;
SoCTargetrealprevious = SoCTargetreal(1);
ChargeModerealprevious = ChargeModereal(1);
CurrentModeVectorrealprevious = CurrentModeVectorreal(1);

for k = 1:size(t,1)-1
    if k >=24*60+11.5*60
        tempb=0;
    end
    if t(k) == t3(IndexTimeInterval(j))
%         SoCTargetreal(k) = SoCTargetBid(IndexTimeInterval(j+1));
        SoCTargetreal(k) = SoCbidscaled(IndexTimeInterval(j+1));
        ChargeModereal(k) = ChargeMode(IndexTimeInterval(j+1)); 
        LengthTimeIntervalreal(previousindex:k-1)=k-previousindex;
        CurrentModeVectorreal(k) = CurrentModeVector(IndexTimeInterval(j+1));
        previousindex = k;
        SoCTargetrealprevious = SoCTargetreal(k);
        ChargeModerealprevious = ChargeModereal(k);
        CurrentModeVectorrealprevious = CurrentModeVectorreal(k);
        j = j+1;
    else
        SoCTargetreal(k) = SoCTargetrealprevious; 
        ChargeModereal(k) = ChargeModerealprevious; %ChargeModereal(k-1);
        CurrentModeVectorreal(k) = CurrentModeVectorrealprevious;%CurrentModeVectorreal(k-1);
    end
end
SoCTargetreal(k+1)= SoCTargetBid(size(SoCTargetBid,1));

LengthTimeIntervalreal(previousindex:k+1) = k-previousindex;
%We now compute the forecast of the SoC for the battery 
SoCreal = zeros(size(t,1),N_Bat);
Ebatreal = zeros(size(t,1),N_Bat);


%% Real Time Operation

%Initialization for the SoC timeframe
RealisedOutputprevious = zeros(1,N_Bat);   
SoCprevious = SoC(1)*ones(1,N_Bat)/N_Bat;
timeinterval = LengthTimeIntervalreal(1);
LastSoCtarget = SoCprevious;
NextSoCtarget = SoCTargetreal(1);
if CurrentModeVectorreal(1)==21 || CurrentModeVectorreal(1)==22
    soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
else 
    soctargettype = 0;
end
StartTime = 0; % time (indice) at which we start this market timeslot 
slopeprevious = zeros(1,N_Bat); %slope of the battery power at the previous time step k
time1previous = zeros(1,N_Bat);  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
RealisedOutput = zeros(1,N_Bat);  
NeededEnergy = SoCTargetBid(1)*ones(1,N_Bat)  - SoCprevious; 



%Initialization for the market time frame
RealisedOutputpreviousmarket = 0;   
% SoCpreviousmarket = SoC(1)*ones(1,N_Bat);
timeintervalmarket = MarketTimeStep;
% LastSoCtargetmarket = SoCpreviousmarket;
% NextSoCtargetmarket = SoCTargetBid(1);
% if CurrentModeVectorbid(1)==21 || CurrentModeVectorbid(1)==22
%     soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
% else 
%     soctargettype = 0;
% end
StartTimemarket = 0; % time (indice) at which we start this market timeslot 
Etmarket=zeros(1,N_Bat);  %energy that is exported from the battery at the time t (j)
slopepreviousmarket = 0; %slope of the battery power at the previous time step k
time1previousmarket = 0;  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
RealisedOutputmarket = 0;  
NeededOutput = To_Grid_bid(1);
% SoCprevious = SoCpreviousmarket;%SoC(1);
Required_Export_bid = max(0,To_Grid_bid-From_Grid_bid);
RealisedOutputMarketVector = zeros(size(t));
j=1;
indexbid = 1;
indexVector = zeros(size(t3));
TargetsBid = zeros(size(t3));

%to avoid the error for the last output
To_Grid_bid = [To_Grid_bid;0];
From_Grid_bid = [From_Grid_bid;0];
boolEt = 0;
Delta_Forecast = max(0,Dbid(indexbid)-Pbid(indexbid));
NeededOutputref = NeededOutput;
NeededOutput = NeededOutput*ones(1,N_Bat);
time1 = zeros(1,N_Bat);
ToGridreal = zeros(1,N_Bat);
FromGridreal = zeros(1,N_Bat);
Et= zeros(1,N_Bat); %energy that is exported from the battery at the time t (j)

EbatrequiredMarket = zeros(1,N_Bat);
RealisedOutputmarketindiv   = zeros(1,N_Bat);
% RealisedOutputpreviousmarketindiv   = zeros(1,N_Bat);
ToGridTotal = zeros(size(t));
FromGridTotal = zeros(size(t));
GridPowerVector = zeros(size(t));
CounterRebound = 0;
boolcounter = 0;
AdditionalEnergy = zeros(1,N_Bat);
PercentExtraprevious = zeros(1,N_Bat);
PercentDistribprevious = 1/N_Bat*ones(1,N_Bat);
VarDebug = zeros(size(t));
        previousredistrib = 0;
   ErequiredMarket = 0;     
redistributedEnergyprevious = 0;
CounterRealTime = -1;


                GridPowerOthers = zeros(1,N_Bat);
                GridPower = zeros(1,N_Bat);
                D_minus_P_aggregated = zeros(1,N_Bat);
                AvailableEnergy_others=zeros(1,N_Bat);
                PossibleAcceptanceofEnergy_others = zeros(1,N_Bat) ;
                RealisedOutputmarket = zeros(1,N_Bat) ;
                AvailableEnergy = zeros(1,N_Bat);
                PossibleAcceptanceofEnergy = zeros(1,N_Bat);
                RealGridPower = zeros(1,N_Bat);
                Percentavail = 1/N_Bat*ones(1,N_Bat);
                Percentacceptance = 1/N_Bat*ones(1,N_Bat);
                RealisedOutputMarket_Real_perBattery = zeros(1,N_Bat);
                Real_RealisedOutputpreviousmarket = 0;
                D_minus_P_houses = zeros(1,N_Bat);
                RealisedOutputMarketVector2 = zeros(size(t));
%                 DistributedEnergyPrevious = 0;
                EnergieAvailable = 0; %SumPorigin(1)-SumDorigin(1);
                EnergieSurplus = 0;
                EnergieSurplusprevious = 0;
                %we need to make sure production and consumption data
                %indexes correspond to the indexes of the SoC targets
 Porigin=Porigininit;               
 Porigin(1:30,:)=[]; 
 Porigin = [Porigin;zeros(30,N_Bat)];

 Dorigin=Dorigininit;               
 Dorigin(1:30,:)=[]; 
 Dorigin = [Dorigin;zeros(30,N_Bat)];  
 PowertargettypeforNextTimeSetp = 0;
%                 N_Bat=N_Bat;
%                 epsilon = epsilon;

CompA = zeros(size(t));
CompB = zeros(size(t));
Household_Power = zeros(size(t,1),N_Bat);
for k = 1:size(t,1)
    CounterRealTime = CounterRealTime+1;
  
   if k>1 
    if mod(CounterRealTime,CommunicationTime)==0  %if we have data from the VPP, we must update the current situation
        % Need to have NeededOutput(d) for d = 1:N_Bat
        % Need to have RealisedOutputmarket(d)
%             SoCrealpreviousall = SoCreal(k-1,:);
            EnergieAvailable = 0;  %corresponds to the sum of all the extra energy (due to SoC already full or too much power) that the households export as they cannot use it locally
            percentage_Allocate_EnergieAvailable = ones(1,N_Bat);
            for d=1:N_Bat
                if Porigin(k-1,d)-Dorigin(k-1,d)>= SoCmax-SoCreal(k-1,d)
                   EnergieAvailable = EnergieAvailable+ Porigin(k-1,d)-Dorigin(k-1,d)- (SoCmax-SoCreal(k-1,d));
                   percentage_Allocate_EnergieAvailable(d)=0;
                elseif Porigin(k-1,d)-Dorigin(k-1,d)>= Pbatmax
                   EnergieAvailable = EnergieAvailable+ Porigin(k-1,d)-Dorigin(k-1,d)- Pbatmax;
                   percentage_Allocate_EnergieAvailable(d)=0;

                end

            end
            %allocation of energy available based on current use of
            %battery, instead of SoC, but could be changed
            percentage_Allocate_EnergieAvailable = (Pbatmax-(Porigin(k-1,1:N_Bat)-Dorigin(k-1,1:N_Bat))).*percentage_Allocate_EnergieAvailable/...
                (0.0001+sum(Pbatmax*percentage_Allocate_EnergieAvailable-(Porigin(k-1,1:N_Bat)-Dorigin(k-1,1:N_Bat)).*percentage_Allocate_EnergieAvailable));


                Real_RealisedOutputmarkettemp = 0;
                for d=1:N_Bat
                  Real_RealisedOutputmarkettemp = Real_RealisedOutputmarkettemp + RealisedOutputMarket_Real_perBattery(d);
                end
                Real_RealisedOutputmarket = Real_RealisedOutputpreviousmarket + Real_RealisedOutputmarkettemp;
                RealisedOutputMarket_Real_perBattery = zeros(1,N_Bat);
                Real_RealisedOutputpreviousmarket = Real_RealisedOutputmarket;
                for d=1:N_Bat
                  RealisedOutputmarket(d) = Real_RealisedOutputmarket;
%                       if ( NeededOutput(d) - RealisedOutputmarket(d)) < -epsilon && soctargettypeMarket ==0 %if there is more energy provided than required %if we already exported too much energy compared to our real, we need to reduce the amount of eergy exported
%                            Etmarket(d) = 0; %Ebatreal(k,d) = Ebatreal(k,d) +  (RealisedOutputmarket - NeededOutput)/a;%max(-Pbatmax, min(Pbatmax,NeededOutput-RealisedOutput)); %so we charge the battery with the surplus of energy. making sure the battery's max power is below Pmax
%                       end  
                end

                for d=1:N_Bat
%                     if (timeintervalmarket-CurrentTimemarket)< 7 ||  ( NeededOutputref - Real_RealisedOutputmarket) < -epsilon%###########  6 min before the end of market time interval = Time at which we consider that Neededoutput = neededoutputref
%                     NeededOutput(d) = NeededOutputref;
%                     elseif NeededOutput(d) > NeededOutputref + epsilon
%                      NeededOutput(d) = max(NeededOutputref, NeededOutput_Real(d) - Sum_D_minus_P_Comm);
%                     end
                end


               NeededOutput_Real=NeededOutput;
               EnergieAvailable = EnergieAvailable+max(0,Real_RealisedOutputmarket-NeededOutputref)/(timeintervalmarket-CurrentTimemarket); %##################


               Sum_D_minus_P_Comm = 0;

              GridPowerReal = 0;


    else  %if we do not have the data from the VPP 
                for d=1:N_Bat
                    s = -sign(Ebatreal(k-1,d));
                    a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
%                     GridPower(d) = Porigin(k-1,d)-Dorigin(k-1,d)- Ebatreal(k-1,d)*a+GridPowerOthers(d);  % we assume during the time we do not have feedback from the VPP, other households export/import the same thing
 % TO BE CHECKED what gives best results . The choice below affects the
 % capability of  exchanging power from one battery to another.
%                     GridPower(d) = Porigin(k-1,d)-Dorigin(k-1,d)- Ebatreal(k-1,d)*a;% 0; %Porigin(k-1,d)-Dorigin(k-1,d)- Ebatreal(k-1,d)*a+GridPowerOthers(d);  % we assume during the time we do not have feedback from the VPP, other households export/import the same thing
                %we could change that by taking into account the forecast. To be
                %determined.
                end

                for d=1:N_Bat
                    if (timeintervalmarket-CurrentTimemarket)< 7   %###########  6 min before the end of market time interval = Time at which we consider that Neededoutput = neededoutputref
                    NeededOutput(d) = NeededOutputref;
                    elseif NeededOutput(d) > NeededOutputref + epsilon
                    end
                end


    end
   end

% there are 2 time lines: the ones from the optimization with the SoC
% changes, and the one from the market where we have to commit to the bids within the
% market time intervals
    if t(k) == t3(IndexTimeInterval(j)) %if we are in a new timestep
        PowertargettypeforNextTimeSetp = 0;
        StartTime = k-1; % time (indice) at which we start this market timeslot 
        Et=zeros(1,N_Bat);  %energy that is exported from the battery at the time t (j)
        slopeprevious = zeros(1,N_Bat); %slope of the battery power at the previous time step k
       time1previous = zeros(1,N_Bat);  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
       RealisedOutput = zeros(1,N_Bat);
       RealisedOutputprevious = zeros(1,N_Bat);        
        timeinterval = LengthTimeIntervalreal(k);
        LastSoCtarget = SoCprevious;
        NextSoCtarget = SoCTargetreal(k);
        if k>1      
         previousSoCtarget = SoCTargetreal(k-1);
        else 
            previousSoCtarget = SoCreal(1,1);
        end
         NeededEnergy = NextSoCtarget*ones(1,N_Bat) - LastSoCtarget; %SoCTargetreal(k)  - SoCTargetreal(k-1); 

        if CurrentModeVectorreal(k)==21 || CurrentModeVectorreal(k)==22 
            soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
            soctargettypeMarket = 1;

        else 
            soctargettype = 0;
        end

        for idx = 0:LengthTimeIntervalBid(j+2)-1  %if the next target will be requiringfull power
            if j+2<size(IndexTimeInterval,1) && (IndexTimeInterval(j+2)+idx)<=size(Extra_Power_Needed,1) && ...
                    Extra_Power_Needed(IndexTimeInterval(j+2)+idx) ==1 && Powertargettypeprevious ==0 
                PowertargettypeforNextTimeSetp = 2;
            end
        end
        Powertargettypeprevious = PowertargettypeforNextTimeSetp;
        j = j+1;
    end

    if t(k) == t3(indexbid)
        NeededOutput = max(0,To_Grid_bid(indexbid+1)-From_Grid_bid(indexbid+1));
        Delta_Forecast = 0;
%         GridExtraPower = 0;
        previousredistrib = 0;
        NeededOutputref = NeededOutput;
        if NeededOutput >epsilon*N_Bat/10
            Delta_Forecast = max(0,Dbid(indexbid+1)-Pbid(indexbid+1));
%             NeededOutput=NeededOutput+Delta_Forecast;
        end
        TargetsBid(indexbid+1) = NeededOutput;
        if NeededOutputref> epsilon*N_Bat/10 || PowertargettypeforNextTimeSetp==2
            Delta_for_Power_Equilibrium = 0;
            if PowertargettypeforNextTimeSetp==2
                 for d = 1:N_Bat
                   Delta_for_Power_Equilibrium= Delta_for_Power_Equilibrium+max(0,NextSoCtarget -  SoCreal(k-1,d)); %we compute the energy we need to fill batteries that need to be filled
                 end
                 NeededOutput = NeededOutput + Delta_for_Power_Equilibrium;%*Percentage_Power_for_Critical_Commitment;
            end                               
        end
        NeededOutput = NeededOutput*ones(1,N_Bat);
        NeededOutput_Real = NeededOutput;
        StartTimemarket = k-1; % time (indice) at which we start this market timeslot 
        Etmarket=zeros(1,N_Bat);  %energy that is exported from the battery at the time t (j)
        slopepreviousmarket = 0; %slope of the battery power at the previous time step k
       time1previousmarket = 0;  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
       RealisedOutputMarket_Real_perBattery = zeros(1,N_Bat);
        Real_RealisedOutputpreviousmarket = 0;
       RealisedOutputmarket = zeros(1,N_Bat);
       Real_RealisedOutputmarket = 0;
       Sum_D_minus_P_Comm=0;
        %######### Very important: because for each Market Time interval
        %the strategy is different, we need to go back to the initial state
        %at each new time interval.
            RealisedOutputMarketVectorprevious= 0;
        indexVector(indexbid+1) = k;
        soctargettypeMarket = 0;
        indexbid = indexbid+1;
        GridPower= zeros(1,N_Bat);
        AdditionalEnergy= zeros(1,N_Bat);
        Percentavail= 1/N_Bat*ones(1,N_Bat);
        Percentacceptance= 1/N_Bat*ones(1,N_Bat);

    end    

         if soctargettype ==1  %the SoC target must be met because of price change or other requirement
            CurrentTime = k-1-StartTime; %the time from which we recompute the energy that we still need
             for d = 1:N_Bat
                time1 = 2*(Pbatmax*Efficiency*(timeinterval-CurrentTime/2+1/2)-Et(d)*(CurrentTime/2+1/2)+abs(RealisedOutput(d))-min(abs(NeededEnergy(d)),Pbatmax*Efficiency*timeinterval))/(Pbatmax*Efficiency-Et(d)+0.000000001); %determines the the time at which we need to use the battery at its maximum power

                     if time1 > timeinterval %if we do not need to be at full power within this time interval
                        slope = 2*(sign(NeededEnergy(d))*min(abs(NeededEnergy(d)),Pbatmax*Efficiency*timeinterval)-RealisedOutput(d) - Et(d)*(timeinterval-CurrentTime))/((timeinterval-CurrentTime)*(timeinterval-CurrentTime+1));
                    else
                        slope = (sign(NeededEnergy(d))*Pbatmax*Efficiency-Et(d))^2/(sign(NeededEnergy(d))*Pbatmax*Efficiency*(2*timeinterval+1-2*CurrentTime)-Et(d)+2*(RealisedOutput(d)-sign(NeededEnergy(d))*min(abs(NeededEnergy(d)),Pbatmax*Efficiency*timeinterval))+0.000001);
                    end 
                    if abs(sign(NeededEnergy(d))*Pbatmax*Efficiency-Et(d) ) < epsilon  %if we already reached the maximum power at previous step while we were still ramping (happens when time1 is an integer
                        time1 = time1previous(d); %we keep the same value as previous iteration, as it will still be the maximum power for the battery
                        slope=slopeprevious(d);
                    end
                    Ebatreal(k,d) = slope+Et(d);
                    if k-StartTime >time1 % if we have already reached the time at which the battery should provide all its energy, then we must set its power to Pmax
                       Ebatreal(k,d)=sign(NeededEnergy(d))*Pbatmax;
                       EbatSoCNeed = sign(NeededEnergy(d))*Pbatmax;

                    end   
                Ebatreal(k,d) = min(Pbatmax, max(-Pbatmax,Ebatreal(k,d)));
                    RealisedOutput(d) = RealisedOutputprevious(d) + Ebatreal(k,d);
                    if sign( NeededEnergy(d) - RealisedOutput(d))*sign(NeededEnergy(d)) <0 && soctargettype==1 %if there is more energy provided than required %if we already exported too much energy compared to our real, we need to reduce the amount of eergy exported
                       Ebatreal(k,d) = Ebatreal(k,d) + NeededEnergy(d) - RealisedOutput(d) ;%max(-Pbatmax, min(Pbatmax,NeededEnergy-RealisedOutput)); %so we charge the battery with the surplus of energy. making sure the battery's max power is below Pmax
                       boolEt = 1;
                    end 
                        Ebatreal(k,d) = min(Pbatmax, max(-Pbatmax,Ebatreal(k,d)));

                    SoCreal(k,d) = SoCprevious(d)+Ebatreal(k,d);
                    [Ebatreal(k,d),SoCreal(k,d)] = coercsoc(SoCreal(k,d),SoCmax,SoCmin,Ebatreal(k,d),1,LastSoCtarget(d), NextSoCtarget,soctargettype,Pbatmax);  %1 because we consider energy and not power.

                    SoCreal(k,d) = SoCprevious(d)+Ebatreal(k,d);
                    EbatSoCNeed = Ebatreal(k,d);
                    time1previous(d) = time1;
                    slopeprevious(d) = slope;    
                     Et(d) = EbatSoCNeed ; %Ebatreal(k);
                     if boolEt ==1
                        Et(d) = 0;
                        boolEt=0;
                     end
                       s = -sign(EbatSoCNeed);
                     a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
                ToGridreal(d) =  max(0,max(0,Porigin(k,d)-Dorigin(k,d)- EbatSoCNeed*a)); %what we export at this time step
                FromGridreal(d) = max(0,max(0,Dorigin(k,d)-Porigin(k,d)+ EbatSoCNeed*a));
             end
         else
             for d = 1:N_Bat
                s = sign(Porigin(k,d)-Dorigin(k,d));
                a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency); 
                Ebatreal(k,d) = min(Pbatmax, max(-Pbatmax,(Porigin(k,d)-Dorigin(k,d))*a));
            SoCreal(k,d) = SoCprevious(d)+Ebatreal(k,d);
            [Ebatreal(k,d),SoCreal(k,d)] = coercsoc(SoCreal(k,d),SoCmax,SoCmin,Ebatreal(k,d),1,LastSoCtarget(d), NextSoCtarget,soctargettypeMarket,Pbatmax);  %1 because we consider energy and not power.
            SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d);

             end


         end

    CurrentTimemarket = k-1-StartTimemarket; %the time from which we recompute the energy that we still need
    EfficiencyA = Efficiency;
    if NeededOutputref >epsilon*N_Bat/10
        EbatSoc = Ebatreal(k,:);
         N_Bat_C = N_Bat;
        if PowertargettypeforNextTimeSetp == 2 %if we need to distribut extra energy to batteries that do not have enough, in preparation of next time interval where there is a big requirement
            N_Bat_C = 0; %Number of batteries that will produce energy to fill other batteries that are empty
            for d=1:N_Bat
               if  SoCreal(k-1,d)>NextSoCtarget
                    N_Bat_C = N_Bat_C+1;
               end
            end
        end


         for d=1:N_Bat
             if PowertargettypeforNextTimeSetp < 2 || SoCreal(k-1,d)>NextSoCtarget

                    time1market = 2*(sign(NeededOutput(d))*Pbatmax*N_Bat_C*EfficiencyA*(timeintervalmarket-CurrentTimemarket/2+1/2)-Etmarket(d)*(CurrentTimemarket/2+1/2)+RealisedOutputmarket(d)-sign(NeededOutput(d))*min(abs(NeededOutput(d)),Pbatmax*N_Bat_C*EfficiencyA*timeintervalmarket))/(sign(NeededOutput(d))*Pbatmax*N_Bat_C*EfficiencyA-Etmarket(d)+0.000000001); %determines the the time at which we need to use the battery at its maximum power
                    if time1market > timeintervalmarket %if we do not need to be at full power within this time interval
                            slopemarket = 2*(sign(NeededOutput(d))*min(abs(NeededOutput(d)),Pbatmax*N_Bat_C*EfficiencyA*timeintervalmarket)-RealisedOutputmarket(d) - Etmarket(d)*(timeintervalmarket-CurrentTimemarket))/((timeintervalmarket-CurrentTimemarket)*(timeintervalmarket-CurrentTimemarket+1));
                    else
                        slopemarket = (sign(NeededOutput(d))*Pbatmax*N_Bat_C*EfficiencyA-Etmarket(d))^2/(sign(NeededOutput(d))*Pbatmax*N_Bat_C*EfficiencyA*(2*timeintervalmarket+1-2*CurrentTimemarket)-Etmarket(d)+2*(RealisedOutputmarket(d)-sign(NeededOutput(d))*min(abs(NeededOutput(d)),Pbatmax*N_Bat_C*EfficiencyA*timeintervalmarket))+0.000001);
                    end 
                    if abs(sign(NeededOutput(d))*Pbatmax*N_Bat_C*EfficiencyA-Etmarket(d) ) < 0  %if we already reached the maximum power at previous step while we were still ramping (happens when time1 is an integer
                        time1market = time1previousmarket; %we keep the same value as previous iteration, as it will still be the maximum power for the battery
                        slopemarket=slopepreviousmarket;
                    end
                    if k-StartTimemarket >time1market-1 && time1market >= 0 % if we have already reached the time at which the battery should provide all its energy, then we must set its power to Pmax
                       slopemarket = sign(NeededOutput(d))*Pbatmax*N_Bat_C;
                    end
                    ErequiredMarket = Etmarket(d) + slopemarket; % this is the quantity of energy we need to export to meet the bid commitment
                    Etmarket(d) =  ErequiredMarket;
                    EbatrequiredMarket(d) = min(0,ToGridreal(d) - ErequiredMarket/N_Bat_C - FromGridreal(d)); % This is what we actually need from the battery as the rest is provided by P-D
                    s = sign(EbatrequiredMarket(d));  %############### EfficiencyWork
                    a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency); 
                    if NeededOutput(d) >0  && EbatrequiredMarket(d)>0 %or if NeededOutput(d) >0 ??  or sign(EbatrequiredMarket)<0 ??
                        Ebatreal(k,d) = min(Pbatmax, max(-Pbatmax, Ebatreal(k,d) + EbatrequiredMarket(d)*a)); %What we need from the battery is what is needed for following the SoC + what is needed to commit to the bid
                    elseif (NeededOutput(d) - RealisedOutputmarket(d))> N_Bat_C*Pbatmax*(timeintervalmarket-CurrentTimemarket)
                        Ebatreal(k,d) = -Pbatmax;
                    else
                        Ebatreal(k,d) = min(Pbatmax, max(-Pbatmax, Ebatreal(k,d) + EbatrequiredMarket(d)*a)); %What we need from the battery is what is needed for following the SoC + what is needed to commit to the bid
                    end
                    s = -sign(Ebatreal(k,d));
                    a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
                    SoCreal(k,d) = SoCprevious(d)+Ebatreal(k,d);
                    [Ebatreal(k,d),SoCreal(k,d)] = coercsoc(SoCreal(k,d),SoCmax,SoCmin,Ebatreal(k,d),1,LastSoCtarget(d), NextSoCtarget,soctargettypeMarket,Pbatmax);  %1 because we consider energy and not power.
                    SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d);
                    time1previousmarket = time1market;
                    slopepreviousmarket = slopemarket; 
             else
                    Ebatreal(k,d) = 0;
             end


         end
    else
        for d = 1:N_Bat
            if sign(Ebatreal(k,d))>0
            s = sign(Ebatreal(k,d));
            s = sign(Porigin(k,d)-Dorigin(k,d));
            a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);                 
                Ebatreal(k,d) = min(Pbatmax,max(-Pbatmax,max(0,max(Ebatreal(k,d),(Porigin(k,d)-Dorigin(k,d))*a))));

            elseif (sign(Ebatreal(k,d))<=0 && CurrentModeVectorreal(k)==12)   % we want to make sure the battery does not stop following the P-D if it reached its SoCtarget
                s = sign(Porigin(k,d)-Dorigin(k,d));
                a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);

                 Ebatreal(k,d) = min(Pbatmax, max(-Pbatmax,min(Ebatreal(k,d),(Porigin(k,d)-Dorigin(k,d))*a)));
            end
            SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d);
            [Ebatreal(k,d),SoCreal(k,d)] = coercsoc(SoCreal(k,d),SoCmax,SoCmin,Ebatreal(k,d),1,LastSoCtarget(d), NextSoCtarget,soctargettype,Pbatmax);  %1 because we consider energy and not power.
            SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d);
        end

    end
    for d = 1:N_Bat
             if(Real_RealisedOutputmarket)>max(NeededOutput(d),NeededOutputref)+epsilon %%&& NeededOutputref <epsilon %&& EnergieAvailable>epsilon
                      if EnergieAvailable>epsilon
                        s = sign(Porigin(k,d)-Dorigin(k,d));
                                            a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
                         Ebatreal(k,d)=min(Pbatmax, max(-Pbatmax,(Porigin(k,d)-Dorigin(k,d))*a+max(0,EnergieAvailable)*Efficiency*percentage_Allocate_EnergieAvailable(d)));
                         SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d);
                         [Ebatreal(k,d),SoCreal(k,d)] = coercsoc(SoCreal(k,d),SoCmax,SoCmin,Ebatreal(k,d),1,LastSoCtarget(d), NextSoCtarget,0,Pbatmax);  %1 because we consider energy and not power.
                         SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d); 
                      else%if PowertargettypeforNextTimeSetp==0
                        s = sign(Porigin(k,d)-Dorigin(k,d));
                        a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
                         Ebatreal(k,d)=min(Pbatmax, max(-Pbatmax,(Porigin(k,d)-Dorigin(k,d))*a));
                         SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d);
                         [Ebatreal(k,d),SoCreal(k,d)] = coercsoc(SoCreal(k,d),SoCmax,SoCmin,Ebatreal(k,d),1,LastSoCtarget(d), NextSoCtarget,0,Pbatmax);  %1 because we consider energy and not power.
                         SoCreal(k,d) =  SoCprevious(d)+Ebatreal(k,d);          
                      end
             end
        VarDebug(k) = max(0,EnergieAvailable);

        s = -sign(Ebatreal(k,d));
        a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
        RealisedOutput(d) = SoCreal(k,d)- LastSoCtarget(d); % RealisedOutputprevious + Ebatreal(k);
        RealisedOutputprevious(d) = RealisedOutput(d);
         SoCprevious(d) = SoCreal(k,d);
         RealisedOutputMarket_Real_perBattery(d) = RealisedOutputMarket_Real_perBattery(d) + Porigin(k,d)-Dorigin(k,d)- Ebatreal(k,d)*a;
        s = -sign(Ebatreal(k,d));
        a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
        RealGridPower(d) = Porigin(k,d)-Dorigin(k,d)- Ebatreal(k,d)*a;
    end
    
    for d = 1:N_Bat
       Household_Power(k,d) =  RealGridPower(d) ;
    end
    ToGridTotal(k) = max(0,sum(RealGridPower)) ; %sum(ToGridreal);
    FromGridTotal(k) = min(0, sum(RealGridPower)); %   sum(FromGridreal);
    if k>1   
      RealisedOutputMarketVector(k) = RealisedOutputMarketVectorprevious+ ToGridTotal(k) + FromGridTotal(k);% sum(RealisedOutputMarket_Real_perBattery); %mean(RealisedOutputmarket); %max(0,Porigin(k)-Dorigin(k)- Ebatreal(k)*a)-max(0,Dorigin(k)-Porigin(k)+ Ebatreal(k)*a);
      RealisedOutputMarketVector2(k) = Real_RealisedOutputmarket;% sum(RealisedOutputMarket_Real_perBattery); %mean(RealisedOutputmarket); %max(0,Porigin(k)-Dorigin(k)- Ebatreal(k)*a)-max(0,Dorigin(k)-Porigin(k)+ Ebatreal(k)*a);

    else 
      RealisedOutputMarketVector(k) = ToGridTotal(k) + FromGridTotal(k);% sum(RealisedOutputMarket_Real_perBattery); %mean(RealisedOutputmarket); %max(0,Porigin(k)-Dorigin(k)- Ebatreal(k)*a)-max(0,Dorigin(k)-Porigin(k)+ Ebatreal(k)*a);
      RealisedOutputMarketVector2(k) = Real_RealisedOutputmarket;% sum(RealisedOutputMarket_Real_perBattery); %mean(RealisedOutputmarket); %max(0,Porigin(k)-Dorigin(k)- Ebatreal(k)*a)-max(0,Dorigin(k)-Porigin(k)+ Ebatreal(k)*a);

    end
    RealisedOutputMarketVectorprevious = RealisedOutputMarketVector(k);
      % GridPowerVector(k) = GridPower(d);
Sum_D_minus_P_Comm = Sum_D_minus_P_Comm+ max(0,SumDorigin(k)-SumPorigin(k));
CompA(k) = Real_RealisedOutputmarket-NeededOutput(d);
CompB(k) = EnergieAvailable;
end
From_Grid_real = - FromGridTotal;
To_Grid_real = ToGridTotal;




%% Computation of benefits using heuristic Battery Control
%Individual batteries case:
for ind =1:size(Porigin,2)
    PH = Porigin(:,ind);
    DH = Dorigin(:,ind);
    MaxBcap = SoCmax;% Maximum battery capcity in WHr [Variable Parameter]
    MinBcap = 0;% Minimum battery capacity at 80% DoD
    IBC = SoCInit; % Initial Battery Capacity (IBC)of 60%MacBcap
    MaxPower = Pbatmax;
    SoCH = zeros(length(PH),1);
    ToGridHeuristic = zeros(length(PH),1);
    FromGrid_Heuristic = zeros(length(PH),1);
    EbatMax = zeros(length(PH),1);
    EBat = zeros(length(PH),1);
    SoCH(1) = IBC;
    T = 1;  %minutely data (1 min)
    Ec = Efficiency;
    Ed = Efficiency;
    for i = 1:1:length(PH)
        if PH(i)>= DH(i)
            if i == 1
                EbatMax(i) = 1/T*min((PH(i)-DH(i))*T, MaxPower);
                EBat(i) = min(EbatMax(i)*Ec, MaxBcap-SoCH(1));
                SoCH(i) = SoCH(1) + EBat(i);
                ToGridHeuristic(i,:) = (PH(i)-DH(i)-EBat(i)/Ec);
            else
                EbatMax(i) = 1/T*min((PH(i)-DH(i))*T, MaxPower);
                EBat(i) = min(EbatMax(i)*Ec, MaxBcap-SoCH(i-1));
                SoCH(i) = SoCH(i-1) + EBat(i);
                ToGridHeuristic(i,:) = (PH(i)-DH(i)-EBat(i)/Ec);
            end
        else
            if i==1
                EbatMax(i) = 1/T*min(abs((PH(i)-DH(i))*T), MaxPower);
                EBat(i) = min(EbatMax(i)/Ed, SoCH(1)-MinBcap);
                SoCH(i) = SoCH(1) - EBat(i);
                FromGrid_Heuristic(i,:) = (PH(i)-DH(i)+ EBat(i)*Ed);
            else
                EbatMax(i) = 1/T*min(abs((PH(i)-DH(i))*T), MaxPower);
                EBat(i) = min(EbatMax(i)/Ed, SoCH(i-1)-MinBcap);
                SoCH(i) = SoCH(i-1) - EBat(i);
                FromGrid_Heuristic(i,:) = (PH(i)-DH(i)+ EBat(i)*Ed);
            end
        end
    end
    
Export_Market_Heuristic = zeros(size(t3));
Import_Market_Heuristic = zeros(size(t3));
j = 1;
for k = 1:size(t3,1)-1
    Export_Market_Heuristic(k+1) = max(0,sum(ToGridHeuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_Heuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    Import_Market_Heuristic(k+1) = min(0,sum(ToGridHeuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_Heuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
end 
    
  Total_Bill_Heuristic(ind+1,1) = sum(-Import_Market_Heuristic.*BP3)-sum(Export_Market_Heuristic(1:size(t3,1)).*SP3);  
  Total_Bill_HeuristicNoSP(ind+1,1) = sum(-Import_Market_Heuristic.*BP3)-sum(Export_Market_Heuristic(1:size(t3,1)).*0);  
end

%If it was only one big prosumer:
    PH = SumPorigin;
    DH = SumDorigin;
    MaxBcap = SoCmax*N_Bat;% Maximum battery capcity in WHr [Variable Parameter]
    MinBcap = 0;% Minimum battery capacity at 80% DoD
    IBC = SoCInit*N_Bat; % Initial Battery Capacity (IBC)of 60%MacBcap
    MaxPower = Pbatmax*N_Bat;
    SoCH = zeros(length(PH),1);
    ToGridHeuristic = zeros(length(PH),1);
    FromGrid_Heuristic = zeros(length(PH),1);
    EbatMax = zeros(length(PH),1);
    EBat = zeros(length(PH),1);
    SoCH(1) = IBC;
    T = 1;  %minutely data (1 min)
    Ec = Efficiency;
    Ed = Efficiency;
    for i = 1:1:length(PH)
        if PH(i)>= DH(i)
            if i == 1
                EbatMax(i) = 1/T*min((PH(i)-DH(i))*T, MaxPower);
                EBat(i) = min(EbatMax(i)*Ec, MaxBcap-SoCH(1));
                SoCH(i) = SoCH(1) + EBat(i);
                ToGridHeuristic(i,:) = (PH(i)-DH(i)-EBat(i)/Ec);
            else
                EbatMax(i) = 1/T*min((PH(i)-DH(i))*T, MaxPower);
                EBat(i) = min(EbatMax(i)*Ec, MaxBcap-SoCH(i-1));
                SoCH(i) = SoCH(i-1) + EBat(i);
                ToGridHeuristic(i,:) = (PH(i)-DH(i)-EBat(i)/Ec);
            end
        else
            if i==1
                EbatMax(i) = 1/T*min(abs((PH(i)-DH(i))*T), MaxPower);
                EBat(i) = min(EbatMax(i)/Ed, SoCH(1)-MinBcap);
                SoCH(i) = SoCH(1) - EBat(i);
                FromGrid_Heuristic(i,:) = (PH(i)-DH(i)+ EBat(i)*Ed);
            else
                EbatMax(i) = 1/T*min(abs((PH(i)-DH(i))*T), MaxPower);
                EBat(i) = min(EbatMax(i)/Ed, SoCH(i-1)-MinBcap);
                SoCH(i) = SoCH(i-1) - EBat(i);
                FromGrid_Heuristic(i,:) = (PH(i)-DH(i)+ EBat(i)*Ed);
            end
        end
    end
    
Export_Market_Heuristic = zeros(size(t3));
Import_Market_Heuristic = zeros(size(t3));
j = 1;
for k = 1:size(t3,1)-1
    Export_Market_Heuristic(k+1) = max(0,sum(ToGridHeuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_Heuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    Import_Market_Heuristic(k+1) = min(0,sum(ToGridHeuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_Heuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
end 
    
  Total_Bill_HeuristicAggregation = sum(-Import_Market_Heuristic.*BP3)-sum(Export_Market_Heuristic(1:size(t3,1)).*SP3);
  
  
  
%% Computation of bills with PV, but no Battery 
%Individual batteries case:
for ind =1:size(Porigin,2)
    PH = Porigin(:,ind);
    DH = Dorigin(:,ind);
    MaxBcap = 0;% Maximum battery capcity in WHr [Variable Parameter]
    MinBcap = 0;% Minimum battery capacity at 80% DoD
    IBC = 0; % Initial Battery Capacity (IBC)of 60%MacBcap
    MaxPower = 0;
    SoCH = zeros(length(PH),1);
    ToGridHeuristicNoBatt = zeros(length(PH),1);
    FromGrid_HeuristicNoBat = zeros(length(PH),1);
    EbatMax = zeros(length(PH),1);
    EBat = zeros(length(PH),1);
    SoCH(1) = IBC;
    T = 1;  %minutely data (1 min)
    Ec = Efficiency;
    Ed = Efficiency;
    for i = 1:1:length(PH)
        if PH(i)>= DH(i)
            if i == 1
                EbatMax(i) = 1/T*min((PH(i)-DH(i))*T, MaxPower);
                EBat(i) = min(EbatMax(i)*Ec, MaxBcap-SoCH(1));
                SoCH(i) = SoCH(1) + EBat(i);
                ToGridHeuristicNoBatt(i,:) = (PH(i)-DH(i)-EBat(i)/Ec);
            else
                EbatMax(i) = 1/T*min((PH(i)-DH(i))*T, MaxPower);
                EBat(i) = min(EbatMax(i)*Ec, MaxBcap-SoCH(i-1));
                SoCH(i) = SoCH(i-1) + EBat(i);
                ToGridHeuristicNoBatt(i,:) = (PH(i)-DH(i)-EBat(i)/Ec);
            end
        else
            if i==1
                EbatMax(i) = 1/T*min(abs((PH(i)-DH(i))*T), MaxPower);
                EBat(i) = min(EbatMax(i)/Ed, SoCH(1)-MinBcap);
                SoCH(i) = SoCH(1) - EBat(i);
                FromGrid_HeuristicNoBat(i,:) = (PH(i)-DH(i)+ EBat(i)*Ed);
            else
                EbatMax(i) = 1/T*min(abs((PH(i)-DH(i))*T), MaxPower);
                EBat(i) = min(EbatMax(i)/Ed, SoCH(i-1)-MinBcap);
                SoCH(i) = SoCH(i-1) - EBat(i);
                FromGrid_HeuristicNoBat(i,:) = (PH(i)-DH(i)+ EBat(i)*Ed);
            end
        end
    end
    
Export_Market_HeuristicNoBatt = zeros(size(t3));
Import_Market_HeuristicNoBatt = zeros(size(t3));
j = 1;
for k = 1:size(t3,1)-1
    Export_Market_HeuristicNoBatt(k+1) = max(0,sum(ToGridHeuristicNoBatt((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_HeuristicNoBat((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    Import_Market_HeuristicNoBatt(k+1) = min(0,sum(ToGridHeuristicNoBatt((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_HeuristicNoBat((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
end 
    
  Total_Bill_HeuristicNoBatt(ind+1,1) = sum(-Import_Market_HeuristicNoBatt.*BP3);%-sum(Export_Market_HeuristicNoBatt(1:size(t3,1)).*SP3);  
end



%% Computation of overall bill by considering  DUoS/TUoS/taxes, ...Metering, taxes ... 
%%Debug
Household_Power_Marketsize = zeros(size(t3,1),N_Bat);
Individual_Bills = zeros(size(t3,1),N_Bat);
Importof_fleet = zeros(size(t3,1),1);
PercentageElecCharge = 45/100;
for k = 1:size(t3,1)-1
    for d = 1:N_Bat
        Household_Power_Marketsize(k,d) = sum(Household_Power((k-1)*MarketTimeStep+1:(k)*MarketTimeStep,d)) ;
    end
    Extra_Energy = 0;
    for d = 1:N_Bat
        if Household_Power_Marketsize(k,d)>0
            Extra_Energy = Extra_Energy+Household_Power_Marketsize(k,d);
        end
    end    
    for d = 1:N_Bat
        if Household_Power_Marketsize(k,d)<0
            ConsumedEnergyfromFleet = Extra_Energy - max(0,(Extra_Energy+Household_Power_Marketsize(k,d)));
            Individual_Bills(k,d) = (abs(Household_Power_Marketsize(k,d))-ConsumedEnergyfromFleet)*BP3(k)+ConsumedEnergyfromFleet*PercentageElecCharge*BP3(k);
            Extra_Energy = Extra_Energy - ConsumedEnergyfromFleet;
            Importof_fleet(k) = Importof_fleet(k)+abs((Household_Power_Marketsize(k,d)));        
        end
    end
end


%% Computation of overall bill WITHOUT BATTERY by considering  DUoS/TUoS/taxes, ...Metering, taxes ... 
%%Debug
Household_Power_Marketsize_NoBAT = zeros(size(t3,1),N_Bat);
Individual_Bills_NoBAT = zeros(size(t3,1),N_Bat);
Importof_fleet_NoBAT = zeros(size(t3,1),1);
PercentageElecCharge = 45/100;
for k = 1:size(t3,1)-1
    for d = 1:N_Bat
        Household_Power_Marketsize_NoBAT(k,d) = sum(Porigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep,d)) - sum(Dorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep,d)) ;
    end
    Extra_Energy = 0;
    for d = 1:N_Bat
        if Household_Power_Marketsize_NoBAT(k,d)>0
            Extra_Energy = Extra_Energy+Household_Power_Marketsize_NoBAT(k,d);
        end
    end    
    for d = 1:N_Bat
        if Household_Power_Marketsize_NoBAT(k,d)<0
            ConsumedEnergyfromFleet_NoBAT = Extra_Energy - max(0,(Extra_Energy+Household_Power_Marketsize_NoBAT(k,d)));
            Individual_Bills_NoBAT(k,d) = (abs(Household_Power_Marketsize_NoBAT(k,d))-ConsumedEnergyfromFleet_NoBAT)*BP3(k)+ConsumedEnergyfromFleet_NoBAT*PercentageElecCharge*BP3(k);
            Extra_Energy = Extra_Energy - ConsumedEnergyfromFleet_NoBAT;
            Importof_fleet_NoBAT(k) = Importof_fleet_NoBAT(k)+abs((Household_Power_Marketsize_NoBAT(k,d)));        
        end
    end
end




%% Computation of overall bill WITHOUT distributed assets by considering VPP, with contracts with households and DUoS/TUoS/taxes, ...Metering, taxes ... 
%%Debug
Household_Power_Marketsize_VPP = zeros(size(t3,1),N_Bat);
VPP_Production_Marketsize_VPP = zeros(size(t3,1),1);
Individual_Bills_VPP = zeros(size(t3,1),N_Bat);
Importof_fleet_VPP = zeros(size(t3,1),1);
PercentageElecCharge = 45/100;
Total_Importof_fleet_VPP = zeros(size(t3,1),1);
Export_VPP = zeros(size(t3));
for k = 1:size(t3,1)-1
    
    
    Export_VPP(k) = (-Pbid(k)-max(0,Powerbat(k))*EfficiencyD+max(0,-Powerbat(k)/EfficiencyC));
    
    for d = 1:N_Bat
        Household_Power_Marketsize_VPP(k,d) = - sum(Dorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep,d)) ;
        VPP_Production_Marketsize_VPP(k) = VPP_Production_Marketsize_VPP(k)+sum(Porigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep,d)) ;
        Total_Importof_fleet_VPP(k) = Total_Importof_fleet_VPP(k)+Household_Power_Marketsize_VPP(k,d);
    end
    Extra_Energy = -Export_VPP(k)  %max(0,VPP_Production_Marketsize_VPP(k)+ Total_Importof_fleet_VPP(k));
%     for d = 1:N_Bat
%         if VPP_Production_Marketsize_VPP(k) >sum(Household_Power_Marketsize_VPP(k,:))
%             Extra_Energy = Extra_Energy+Household_Power_Marketsize_VPP(k,d);
%         end
%     end    
    for d = 1:N_Bat
        if Household_Power_Marketsize_VPP(k,d)<0
            ConsumedEnergyfromFleet_VPP = Extra_Energy - max(0,(Extra_Energy+Household_Power_Marketsize_VPP(k,d)));
            Individual_Bills_VPP(k,d) = (abs(Household_Power_Marketsize_VPP(k,d))-ConsumedEnergyfromFleet_VPP)*BP3(k)+ConsumedEnergyfromFleet_VPP*PercentageElecCharge*BP3(k);
            Extra_Energy = Extra_Energy - ConsumedEnergyfromFleet_VPP;
            Importof_fleet_VPP(k) = Importof_fleet_VPP(k)+abs((Household_Power_Marketsize_VPP(k,d)));        
        end
    end
end



%% Synthesis

% Export_Market_Heuristic = zeros(size(t3));
% Import_Market_Heuristic = zeros(size(t3));
Export_Market_real = zeros(size(t3));
Import_Market_real = zeros(size(t3));
Export_Market_without_Battery = zeros(size(t3));
Import_Market_without_Battery = zeros(size(t3));
Export_Market_PV = zeros(size(t3));
Export_Market_Wind = zeros(size(t3));
Export_all_VPP = zeros(size(t3));
DoriginMarketsize = zeros(size(t3));
j = 1;
for k = 1:size(t3,1)-1
    Export_Market_real(k+1) = max(0,sum(ToGridTotal((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGridTotal((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    Import_Market_real(k+1) = min(0,sum(ToGridTotal((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGridTotal((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    Export_Market_PV(k+1) = max(0,sum(Porigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep,size(Porigin,2)-1))) ;
    Export_Market_Wind(k+1) = max(0,sum(Porigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep,size(Porigin,2))));
    Export_all_VPP(k+1) = max(0,sum(SumPorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    DoriginMarketsize(k+1) = max(0,sum(SumDorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    Export_Market_without_Battery(k+1) = max(0,sum(SumPorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) - sum(SumDorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
    Import_Market_without_Battery(k+1) = min(0,sum(SumPorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) - sum(SumDorigin((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
%     Export_Market_Heuristic(k+1) = max(0,sum(ToGridHeuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_Heuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
%     Import_Market_Heuristic(k+1) = min(0,sum(ToGridHeuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)) + sum(FromGrid_Heuristic((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)));
end 

houseNumber = 1;
% Import_Cost(houseNumber,1) = sum(-Import_Market_real.*BP3);

Grid_ImbalanceCharges = 0;


Export_Market_real_Counting = zeros(size(Export_Market_real));
for k = 1:size(t3,1)
    if  To_Grid_bid(k)<=5
%        Export_Market_real_Counting(k) = 0;
    else
        Grid_ImbalanceCharges = Grid_ImbalanceCharges - (Export_Market_real(k)-To_Grid_bid(k))*SystemBuyPriceMWh(k)*sign(NetImbalanceVolumeMWh(k));
        Export_Market_real_Counting(k) =  To_Grid_bid(k);% Export_Market_real(k);
%         Export_Market_real_Counting(k) =min(Export_Market_real_Counting(k) ,To_Grid_bid(k)); %we cannot get more money than the commitment
    end
end
% Export_Revenue(houseNumber,1) = sum(Export_Market_real_Counting(1:size(t3,1)).*SP3);
% Total_Bill_RT(houseNumber,1) = sum(-Import_Market_real.*BP3)-sum(Export_Market_real_Counting(1:size(t3,1)).*SP3)+Grid_ImbalanceCharges;
% Total_Bill_No_Battery(houseNumber,1) = sum(-Import_Market_without_Battery.*BP3);%-sum(Export_Market_without_Battery(1:size(t3,1)).*SP3);
% Total_Bill_Heuristic(houseNumber,1) = sum(Total_Bill_Heuristic);
% Expected_Bill_After_Optim(houseNumber,1) = sum(From_Grid_bid(1:size(t3,1)).*BP3)-sum(To_Grid_bid(1:size(t3,1)).*SP3); %Beware it accurate only if there is no error forecast. To get real Value run it without uncertainty and 1/2 hour forecast)
fprintf("\n\n\n\n\n");
fprintf("Bill with real time control: %0.1f \n", sum(-Import_Market_real.*BP3)-sum(Export_Market_real_Counting(1:size(t3,1)).*SP3)+Grid_ImbalanceCharges);
fprintf("Equivalent Bill with real time control including DUoS and wholesale beenfits: %0.1f \n", sum(sum(Individual_Bills))-sum(Export_Market_real_Counting(1:size(t3,1)).*SP3)+Grid_ImbalanceCharges);
fprintf("S5:Revenues for the aggregator: %0.1f \n", sum(DoriginMarketsize.*BP3) - (sum(sum(Individual_Bills))-sum(Export_Market_real_Counting(1:size(t3,1)).*SP3)+Grid_ImbalanceCharges));
fprintf("S5: Bill with real time control including DUoS: %0.1f \n", sum(sum(Individual_Bills)));
fprintf("S5:Revenues for the aggregator from Market: %0.1f \n", (sum(Export_Market_real_Counting(1:size(t3,1)).*SP3)-Grid_ImbalanceCharges));

% fprintf("Bill without Battery (but with selling price): %0.1f \n",  sum(-Import_Market_without_Battery.*BP3)-sum(Export_Market_without_Battery(1:size(t3,1)).*SP3));
% Total_Bill_No_Battery(houseNumber,1) = sum(-Import_Market_without_Battery.*BP3)-sum(Export_Market_without_Battery(1:size(t3,1)).*0);
fprintf("Bill with aggregation of the fleet but without Battery (without selling price): %0.1f \n",  sum(sum(Individual_Bills_NoBAT)));

fprintf("Bill with heuristic control and aggregation  (only one big consumer): %0.1f \n", Total_Bill_HeuristicAggregation);
fprintf("Bill with heuristic control with SP: %0.1f \n", sum(Total_Bill_Heuristic));
fprintf("S2: Bill with heuristic control without SP: %0.1f \n", sum(Total_Bill_HeuristicNoSP));
fprintf("S1: Bill with PV, without battery and no SP: %0.1f \n", sum(Total_Bill_HeuristicNoBatt));
fprintf("S0: Bill without PV, without battery and no SP: %0.1f \n", sum(DoriginMarketsize.*BP3));

fprintf("Expected max Bill with optim control: %0.1f \n", sum(From_Grid_bid(1:size(t3,1)).*BP3)-sum(To_Grid_bid(1:size(t3,1)).*SP3));
Revenue_from_VPP_Before_Aggregation = sum(Export_Market_PV.*SP3)+sum(Export_Market_Wind.*SP3);  %under the assumption of perfet forecast
Revenue_AllAssets_VPP = sum(Export_all_VPP.*SP3);
fprintf("S3:Max expected Revenues for VPP with assets in one place without battery : %0.1f \n", Revenue_AllAssets_VPP);
VPP_Revenues = sum(ToGrid.*SP3);
Bills_of_Houseolds = sum(sum(Individual_Bills_VPP));
RevenueVPP = sum(DoriginMarketsize.*BP3) - sum(sum(Individual_Bills_VPP)) + VPP_Revenues;
fprintf("S6:Max Revenues for VPP with assets in one place with battery with households: %0.1f \n", RevenueVPP);
fprintf("S6:Households Bills: %0.1f \n", Bills_of_Houseolds);
fprintf("S6: Market Revenues: %0.1f \n", VPP_Revenues);

fprintf("S0: original Bill: %0.1f £\n", sum(DoriginMarketsize.*BP3));
fprintf("S1: original Bill: %0.1f £\n", sum(DoriginMarketsize.*BP3));
Error =zeros(size(t3));
countE = 0;
for k = 1:size(t3,1)
    if To_Grid_bid(k)>10000 % we only consider bids above 10 kWh
        Error(k) = abs(To_Grid_bid(k)-Export_Market_real(k))/ abs(To_Grid_bid(k));
        countE = countE+1;
    end
end
fprintf("Average Export error: %0.3f percent \n\n", sum(Error)*100/countE);
fprintf("Max Export error: %0.3f percent \n\n", max(Error)*100);

Self_consump_rateBefore = zeros(N_Bat,1);
Self_consump_rateAfter = zeros(N_Bat,1);

for d=1:N_Bat-2
    for k=1:size(t,1)
        if Porigin(k,d)>=Dorigin(k,d)
            Self_consump_rateBefore(d) = Self_consump_rateBefore(d)+1;
        end
        if Household_Power(k,d)<=0
            Self_consump_rateAfter(d) = Self_consump_rateAfter(d)+1;
        end        
    end
end
Self_consump_rateBefore = Self_consump_rateBefore/size(t,1);
Self_consump_rateAfter = Self_consump_rateAfter/size(t,1);
fprintf(" Old self consumption rate: %0.1f \n", mean(Self_consump_rateBefore)*100);
fprintf(" New self consumption rate: %0.1f \n", mean(Self_consump_rateAfter)*100);


save('Results_Scenario3_VPP_Agile_and_Nordpoolfull_No_forecast_error.mat','Total_Bill_RT','Total_Bill_No_Battery','Total_Bill_Heuristic','Expected_Bill_After_Optim','Export_Revenue','Import_Market_real','Export_Market_real_Counting','Import_Market_without_Battery','Export_Market_without_Battery','Import_Market_Heuristic','Export_Market_Heuristic','From_Grid_bid','To_Grid_bid'); 

fprintf("\n\n\n Monthly import Before PV and Batteries: %0.1f MWh \n", sum(SumDorigin)/1000000);
fprintf("Monthly import after PV and Batteries: %0.1f MWh \n", sum(From_Grid_real)/1000000);
fprintf("Monthly import before Batteries: %0.1f MWh \n", sum(max(0,SumDorigin-SumPorigin))/1000000);
sum(From_Grid_real)/1000 - sum(max(0,SumDorigin-SumPorigin))/1000;

close all

%% Display
Households_list = 1:houseNumber;

        figure4 = figure('Name','Bills','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(Households_list,Total_Bill_RT,Households_list,Total_Bill_No_Battery,Households_list,Expected_Bill_After_Optim, Households_list, Total_Bill_Heuristic ,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','Bill RT Optim','LineWidth',3,...
            'Color',[0.929411764705882 0.694117647058824 0.125490196078431]);
        set(plot1(2),'DisplayName','Bill wihout Battery','LineWidth',3,...
            'Color',[0.85 0.325 0.098]);
         set(plot1(3),'DisplayName','Optimal Expected Bill','Linewidth',1,...
            'Color',[0 0.447058823529412 0.741176470588235]);
           set(plot1(4),'DisplayName','Bill with Heuristic Control','linewidth',2,...
            'Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);



        % Create ylabel
        ylabel('Bill','FontWeight','bold');
        % Create xlabel
        xlabel('Number of Household','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        
% 
% fprintf("Total Bill without Batteries: %.1f \n", sum(Pbid.*BP3)-sum(Dbid(1:size(t3,1)).*SP3));
% 
% fprintf("Total Bill with Optimal Battery (Beware it is based on inaccurate forecast. To get real Value run it without uncertainty and 1/2 hour forecast): %.1f \n", ...
%     sum(From_Grid_bid(1:size(t3,1)).*BP3)-sum(To_Grid_bid(1:size(t3,1)).*SP3));
% Export_Market_real_Counting = Export_Market_real;
% for k = 1:size(t3,1)
%    if  To_Grid_bid(k)<=0.1
%        Export_Market_real_Counting(k) = 0;
%    end
% end
% fprintf("Total Bill with Real Time Batteries: %.1f \n", sum(-Import_Market_real.*BP3)-sum(Export_Market_real_Counting(1:size(t3,1)).*SP3));






 %% figure

if N_Bat==1
    
        figure4 = figure('Name','Timeseries','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(t3,SoC(1:size(t3,1))/N_Bat,t3, SoCTargets,t3,BP*10000000+2000,t3,SP*10000000+2000,t,CurrentModeVectorreal(1:size(t,1))*100,t,SoCTargetreal,...
            t3,SoCbid(1:size(t3,1))/N_Bat,t3,Pforecast,t3,Dforecast,t,SoCreal(:,1),t3,Required_Export_bid, t3, Export_Market_real, t,RealisedOutputMarketVector,t3, Import_Market_real,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','SoC','LineWidth',3,...
            'Color',[0.929411764705882 0.694117647058824 0.125490196078431]);
        set(plot1(2),'DisplayName','SoC Target','LineWidth',3,...
            'Color',[0.85 0.325 0.098]);
         set(plot1(3),'DisplayName','Buying Price','Linewidth',1,...
            'Color',[0 0.447058823529412 0.741176470588235]);
        set(plot1(4),'DisplayName','Selling Price','linewidth',1,...
            'Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);
        set(plot1(5),'DisplayName','Vector mode');
        set(plot1(6),'DisplayName','SoCTarget real');
        set(plot1(7),'DisplayName','SoC Bid');
        set(plot1(8),'DisplayName','P ');
        set(plot1(9),'DisplayName','D');
        set(plot1(10),'DisplayName','SoC Real 2',  'color', '#0A9BD2', 'linewidth',3);
        set(plot1(11),'DisplayName','Required Export bid', 'linewidth',4.5);
        set(plot1(12),'DisplayName',' Real Export ', 'linewidth',3.5);
        set(plot1(13),'DisplayName',' Realised Output ', 'linewidth',0.5, 'color',[0,0,0]);
        set(plot1(14),'DisplayName',' Real Import ', 'linewidth',3);


        % Create ylabel
        ylabel('Power','FontWeight','bold');
        % Create xlabel
        xlabel('Time','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        plot(t,Ebatreal*30,'LineWidth',1.5,'DisplayName','EbatReal')


        plot(t,(SumPorigin-SumDorigin)*30+1000,'DisplayName','SumPorigin - SumDorigin')
        % plot(t3,CurrentModeVectorbid*100)
        plot(t,-From_Grid_real(1:size(t,1))*30,'LineWidth',2,'DisplayName',' - From Grid ')
        yline(1000)
        plot(t3,-FromGrid,'LineWidth',1.5,'DisplayName','FromGrid After optim')
        % plot(t,GridPowerVector*30,'DisplayName','GridPower')
elseif N_Bat ==2
%     close all
            figure4 = figure('Name','Timeseries','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(t3,SoC(1:size(t3,1))/N_Bat,t3, SoCTargets,t3,BP*10000000+2000,t3,SP*10000000+2000,t,CurrentModeVectorreal(1:size(t,1))*100,t,SoCTargetreal,...
            t3,SoCbid(1:size(t3,1))/N_Bat,t3,P,t3,D,t,SoCreal(:,2),t3,Required_Export_bid, t3, Export_Market_real, t,RealisedOutputMarketVector,t,SoCreal(:,2),t3, Import_Market_real,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','SoC','LineWidth',3,...
            'Color',[0.929411764705882 0.694117647058824 0.125490196078431]);
        set(plot1(2),'DisplayName','SoC Target','LineWidth',3,...
            'Color',[0.85 0.325 0.098]);
         set(plot1(3),'DisplayName','Buying Price','Linewidth',1,...
            'Color',[0 0.447058823529412 0.741176470588235]);
        set(plot1(4),'DisplayName','Selling Price','linewidth',1,...
            'Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);
        set(plot1(5),'DisplayName','Vector mode');
        set(plot1(6),'DisplayName','SoCTarget real');
        set(plot1(7),'DisplayName','SoC Bid');
        set(plot1(8),'DisplayName','P ');
        set(plot1(9),'DisplayName','D');
        set(plot1(10),'DisplayName','SoC Real 2',  'color', '#0A9BD2', 'linewidth',3);
        set(plot1(11),'DisplayName','Required Export bid', 'linewidth',4.5);
        set(plot1(12),'DisplayName',' Real Export ', 'linewidth',3.5);
        set(plot1(13),'DisplayName',' Realised Output ', 'linewidth',0.5, 'color',[0,0,0]);
        set(plot1(14),'DisplayName','SoC Real 2',  'color', [1,0,1], 'linewidth',1);
        set(plot1(15),'DisplayName',' Real Import ', 'linewidth',3);


        % Create ylabel
        ylabel('Power','FontWeight','bold');
        % Create xlabel
        xlabel('Time','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar


        plot(t,(SumPorigin-SumDorigin)*30+1000,'DisplayName','SumPorigin - SumDorigin')
        plot(t,Ebatreal(:,2)*30, 'DisplayName','Ebatreal 2')
        % plot(t3,CurrentModeVectorbid*100)
        plot(t,-From_Grid_real(1:size(t,1))*30,'LineWidth',2,'DisplayName',' - From Grid ')
        yline(1000)
        plot(t,SoCreal(:,1),'DisplayName','SoC real 1')
        plot(t,(Ebatreal(:,1))*30,'DisplayName','Ebatreal1')
%         plot(t,(Ebatreal(:,3))*30,'DisplayName','Ebatreal3')
        plot(t3,-FromGrid,'LineWidth',1.5,'DisplayName','FromGrid After optim','color', '#F7DA06')
        
else
    
    
%          close all
        figure4 = figure('Name','Timeseries','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(t3,SoC(1:size(t3,1))/N_Bat,t3, SoCTargets,t3,BP*10000000+2000,t3,SP*10000000+2000,t,CurrentModeVectorreal(1:size(t,1))*100,t,SoCTargetreal,...
            t3,SoCbid(1:size(t3,1))/N_Bat,t3,P,t3,D,t,SoCreal(:,2),t3,Required_Export_bid, t3, Export_Market_real, t,RealisedOutputMarketVector,t,SoCreal(:,3),t3, Import_Market_real,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','SoC','LineWidth',3,...
            'Color',[0.929411764705882 0.694117647058824 0.125490196078431]);
        set(plot1(2),'DisplayName','SoC Target','LineWidth',3,...
            'Color',[0.85 0.325 0.098]);
         set(plot1(3),'DisplayName','Buying Price','Linewidth',1,...
            'Color',[0 0.447058823529412 0.741176470588235]);
        set(plot1(4),'DisplayName','Selling Price','linewidth',1,...
            'Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);
        set(plot1(5),'DisplayName','Vector mode');
        set(plot1(6),'DisplayName','SoCTarget real');
        set(plot1(7),'DisplayName','SoC Bid');
        set(plot1(8),'DisplayName','P ');
        set(plot1(9),'DisplayName','D');
        set(plot1(10),'DisplayName','SoC Real 2',  'color', '#0A9BD2', 'linewidth',3);
        set(plot1(11),'DisplayName','Required Export bid', 'linewidth',4.5);
        set(plot1(12),'DisplayName',' Real Export ', 'linewidth',3.5);
        set(plot1(13),'DisplayName',' Realised Output ', 'linewidth',0.5, 'color',[0,0,0]);
        set(plot1(14),'DisplayName','SoC Real 3',  'color', [1,0,1], 'linewidth',1);
        set(plot1(15),'DisplayName',' Real Import ', 'linewidth',3);


        % Create ylabel
        ylabel('Power','FontWeight','bold');
        % Create xlabel
        xlabel('Time','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar


        plot(t,(SumPorigin-SumDorigin)*30+1000,'DisplayName','SumPorigin - SumDorigin')
        plot(t,Ebatreal(:,2)*30, 'DisplayName','Ebatreal 2')
        % plot(t3,CurrentModeVectorbid*100)
        plot(t,-From_Grid_real(1:size(t,1))*30,'LineWidth',2,'DisplayName',' - From Grid ')
        yline(1000)
        plot(t,SoCreal(:,1),'DisplayName','SoC real 1')
        plot(t,(Ebatreal(:,1))*30,'DisplayName','Ebatreal1')
        plot(t,(Ebatreal(:,3))*30,'DisplayName','Ebatreal3')
        plot(t3,-FromGrid,'LineWidth',1.5,'DisplayName','FromGrid After optim','color', '#F7DA06')
        % plot(t,GridPowerVector*30,'DisplayName','GridPower')
 
end

   
%% Impact on grid
% close all

index3s = 24*48+1;
index3e = size(t3,1);
indexs = 24*60*24+1;
indexe = size(t,1);
        figure5 = figure('Name','Export Import 1','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure5);
        hold(axes1,'on');
%         yyaxis left
        plot1 = plot(t3(index3s:index3e),(Export_Market_real(index3s:index3e)*2+Import_Market_real(index3s:index3e)*2)/1000,...
            t3(index3s:index3e),(Export_Market_without_Battery(index3s:index3e)*2+Import_Market_without_Battery(index3s:index3e)*2)/1000,'Parent',axes1,'LineWidth',3);
         set(plot1(1),'DisplayName','Scenario 5 with RT algorithm',...
            'Color','#636363','LineStyle','-','LineWidth',2);      
         set(plot1(2),'DisplayName','Scenario 1 without batteries',...
            'Color','#71d160','LineStyle','-','LineWidth',4);
set(axes1,'ycolor','#000000');
ylim([-200 450])
 yticks([0 400])
datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Power (kW)');
 
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off
% set(legend1,...
%     'Position',[0.554861120765821 0.753682773033207 0.350520823678623 0.17163411631246]);

        figure6 = figure('Name','To Grid - from grid','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure6);
        hold(axes1,'on');
%         yyaxis left
        plot1 = plot(t(indexs:indexe),(To_Grid_real((indexs:indexe))*60-From_Grid_real((indexs:indexe))*60)/1000,t(indexs:indexe),(SumPorigin(indexs:indexe)*60 - SumDorigin(indexs:indexe)*60)/1000,'Parent',axes1,'LineWidth',3);
         set(plot1(1),'DisplayName','Scenario 5 with RT algorithm',...
            'Color','#636363','LineStyle','-','LineWidth',2);      
         set(plot1(2),'DisplayName','Scenario 1 without batteries',...
            'Color','#71d160','LineStyle','-','LineWidth',4);
set(axes1,'ycolor','#000000');
ylim([-200 450])
 yticks([0 400])
 datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Power (kW)');
 
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off






%% idem

        Housenumber =  44;


      figure6 = figure('Name','Household 44','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure6);
        hold(axes1,'on');
%         yyaxis left
        plot1 = plot(t(indexs:indexe),Household_Power((indexs:indexe),Housenumber)*30/1000, t(indexs:indexe),Porigin((indexs:indexe),Housenumber)*30/1000-Dorigin((indexs:indexe),Housenumber)*30/1000,'Parent',axes1,'LineWidth',2);
         set(plot1(1),'DisplayName','Scenario 5 with RT algorithm',...
            'Color','#636363','LineStyle','-','LineWidth',2);      
         set(plot1(2),'DisplayName','Scenario 1 without batteries',...
            'Color','#71d160','LineStyle','-','LineWidth',4);
set(axes1,'ycolor','#000000');
% ylim([-200 450])
%  yticks([0 400])
 datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Power (kW)');
 
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off








       figure4 = figure('Name','All households','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(t,To_Grid_real(1:size(t,1))*30-From_Grid_real(1:size(t,1))*30, t,SumPorigin(1:size(t,1))*30 - SumDorigin(1:size(t,1))*30,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','Import after algorithm','LineWidth',2);
        set(plot1(2),'DisplayName','Import before algorithm','LineWidth',2);

        % Create ylabel
        ylabel('Power','FontWeight','bold');
        % Create xlabel
        xlabel('Time','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar


        Housenumber =  44;
       figure4 = figure('Name','Household 44 bis','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(t,Household_Power(1:size(t,1),Housenumber)*30, t,Porigin(1:size(t,1),Housenumber)*30-Dorigin(1:size(t,1),Housenumber)*30,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','Import after algorithm','LineWidth',2);
        set(plot1(2),'DisplayName','Import before algorithm','LineWidth',2);

        % Create ylabel
        ylabel('Power','FontWeight','bold');
        % Create xlabel
        xlabel('Time','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar



Average_after = zeros(N_Bat,1);
Average_before = zeros(N_Bat,1);
PAR_after  = zeros(N_Bat,1);
PAR_before  = zeros(N_Bat,1);
for d = 1:N_Bat
    Average_after(d) = mean(Household_Power(:,d))*60/1000;
Average_before(d) = mean(Porigin(:,d)-Dorigin(:,d))*60/1000;
PAR_after(d) = max(Household_Power(:,d))/(Average_after(d)+0.000001)*60/1000;
PAR_before(d) = max(Porigin(:,d)-Dorigin(:,d))*60/1000/Average_before(d);
    
end
PAR_after(6) = 0;
AggregatedLoad_After = sum(Household_Power,2);
Average_Aggregated_After = mean(AggregatedLoad_After)*60/1000
Average_Aggregated_Before = mean(SumPorigin-SumDorigin)*60/1000
PAR_Before_Aggr = max(SumPorigin-SumDorigin)*60/1000/Average_Aggregated_Before
PAR_After_Aggr= max(AggregatedLoad_After)*60/1000/Average_Aggregated_After
(mean(To_Grid_real(1:size(t,1))*1-From_Grid_real(1:size(t,1))*1))*60/1000
mean(SumPorigin(1:size(t,1))*1 - SumDorigin(1:size(t,1))*1)*60/1000

     figure4 = figure('Name','PAR Average','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(1:N_Bat,Average_after,1:N_Bat,Average_before,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','Average After','LineWidth',2);
        set(plot1(2),'DisplayName','Average before algorithm','LineWidth',2);

        % Create ylabel
        ylabel('Power','FontWeight','bold');
        % Create xlabel
        xlabel('Time','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        
       figure4 = figure('Name','Par individual','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
        % Create multiple lines using matrix input to plot
        plot1 = plot(1:N_Bat,PAR_after,1:N_Bat,PAR_before,'Parent',axes1,'LineWidth',2);
        set(plot1(1),'DisplayName','PAR After','LineWidth',2);
        set(plot1(2),'DisplayName','PAR before algorithm','LineWidth',2);

        % Create ylabel
        ylabel('Power','FontWeight','bold');
        % Create xlabel
        xlabel('Time','FontWeight','bold');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',12,'FontWeight','bold','GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
  
        
        
        
        
      figure6 = figure('Name','Abs(PAR)','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure6);
        hold(axes1,'on');
%         yyaxis left
        plot1 = plot(1:N_Bat,abs(PAR_after),1:N_Bat,abs(PAR_before),'Parent',axes1,'LineWidth',2);
         set(plot1(1),'DisplayName','Scenario 5 with RT algorithm',...
            'Color','#636363','LineStyle','-','LineWidth',4);      
         set(plot1(2),'DisplayName','Scenario 1 without batteries',...
            'Color','#71d160','LineStyle','-','LineWidth',4);
set(axes1,'ycolor','#000000');
% ylim([-200 450])
%  yticks([0 400])
%  datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Power (kW)');
 xlabel('Household Number');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off






%% Paper figures

index3s = 24*48+1;
index3e = size(t3,1);
indexs = 24*60*24+1;
indexe = size(t,1);
%          close all
        figure4 = figure('Name','Export/Power','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure4);
        hold(axes1,'on');
%         yyaxis left
        plot1 = plot(t3(index3s:index3e),Required_Export_bid(index3s:index3e)*2/1000,t3(index3s:index3e),P(index3s:index3e)*2/1000,t3(index3s:index3e),D(index3s:index3e)*2/1000,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','Contractual Export Agreement',...
            'Color','#D8D8D8','LineStyle','-','LineWidth',5);
        set(plot1(2),'DisplayName','Forecast Aggregated Production',...
            'Color','#E97F00','LineStyle','-');
        set(plot1(3),'DisplayName','Forecast Aggregated Demand',...
            'Color','#000000','LineStyle','-');
      
set(axes1,'ycolor','#000000');
yticks([0 200 400])
datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Energy (kWh)');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
 legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off



index3s = 24*48+1;
index3e = size(t3,1);
indexs = 24*60*24+1;
indexe = size(t,1);
for i = index3s:index3e
    if Required_Export_bid(i)<1
        Export_Market_real(i) = 0; % we do not count what was not in a commitment
    end
end

figure5 = figure('Name','Selling Price','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure5);
        hold(axes1,'on');
%         yyaxis left
        plot1 = plot(t3(index3s:index3e),SP3(index3s:index3e)*10000/1.1,'Parent',axes1,'LineWidth',3);
         set(plot1(1),'DisplayName','Market Export Price',...
            'Color','#929292','LineStyle','-','LineWidth',3);
      
set(axes1,'ycolor','#000000');
 yticks([0.3 0.5 0.7])
datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('€/MWh');
 
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
        legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off
% set(legend1,...
%     'Position',[0.554861120765821 0.753682773033207 0.350520823678623 0.17163411631246]);

         figure6 = figure('Name','Real Export','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure6);
        hold(axes1,'on');
%         yyaxis left
        plot1 = plot(t3(index3s:index3e),Required_Export_bid(index3s:index3e)*2/1000,t3(index3s:index3e),Export_Market_real(index3s:index3e)*2/1000,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','Contractual Export Agreement',...
            'Color','#CCCCCC','LineStyle','-','LineWidth',4);
        set(plot1(2),'DisplayName','Realised Export',...
            'Color','#151E45','LineStyle','-','LineWidth',2);

      
set(axes1,'ycolor','#000000');
yticks([0 200 400])
datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Energy (kWh)');
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off
       





%% Households graphs
numhouse = 11;
        figure6 = figure('Name','Household 11','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure6);
        hold(axes1,'on');
                
        yyaxis right
        plot1 = plot(t(indexs:indexe),SoCreal(indexs:indexe,numhouse)/max(SoCreal(:,numhouse))*100,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','SoC',...
            'Color','#9AEBFF','LineStyle','-','LineWidth',6);
         ylabel('SoC (%)');
           yticks([0 50 100])
  
        set(axes1,'ycolor','#000000');
        
         yyaxis left
        plot1 = plot(t(indexs:indexe),Porigin(indexs:indexe,numhouse)*60/1000,t(indexs:indexe),Dorigin(indexs:indexe,numhouse)*60/1000,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','Production PV',...
            'Color','#E97F00','LineStyle','-','LineWidth',2);
        set(plot1(2),'DisplayName','Demand',...
            'Color','#000000','LineStyle','-','LineWidth',2);

      
        set(axes1,'ycolor','#000000');
        yticks([0 5 10])
datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Power(kW)');

                 
                 
                 
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off

numhouse = 70;
        figure6 = figure('Name','Household 70','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure6);
        hold(axes1,'on');
                
        yyaxis right
        plot1 = plot(t(indexs:indexe),SoCreal(indexs:indexe,numhouse)/max(SoCreal(:,numhouse))*100,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','SoC',...
            'Color','#9AEBFF','LineStyle','-','LineWidth',6);
         ylabel('SoC (%)');
           yticks([0 50 100])
  
        set(axes1,'ycolor','#000000');
        
         yyaxis left
        plot1 = plot(t(indexs:indexe),Porigin(indexs:indexe,numhouse)*60/1000,t(indexs:indexe),Dorigin(indexs:indexe,numhouse)*60/1000,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','Production PV',...
            'Color','#E97F00','LineStyle','-','LineWidth',2);
        set(plot1(2),'DisplayName','Demand',...
            'Color','#000000','LineStyle','-','LineWidth',2);

      
        set(axes1,'ycolor','#000000');
        yticks([0 5 10])
datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Power(kW)');

                 
                 
                 
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off


numhouse = 68;
        figure6 = figure('Name','Household 68','Color',[1 1 1]);
        % Create axes
        axes1 = axes('Parent',figure6);
        hold(axes1,'on');
                
        yyaxis right
        plot1 = plot(t(indexs:indexe),SoCreal(indexs:indexe,numhouse)/max(SoCreal(:,numhouse))*100,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','SoC',...
            'Color','#9AEBFF','LineStyle','-','LineWidth',6);
         ylabel('SoC (%)');
           yticks([0 50 100])
  
        set(axes1,'ycolor','#000000');
        
         yyaxis left
        plot1 = plot(t(indexs:indexe),Porigin(indexs:indexe,numhouse)*60/1000,t(indexs:indexe),Dorigin(indexs:indexe,numhouse)*60/1000,'Parent',axes1,'LineWidth',5);
         set(plot1(1),'DisplayName','Production PV',...
            'Color','#E97F00','LineStyle','-','LineWidth',2);
        set(plot1(2),'DisplayName','Demand',...
            'Color','#000000','LineStyle','-','LineWidth',2);

      
        set(axes1,'ycolor','#000000');
        yticks([0 5 10])
datetick('x','mmm dd','keeplimits','keepticks')
        ylabel('Power(kW)');

                 
                 
                 
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',40,'GridAlpha',0.1,'GridLineStyle',':','LineWidth',1.5,'XGrid','on','YGrid','on');
        % Create legend
legend(axes1,'show');
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        H=gca;
H.LineWidth=1; %change to the desired value 
set(gca, 'FontName', 'Times New Roman')
grid off



%% end of program



function [pbat,soc] = coercsoc(soc,socmax,socmin,pbat,timestep,lastsoc, nextsoc,soctargetstype,pbatmax)
bool1 = 1;
bool2 = 1;
 
    if (lastsoc - nextsoc < 0.001) && soc>nextsoc && soctargetstype==1  %if we already charged too much
        bool1=0;
        pbat = min(pbatmax, max(-pbatmax,pbat - (soc-nextsoc)*timestep)); %pbat + (soc-nextsoc)*timestep;
%         soc = nextsoc;
    elseif    (lastsoc - nextsoc >0.001) && soc< nextsoc && soctargetstype==1  %if we already discharged too much
        bool2 = 0;
        pbat =  min(pbatmax, max(-pbatmax, pbat-(soc - nextsoc)*timestep)); %pbat +(soc - nextsoc)*timestep; 
%         pbat =  min(pbatmax, max(-pbatmax,pbat +(soc - nextsoc)*timestep)); %pbat +(soc - nextsoc)*timestep; 

%         soc = nextsoc;
    end

    if soc > socmax && bool1
        pbat = min(pbatmax, max(-pbatmax,pbat - (soc-socmax)*timestep)); %pbat + (soc-socmax)*timestep;
%         soc = socmax;
    elseif soc<socmin && bool2
        pbat = min(pbatmax, max(-pbatmax,pbat -(soc - socmin)*timestep)); %pbat +(soc - socmin)*timestep;
%         soc = socmin;
    end
    

end



% t=1:size(P,1);
% plot(t,P,t,D,t,SoC(1:size(P,1)))%,t,x(LengthOptim*9+1:LengthOptim*10+1))
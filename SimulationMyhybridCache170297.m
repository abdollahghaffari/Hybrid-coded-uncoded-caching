% By Abdollah Ghaffari Sheshjavani
% the code is for the paper in title "Content caching for shared medium networks under heterogeneous users’ behaviors" Computer Networks, 2021
clear all
Policy=1; % 1=duplicate select files is consider for not cached files(uncoded requests)   else=  duplicate select files is not considered for not cached files
ShowDetailResults=false;
Method=1; %1=our hybrid   2=purecoded  3=pureUncoded
Time=2000; % rounds of runs

% input parameters ///////////////////////////////////////////////////////////////////////////////////////
ziph_parameter=1;  % content popularity parameter (Zipf distribution)
K=10; % number of SBS
N=1000; % number of video files
M=100; % SBS cache size

% these bellow parameters used for evaluating different distribution of users inside SBSs (here is 10 SBSs)
Z1=[1 1 1 1 1 5 15 20 25 30 ];
Z2=[0 2 2 3 7 11 14 16 20 25 ];
Z3=[0 2 4 6 9 11 14 16 18 20 ];
Z4=[1 3 5 7 9 11 13 15 17 19 ];
Z5=[2 4 6 8 9 11 12 14 16 18 ];
Z6=[3 5 7 9 9 11 11 13 15 17 ];
Z7=[4 6 9 9 9 10 11 12 14 16 ];
Z8=[5 7 9 9 9 10 11 12 13 15 ];
Z9=[6 8 9 9 9 10 11 12 12 14 ];
Z10=[8 9 9 9 9 10 11 11 12 12 ];
Z11=[10 10 10 10 10 10 10 10 10 10 ];
%z112=[1 1 1 1 1 1 1 1 1 1 ];
%z22=[0 0 0 0 1 1 1 2 2 3 ];
Z=Z11;      % selecting distribution of users (the Number of users in the coverage of each SBS) for simulation and analysis
Zmax=max(Z);    % max number of users in a SBS
ss=sum(Z);      % total number of all users in the system
Zvar=0;         % Variance of user distribution in SBSs
RateResults=zeros(Time);
for i=1:K
   Zvar=Zvar+(Z(i)-(ss/K))^2;
end
Zvar=Zvar/(K-1);
ZstandardDev=sqrt(Zvar);
% The end of input parameters////////////////////////////////////////////////////////////////////////////

out=optimizationDifferentZFunction(Z,K,N,M,ziph_parameter,Method);% 1 = hybrid coded-uncoded, 2=pure coded, 3=pure uncoded
if(out(1)==M)
    out(2)=M;
end
N1=out(1); % it is N*
M1=out(2); % it is it must lower than M
M2=M-M1;
T=(K*M2)/(N1-M1);
for l=1 :K
    a(l)=l;
end
if M2>0
    ComKT=combination(K,T);
    CombKT1=combination(K,T+1);

    %Coded cache placment phase

    slices=nchoosek(a,T);
    sizeTchooseK=size(slices(:,1));

    for l=1:K
        index=1;
        for j=1 : sizeTchooseK
            for i=1 :T
                if (slices(j,i)==l)
                   for h=1:T
                    SBSCodedcache(index,h,l)=slices(j,h);
                   end
                   index=index+1;
                end
            end
        end
    end
if ShowDetailResults ==true
    size(slices)
    slices
    SBSCodedcache
end
end
% end of coded cache placement Phase
allratecoded=0;
allratecodedFormula=0;
allrateuncoded=0;
allrate=0;
allrateFormula=0;
%EKi=zeros(1,Z);  % Expectation Ki - only for evaluation
%AvgEKi=0;
% it is for ziph distribution
Prob=zeros(1,N);
for idx=1:N
    Prob(idx)=(1/idx)^ziph_parameter;
end
sumprobe=zeros(1,N);
Prob=Prob/sum(Prob);
for idp=1:N
    sumprobe(1,idp)=sum(Prob(1:idp));
end
% the end for ziph distribution
for t=1 : Time
    SBSRequests=zeros(K,Zmax);
    SBScodedQueu=zeros(K,Zmax);
    requestQueusize=zeros(1,K);
    uncodedRequest=zeros(1,K);
    cacheMiss=zeros(1,Zmax);
    rateuncoded=0;
    ratecoded=0;
    ratecoded2=0;
    for i=1 : K
           
        % it is for ziph distribution
            for idx2=1: Z(i)
               temp=rand;
               for idx=1: N
                   if temp<=sumprobe(1,idx);
                        SBSRequests(i,idx2)=idx;
                        break;
                    end
                end
            end
       % the end for ziph distribution  ==    SBSRequests(i,:)=zipfrnd(ziph_parameter,N,Z);
      
       for j=1 : Z(i)
          SBSduplicateRequest=false;
          MBSUncodedduplicateRequest=false;
          if j>1 
              for n=1 : j-1
                   if SBSRequests(i,j)==SBSRequests(i,n)
                       SBSduplicateRequest=true;
                   end
              end
          end
          if i>1
              for n=1 : i-1
                  for l=1: Z(i)
                       if SBSRequests(i,j)==SBSRequests(n,l)
                           MBSUncodedduplicateRequest=true;
                       end
                  end
              end
          end
           if SBSRequests(i,j)>M1 
              % cacheMiss(j)=cacheMiss(j)+1;
               if SBSRequests(i,j)<=N1 &&  SBSduplicateRequest==false
                 requestQueusize(i)= requestQueusize(i)+1;
                 SBScodedQueu(i,requestQueusize(i))=SBSRequests(i,j);
               else if SBSRequests(i,j)>N1 
                      if Policy==1
                          if SBSduplicateRequest==false && MBSUncodedduplicateRequest==false;
                             uncodedRequest(i)=uncodedRequest(i)+1; 
                          end
                      else
                          uncodedRequest(i)=uncodedRequest(i)+1;
                      end
                    end
               end
           end 
       end
    end
    for p=1: K
    rateuncoded=rateuncoded+uncodedRequest(p);  
    end
if M2>0
    allXOR=nchoosek(a,T+1);
    for q=1: Zmax
       ki=0;
       index=0;
       for p=1 : K
            if requestQueusize(p)>0
                ki=ki+1;
                requestQueusize(p)=requestQueusize(p)-1;
            end
       end
       if ki>0
        % calculating rate
            for i=1: size(allXOR(:,1))
               flag=false;
               for j=1: T+1
                   if SBScodedQueu(allXOR(i,j),q)>0
                       flag=true;
                   end
               end
               if flag==true
                  index=index+1;
                  sendSlice(index,:)=allXOR(i,:);
                  ratecoded=ratecoded+(1/sizeTchooseK(1));
               end
           end
            showSend=zeros(index,T+1);
           for i=1: index
             showSend(i,:)=sendSlice(i,:);
           end
           if ShowDetailResults ==true
              showSend
           end
           % The end of calculating rate 
       
       % calculating rate with use formula
        notsending=0;
        if K-ki>=T+1
            notsending=combination(K-ki,T+1);
        end
        %EKi(q)=EKi(q)+ki;
        ratecoded2=ratecoded2+min ((CombKT1-notsending)/ComKT,(ki-((ki*M2)/(N1-M1)))*((N1-M1)/ki));% CombKT1 is combination(K,T+1)
      % The end of calculating rate with use formula
       end
    end
end
if ShowDetailResults ==true
    SBSRequests
    SBScodedQueu
end
RateResults(t)=ratecoded+rateuncoded;
allratecoded=allratecoded+ratecoded;
allratecodedFormula=allratecodedFormula+ratecoded2;
allrateuncoded=allrateuncoded+rateuncoded;
allrate=allratecoded+allrateuncoded;
allrateFormula=allratecodedFormula+allrateuncoded;
end
% calculating Average of results
allratecoded=allratecoded/Time;
allratecodedFormula=allratecodedFormula/Time;
allrateuncoded=allrateuncoded/Time;
allrate=allrate/Time;
allrateFormula=allrateFormula/Time;
standardDeviationR=0;
for i=1:Time
   standardDeviationR=standardDeviationR+(RateResults(i)-allrate)^2; 
end
standardDeviationR=sqrt(standardDeviationR/(Time-1));
%for i=1:Z
%EKi(i)=EKi(i)/Time;
%AvgEKi=AvgEKi+EKi(i);
%end
%AvgEKi=AvgEKi/Z;
%show results
%Time
%Z
Method
ziph_parameter
Z
Zvar
ZstandardDev
M1
N1
AnalyticalRate=out(3)
SimulationRate=allrate
standardDeviationR
confidenceInt95L=SimulationRate-1.96*(standardDeviationR/sqrt(Time))
confidenceInt95H=SimulationRate+1.96*(standardDeviationR/sqrt(Time))
%EKi
%AvgEKi
%allrateuncoded
%allratecoded
%allratecodedFormula
%allrateFormula
x = '=======================================';
x

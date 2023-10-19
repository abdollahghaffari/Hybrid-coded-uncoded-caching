% Abdollah Ghaffari sheshjavani 7/22/2018 , 30/10/1398
function out=newOptimizationDifferentZFunction(Z,K,N,M,ziph_parameter,Method)

%Method=1; % 1=our hybrid   2=purecoded  3=pureUncoded

cachceckstart=0;
cachchecksize=0;

Zmax=max(Z);
sigmaZ=sum(Z);

if Method==1
    cachchecksize=M-1;
end
if Method==3
    cachceckstart=M;
    cachchecksize=M;
end
% this is variables
BestBandwidth=-1;
bestN1=-1; % it is the best N*
bestM1=-1;
bestbase=-1;
bestuncached=-1;
%
allziph=0;
for l=1 : N
    allziph=allziph+ ((1/l)^ziph_parameter);
end
myfactorial(1)=1;
for l=1 : K
    myfactorial(l+1)=l*myfactorial(l); % notice ==> factorial(k)=myfactorial(k+1)
end
% this section compute combination(Z,j) or factorial(n)/(factorial(k)*factorial(n-k))with dynomic programing
for i=1 : K   
  temp1(i)=factorial(Z(i));
 % temp1(i)
end
combZ=zeros(K,Zmax);
for i=1 : K
    temp2=1;
    temp3=temp1(i);
    combZ(i,1)=1;
    for j=1 : Z(i)
        temp2=temp2*j;
        temp3=temp3/(Z(i)-j+1);
        combZ(i,(j+1))=temp1(i)/(temp2*temp3);%combination(Z,j);
    end
end
% the end of computing combination(Z,j)
MBSBackhaulOverhead=0;
SumofCodedZiph1=0;
SumofCodedZiphNew1=0;
%uncached1=0;
uncached2=0;
 
for N1=M : N 
%else if Policy2==1
    if uncached2==0
        for n=N1+1 : N
            uncached2=uncached2+ (1-(1-(((1/n)^ziph_parameter)/allziph))^sigmaZ); %(Z*K));
            %order1=order1+1;
        end
    else
            uncached2=uncached2 - (1-(1-(((1/N1)^ziph_parameter)/allziph))^ sigmaZ); %(Z*K)); % with dynomic programing we reduce the order of this section from N^2 to 2N 
            %order1=order1+1;
    end
       uncached=uncached2;
    lastM1=0;
    NewlastM1=0;
    %SumofCodedZiph2=SumofCodedZiph1;
    SumofCodedZiphNew2=SumofCodedZiphNew1;
    SumofCodedZiph=0;     
    zarib1=0;
    for M1=cachceckstart : cachchecksize
        Flag=false;
        if M1==M && N1==M
          MBSBackhaulOverhead=uncached;
          Flag=true;
        end
        if mod(K*(M-M1),(N1-M1))==0 || Flag==true  
          if Flag==true
              T=0;
          else
           T=K*(M-M1)/(N1-M1);
           base=Zmax*(K*(N1-M)/((N1-M1)+K*(M-M1))); % K-T/T+1
           MBSBackhaulOverhead=base+uncached;
           end
           expectation=0;
           if T>0
               %zarib1=1/(combination(K,T));
               zarib1=1/ (myfactorial(K+1)/(myfactorial(T+1)*myfactorial(K-T+1)));
%\\\\\\\\\\\\\\\\\\\\\  
               SumofP_Codedfiles=0;
               for n=M1+1 : N1
                   SumofP_Codedfiles=SumofP_Codedfiles+ ((1/n)^ziph_parameter);
                   %order2=order2+1;
               end  
               SumofP_Codedfiles=SumofP_Codedfiles/allziph;
               EP_codedfiles=SumofP_Codedfiles/(N1-M1);
               qcn(1)= SumofP_Codedfiles;
               for j=2: Zmax+1
                   qcn(j)= SumofP_Codedfiles-((j-1)*EP_codedfiles); 
                   if (j>N1-M1) 
                        qcn(j)=0;  
                   end
               end
 %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\            
               expectation=0;
               PZ= zeros(Zmax+1, Zmax+1);
               PZ(1,1)= 1.0;
               for s=2: Zmax+1
                    for qipz=1: s
                        if qipz == 1 
                            PZ(s,qipz)=PZ(s-1,qipz)*(1 - qcn(1));
                        else 
                            PZ(s,qipz) = (PZ(s-1,qipz) * (1 -qcn(qipz))+ PZ(s-1,qipz-1) * qcn(qipz-1)); 
                        end
                    end
               end
     %  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
               for i=1: Zmax 
                   for ki=0: K
                       PQ = zeros(K, K+1);
                       PQ(1,1)=1.0;
                       for c=1: K  
                            Pzi(c)=0;
                            if i<=Z(c)
                                for j=i :Z(c)
                                    Pzi(c)=Pzi(c)+PZ(Z(c)+1,j+1);
                                end
                            end
                %  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                            for qi=0: c
                                if qi == 0 
                                    PQ(c+1,qi+1)=PQ(c,qi+1)*(1 - Pzi(c));
                                else 
                                    PQ(c+1,qi+1) = (PQ(c,qi+1) * (1 - Pzi(c))+ PQ(c,qi) * Pzi(c)); 
                                end
                            end
                        end
                        if K-ki >=T+1
                             expectation=expectation+PQ(K+1,ki+1)*(myfactorial(K-ki+1)/(myfactorial(T+1+1)*myfactorial(K-ki-T)));
                        end
                   end
               end
           end

           %N1
           %M1
          % expectation
          % zarib1
          % MBSBackhaulOverhead
           MBSBackhaulOverhead=MBSBackhaulOverhead-(zarib1*expectation);
           %MBSBackhaulOverhead
          if BestBandwidth==-1
              BestBandwidth=MBSBackhaulOverhead;
              bestN1=N1;
              bestM1=M1;
              %bestN1
              %bestM1
              %BestBandwidth
          else if MBSBackhaulOverhead<=BestBandwidth
               BestBandwidth=MBSBackhaulOverhead;
               bestN1=N1;
               bestM1=M1;
               bestuncached=uncached;
               bestbase=base;
               %bestN1
               %bestM1
               %BestBandwidth
              end
          end
        end
     if N1==100
       	break;
     end
    end
end
% Show Results
%Zmax
%ziph_parameter
%bestN1
%bestM1
%bestuncached
%bestbase
%BestBandwidth
out(1)=bestN1;
out(2)=bestM1;
out(3)=BestBandwidth;
%x = '=======================================';
%x
end
%order1
%order2
%order3
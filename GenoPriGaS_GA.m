%-------------------[Genomic Privacy Game Solver] (v0.02)-----------------
%This program finds the best solution for sharing genomic summary statistics (MAF of SNPs) 
%under an economically motivated recipient(adversary)'s inference attack based on 
%a Stackelberg game model named Genomic Privacy Game.
%Inference attack is to infer if a targeted DNA is in a genome pool with
%publised summary statistics (i.e. minor allele frequency of SNPs)
%The sharer's strategy set is (relaxed) not constraint to releasing only 
%top SNPs.
%Genetic Algorithm (GA) is used to find the best solution, especially for large scale datasets.
%Copyright 2016 Zhiyu Wan
%HIPLAB, Department of Biomedical Informatics, Vanderbilt University
%Oct 26, 2016
%-------------------------------------------------------------------------

warning off %#ok<*WNOFF>
close all % delete all the figures whos handles are not hidden
clear
format compact
rng(0);
tic %start a stopwatch timer

%--Sensitivity analysis on the expected cost of penalty to the adversary--%

% % Game Configuration
Homer_Attack=0; % 1=Homer's attack; 0=Jordan's attack (default)
relax2=0; % 1=Adversary only uses released independent SNPs; 0=Adversry uses all released SNPs (default)
Demo=0; % 1=Show the GA running process; 0=Hide the GA running process (default)
N_Class=2; % 1=1-class classification (LR test: in pool or in population); 
           % 2=2-classes classification (LR test: in pool or in reference) (default)
Attacker_Baseline=0;% 1=adversary always attacks all individuals; 
                    % 0=adversary is acting rationally (economically motivated) (default)

N_T=2500;%size of the target set (default:2500)
U_factor=1;%The weight of utility / the weight of privacy (default:1)
L=360;% Loss to the sharer per successful attack ($) (default:360)
G=L; % Gain to the recipient per successful attack ($) (default:360)
Cm=0;% Cost of marketing (access each target) (default:0)
Cp=240;% Cost of expected penalty from Data Use Agreement (default:240)
K=36;% # of measure points in the sensitivity analysis (default:36)
sample_rate=0.05;%the prior probability a target is in the pool (default:0.05)

mafcutoff=0.0001;% cutoff on minor allele frequency (MAF) in the filtering process (default:0.0001)
ldcutoff=0.05;% cutoff on linkage disequilibrium (r2) in the filtering process (default:0.05)
Race=0;%1=EastAsian, 2=SouthAsian, 3=African, 4=European, 5=American, 0=All; (default:0) (Not effective in this version)
power=1;% maximal allowable detection power (true positive rate) for the adversary (default:1)
theta_miss=0.1;% maximal allowable missing rate for each SNP in the filtering process (default:0.1)

% % Genetic Algorithm Configuration
n1=268;% size of subpopulation 1 (better to be larger than the # of independt SNPs after filtering) (default:268)
n2=n1;% size of subpopulation 2 (better to be larger than the # of independt SNPs after filtering) (default:268)
T=100;% # of iterations (default:100)

p_E=0.2;%portion of elite (default:0.2)
p_m=0.05;%upper bound of mutation probability (default:0.05)
p_p=0.8;%production proportion (default:0.8)
c_f=0.4;%Crossover Fraction (default:0.4)
Selection_function=1;%1=Roulette;2=Tournament (default:1)
Crossover_function=3;%0=half-half;1=one point;2=two points;3=scatters (default:3)
Tournament_size=10;%Tournament size (default:10)
Interval=10;%Immigration Interval (default:10)
Direction=2;%Immigration Direction. 1=one-way;2=two-way; (default:2)
Frac=0.2;%Immigration Fraction (default:0.2)

load('Sphinx/Pool+Ref')%load input datasets: Pool (eMERGE network SPHINX dataset)
                       %and Reference (1000Genomes Phase 3 dataset)
Pool=Pool(:,1:11574);%keep only autosomal SNPs (Not necessary if the input dataset does not have autosomal SNPs)
Ref=Ref(:,1:11574);%keep only autosomal SNPs (Not necessary if the input dataset does not have autosomal SNPs)
Pop=[Pool;Ref];%Combine two dataset
pMark=[true(size(Pool,1),1);false(size(Ref,1),1)];%mark the individuals in pool
rMark=~pMark;%mark the individuals in reference

n_pool=sum(pMark);%calculate the # of individuals in pool
n_pop=length(pMark);%calculate the total # of individuals

% % % Filtering
% % Handling missing values
MissingSNP=sum(Pool==-1);%A vector, in which each element is # of individuals in pool with missing values for a SNP
RetainedSNP=MissingSNP<=n_pool*theta_miss;%A binary vector, in which each element indicates a SNP will be retained or not
RetainedID=1:size(Pool,2);%Create a ID for each individual in pool
Pool=Pool(:,RetainedSNP);%Discard all the SNPs that should not be retained from the pool
Ref=Ref(:,RetainedSNP);%Discard all the SNPs that should not be retained from the reference
RetainedID=RetainedID(RetainedSNP);%Record the IDs of those retained SNPs

% % Compute MAF
MissPool=(Pool==-1);%A binary matrix, in which each element indicates a SNP for an individual in pool is missing or not
MissRef=(Ref==-1);%A binary matrix, in which each element indicates a SNP for an individual in reference is missing or not
Poolfreq=sum(Pool+MissPool)/2./sum(1-MissPool);%Vector of MAFs in the pool
Reffreq=sum(Ref+MissRef)/2./sum(1-MissRef);%Vector of MAFs in the reference
if N_Class==1 % LR test: in pool or in population (can not tell)
    Allfreq=(Poolfreq.*sum(1-MissPool)+Reffreq.*sum(1-MissRef))./(sum(1-MissPool)+sum(1-MissRef));%Vector of MAFs in population
    Reffreq=Allfreq; %Rename
end

% % Compute utility score (weight)
Score=(abs(Poolfreq-Reffreq));%Compute the utility weight for each SNP as the absolute difference 
                              %between the MAF in the pool and MAF in the reference

% % Remove SNPs with MAF < mafcutoff
Reject=Poolfreq<mafcutoff|Poolfreq>(1-mafcutoff);%Mark those SNPs with too small MAf that should be rejected

% Retain SNPs
RetainedPool=Pool(:,~Reject);% Discard SNPs with too small MAF from the pool
m=size(RetainedPool,2);% # of SNPs in the new pool
RetainedScore=Score(~Reject);% Discard SNPs with too small MAF in the utilty weight vector
RetainedRef=Ref(:,~Reject);% Discard SNPs with too small MAF freom the reference
RetainedPoolfreq=Poolfreq(:,~Reject);% Discard SNPs with too small MAF from the poof MAF vector
RetainedReffreq=Reffreq(:,~Reject);% Discard SNPs with too small MAF from the reference MAF vector
RetainedID=RetainedID(~Reject);%Record the IDs of remaining SNPs

% % Sort rank
[sortedscore,sortedindex]=sort(RetainedScore,'descend');%Sort SNPs according to their utility weight in a descending order
SortedPool=RetainedPool(:,sortedindex);%Adjust the order of SNPs in the pool accordingly
SortedRef=RetainedRef(:,sortedindex);%Adjust the order of SNPs in the reference accordingly
SortedPoolfreq=RetainedPoolfreq(:,sortedindex);%Adjust the order of SNPs in the pool-MAF vector accordingly
SortedReffreq=RetainedReffreq(:,sortedindex);%Adjust the order of SNPs in the reference-MAF vector accordingly
SortedID=RetainedID(sortedindex);%Adjust the IDs of remaining SNPs accordingly

% % Compute the correlation matrix
AbnormalRef=(SortedRef==-1);%A binary matrix, each element indicates a SNP for an individual in reference is missing or not
AdRef=SortedRef;%Create a copy of the reference dataset
AdRef(logical(AbnormalRef))=NaN;%Change the way to mark missing values in the reference dataset 
%[RHO, PVAL] = chi2p(AdRef,3);%Pearson's chi-squared test for a matrix.
% Because this step is very time-comsuming, we pre computed correlation matrix below
load('Sphinx/RHO+PVAL_chi2p_noNaN')% Load the pre-computed correlation matrix
BB_PVAL=(PVAL+PVAL')<ldcutoff;% If correlation coefficient is smaller than a threshold ldcutoff, regard as correlated
BB=BB_PVAL-eye(m);% BB_ij is 1 if SNP i and SNP j is correlated and not the same SNP

% % Exclude the SNP outliers
Selection= true(1,m); % Initialize a Selection vector
n_drop=sum(sortedscore>(6*std(Poolfreq-Reffreq))); % Compute the number of SNPs that are outliers 
       %(i.e. SNPs with too larger utility weight) Normal Distribution within 6 sigma Specification Limits
Selection(1:n_drop)=false; % Drop those SNPs

% % Select a subset of independent published SNPs for the publisher
for iii=1:(m-1) %for each SNPs except the last one
    if Selection(iii)==false % if the SNP has already marked as "no selection"
        continue
    end
    Selection(iii:m)=Selection(iii:m) & not(BB(iii,iii:m));% Mark all SNPs that are correlated with current SNP as "no selection"
end

SelectSortedPool=SortedPool(:,Selection);% Trim the pool according to the Selection vector
SelectSortedRef=SortedRef(:,Selection);% Trim the reference according to the Selection vector
m=size(SelectSortedPool,2) %  # of SNPs in the new pool
SelectSortedPoolfreq=SortedPoolfreq(:,Selection);% Trim the pool-MAF vector according to the Selection vector
SelectSortedReffreq=SortedReffreq(:,Selection);% Trim the reference-MAF vector according to the Selection vector
poolfreq=SelectSortedPoolfreq;% Rename (Create a copy of) the pool-MAF
nonpoolfreq=SelectSortedReffreq;% Rename (Create a copy of) the reference-MAF
SelectSortedID=SortedID(Selection);% Adjust the IDs of remaining SNPs accordingly
selectsortedscore=sortedscore(Selection); % Trim the utility weight vector according to the Selection vector
utility=selectsortedscore;% Rename (Create a copy of) the utility weight vector
sum_utility=sum(utility);% Compute the sum of the utility weight vector

% % % Game Configuration (part 2) - regarding the target set
select_rate=N_T*sample_rate/n_pool;% Calculate the selection rete in the pool
select_rate_ref=N_T*(1-sample_rate)/(n_pop-n_pool);% Calculate the selection rete in the reference
U=N_T*0.05*L*U_factor;% Set (or Calculate) the worth of the data to the data sharer (i.e., his or her maximal benefit)

% % Create the target set
ptID=1:n_pool; %Targets include all individuals in the pool
n=length(ptID); % # of targets in pool
Pool=Pop(ptID,SelectSortedID);% Targets in pool
nptID=(n_pool+1):n_pop; %Targets include all individuals in the reference
n_ref=length(nptID); % # of targets in reference
Ref=Pop(nptID,SelectSortedID);% Targets in reference

% % % Compute LR statistics
lralternate = zeros(n,m); % Initialize the LR statistics in the pool
lralternate_ref = zeros(n_ref,m); % Initialize the LR statistics in the reference
for j = 1:m %j-th SNP
    for i=1:n %i-th individual in pool
        alternatepoolfreq = poolfreq(j); %MAF in pool
        alternatereffreq = nonpoolfreq(j);%MAF in reference
        lralternate(i,j) = singlesnplr (Pool(i,j), alternatepoolfreq, alternatereffreq);%Calculate the LR statistic
    end
    for i=1:n_ref %for the individuals not in the pool
        alternatepoolfreq = poolfreq(j); %MAF in pool
        alternatereffreq = nonpoolfreq(j);%MAF in reference
        lralternate_ref(i,j) = singlesnplr (Ref(i,j), alternatepoolfreq, alternatereffreq);%Calculate the LR statistic
    end
end
lralternate=lralternate';% Transpose
lralternate_ref=lralternate_ref';% Transpose

% % % Start of the Sensitivity analysis % % %
    
% Initialization or recording vectors (or matrix)
A1=[];% Publisher's payoff (Publisher=Sharer)
A2=[];% # of released SNPs
A4=[];% true positive count
A5=[];% false positive count
A8=[];% Adversary's payoff (Adversary=Recipient)
A9=[];% Adversary's estimated payoff
A10=[];% Publisher's benefit
A11=[];% Publisher's cost
AA=[];% Publisher's strategy
%Loop
for k=0:K%For (k+1)-th measure points
    Cp=k/K*G;% Set a value for the cost of expected panalty according to k
    
    % % Initialization of populations in GA algorithm    
    if k==0%
        P2=logical(randi(2,m,n2)-1);% population 2 is randomly generated
    else
        P2(:,1)=P1(:,1);%population 2 is population 1 in the previous run
    end
    P1=true(m,n1);
    for i=1:n1
        P1(end-i+2:end,i)=false(i-1,1);%population 1 is a binary triangular matrix
    end

    % % Iterations in GA algorithm
    for t=1:T
        %To compute F1=Func(P1);
        F1=zeros(n1,1);
        for j=1:n1
            RH_select=P1(:,j);%select j-th column in population 1 as a SNP-sharing strategy
            if sum(RH_select)==0
                F1(j)=0;
            else
                benefit=U*sum(utility(RH_select))/sum_utility;% calculate the benefit
                Selection= RH_select;%Rename (Make a copy)
                if relax2==1 %if the adversary only uses released independent SNPs
                    % Select a subset of independent SNPs for the adversary
                    for iii=1:(m-1)
                        if Selection(iii)==false
                            continue
                        end
                        Selection(iii:m)=Selection(iii:m) & not(BB(iii:m,iii));
                    end
                end
                %For targets in the pool
                sllr=lralternate(Selection,:);%vector of log-likelihood ratios of SNPs to be shared
                lr=exp(sum(sllr,1)');%calculate the likelihood ratio for the sharing strategy
                if Homer_Attack==1 %if the adversary is using the Homer's attack
                    est_dn_tp_i=1./(1+(1-sample_rate)/(sample_rate)./lr);%the true positive count vecotr estimated by the adversary
                else %if the adversary is using the Jordan's attack           
                    est_dn_tp_i=lr*sample_rate;
                    est_dn_tp_i(est_dn_tp_i>1)=1;%the true positive count vecotr estimated by the adversary
                    %est_dn_tp_i=(lr*sample_rate>1)*1+(lr*sample_rate<=1).*lr*sample_rate;%just another way
                end
                est_dn_fp_i=1-est_dn_tp_i;%the false positive count vecotr estimated by the adversary
                if Attacker_Baseline ==0% if the adversary is acting rationally (economically motivated); 
                    attack=(G*est_dn_tp_i-(Cp+Cm))>0; %the attack decision
                else % if the adversary always attacks all individuals; 
                    attack= true(n,1); %the attack decision
                end
                est_dn_fp=attack.*est_dn_fp_i;%the false positive count estimated by the adversary
                est_dn_tn=(~attack).*est_dn_fp_i;%the true negative count estimated by the adversary
                
                F1(j)=benefit-L*sum(attack)*select_rate;%publisher's expected payoff = benfit - expected cost(Bayes' estimation)
                if sum(attack)*select_rate/n>power %if the adversary's detection power (true positive rate) is too large
                    if F1(j)>0
                        F1(j)=-F1(j);% the publisher will tend to not use this j-th sharing strategy
                    end
                end
            end
        end
        
        %To compute F2=Func(P2); exactly same as the way to compute F1 from P1
        F2=zeros(n2,1);
        for j=1:n2 %jth row in P2
            RH_select=P2(:,j);% select j-th column in population 1 as a SNP-sharing strategy
            if sum(RH_select)==0
                F2(j)=0;
            else
                benefit=U*sum(utility(RH_select))/sum_utility;% calculate the benefit
                Selection= logical(RH_select); %Rename (Make a copy)
                if relax2==1 %if the adversary only uses released independent SNPs
                    % Select a subset of independent SNPs for the adversary
                    for iii=1:(m-1)
                        if Selection(iii)==false
                            continue
                        end
                        Selection(iii:m)=Selection(iii:m) & not(BB(iii:m,iii));
                    end
                end
                %For targets in the pool
                sllr=lralternate(Selection,:); %vector of log-likelihood ratios of SNPs to be shared
                lr=exp(sum(sllr,1)');%calculate the likelihood ratio for the sharing strategy
                if Homer_Attack==1 %if the adversary is using the Homer's attack
                    est_dn_tp_i=1./(1+(1-sample_rate)/(sample_rate)./lr); %the true positive count vecotr estimated by the adversary
                else %if the adversary is using the Jordan's attack
                    est_dn_tp_i=lr*sample_rate;
                    est_dn_tp_i(est_dn_tp_i>1)=1; %the true positive count vecotr estimated by the adversary
                    %est_dn_tp_i=(lr*sample_rate>1)*1 + (lr*sample_rate<=1).*lr*sample_rate; %just another way
                end
                est_dn_fp_i=1-est_dn_tp_i; %the false positive count vecotr estimated by the adversary
                if Attacker_Baseline ==0 % if the adversary is acting rationally (economically motivated); 
                    attack=(G*est_dn_tp_i-(Cp+Cm))>0;%the attack decision
                else % if the adversary always attacks all individuals; 
                    attack= true(n,1); %the attack decision
                end
                est_dn_fp=attack.*est_dn_fp_i;%the false positive count estimated by the adversary
                est_dn_tn=(~attack).*est_dn_fp_i;%the true negative count estimated by the adversary
                
                F2(j)=benefit-L*sum(attack)*select_rate;%publisher's expected payoff = benfit - expected cost(Bayes' estimation)
                if sum(attack)*select_rate/n>power %if the adversary's detection power (true positive rate) is too large
                    if F2(j)>0
                        F2(j)=-F2(j); % the publisher will tend to not use this j-th sharing strategy
                    end
                end
            end
        end
        % % Sort
        [a1,b1]=sort(F1,'descend'); % sort F1 in a descending order
        P1=P1(:,b1); % Adjust P1 accordingly
        F1=a1;
        [a2,b2]=sort(F2,'descend'); % sort F2 in a descending order
        P2=P2(:,b2);% Adjust P2 accordingly
        F2=a2;
        % % Immigration
        if mod(t,Interval)==0 % Every interval
            if Direction==1%one-way from P2 to P1
                % Copy the 1st immigration fraction of population 2 to the
                % end immigration fraction of population 1
                P1(:,(n1-Frac*min(n1,n2)+1):n1)=P2(:,1:Frac*min(n1,n2));
                F1((n1-Frac*min(n1,n2)+1):n1)=F2(1:Frac*min(n1,n2));
            elseif Direction==2%two-way
                % Copy the 1st immigration fraction of population 2 to the
                % end immigration fraction of population 1
                P1(:,(n1-Frac*min(n1,n2)+1):n1)=P2(:,1:Frac*min(n1,n2));
                F1((n1-Frac*min(n1,n2)+1):n1)=F2(1:Frac*min(n1,n2));
                % Copy the 1st immigration fraction of population 1 to the
                % end immigration fraction of population 2
                P2(:,(n2-Frac*min(n1,n2)+1):n2)=P1(:,1:Frac*min(n1,n2));
                F2((n2-Frac*min(n1,n2)+1):n2)=F1(1:Frac*min(n1,n2));
            end
            %resort
            [a1,b1]=sort(F1,'descend'); % sort F1 in a descending order
            P1=P1(:,b1); % Adjust P1 accordingly
            F1=a1;
            [a2,b2]=sort(F2,'descend'); % sort F2 in a descending order
            P2=P2(:,b2); % Adjust P2 accordingly
            F2=a2;
        end
        
        if Demo==1 % if in demonstration mode
            % % plot
            subplot(2,2,1)
            imagesc(P1)
            title(['t=',num2str(t)])
            
            subplot(2,2,2)
            plot(F1)
            title([num2str(F1(1)),' > ',num2str(old_A1)])
            
            subplot(2,2,3)
            imagesc(P2)
            title(['Cp=',num2str(Cp)])
            
            subplot(2,2,4)
            plot(F2)
            title([num2str(F2(1)),' > ',num2str(old_A1)])
        end
        
        % % % Evolution process for population 1
        % % Elite
        P1_new=P1; %rename (make a copy)
        a1=(a1+L*n).^2;%scaling function
        pc=[0;cumsum(a1(1:ceil(n1*p_p)))];%production chance
        % % Reproduction
        for i=1:ceil(n1*c_f)
            if Selection_function==1%Roulette
                r1=rand*pc(end);%prob. proportional to Fitness function score
                [~,order1]=sort([pc;r1]);
                s1=find(order1==(ceil(n1*p_p)+2))-1;

                r2=rand*pc(end);
                [~,order2]=sort([pc;r2]);
                s2=find(order2==(ceil(n1*p_p)+2))-1;

            elseif Selection_function==2%Tounament
                r1=randperm(n1,Tournament_size);
                r2=randperm(n2,Tournament_size);
                s1=min(r1);
                s2=min(r2);
            end
            %Crossover
            if Crossover_function==0%half-half
                col1=[ones(ceil(m/2),1),zeros(floor(m/2),1)];
            elseif Crossover_function==1%one point
                cut=randi(m-1);%uniform distributed
                col1=[ones(cut,1),zeros(m-cut,1)];
            elseif Crossover_function==2%two points
                cut=randperm(m,2);
                col1=[ones(min(cut),1),zeros(max(cut)-min(cut),1),ones(m-max(cut),1)];
            elseif Crossover_function==3%scatters
                col1=randi(2,m,1)-1;
            end
            child=logical(P1(:,s1).*col1+P1(:,s2).*~col1);% generate children
            P1_new(:,floor(n1*p_E)+i)=child;% Keep the elite and save the children
        end
        % % mutation
        for i=1:(n1-floor(n1*p_E)-ceil(n1*c_f))
            r1=randi(n1);
            Mutant=P1(:,r1);
            r2=randi(ceil(m*p_m));
            pr=randperm(m,r2);
            Mutant(pr)=~Mutant(pr);
            P1_new(:,ceil(n1*c_f)+floor(n1*p_E)+i)=Mutant;% Save the mutant
        end
        P1=P1_new;%Save the new generation
        
        % % % Evolution process for population 2 (exactly same as population 1)
        % % Elite
        P2_new=P2; %rename (make a copy)
        a2=((a2+L*n));%scaling function
        pc=[0;cumsum(a2(1:ceil(n2*p_p)))];%production chance
        % % Reproduction
        for i=1:ceil(n2*c_f)
            if Selection_function==1%Roulette
                r1=rand*pc(end);%prob. proportional to Fitness function score
                [~,order1]=sort([pc;r1]);
                s1=find(order1==(ceil(n2*p_p)+2))-1;

                r2=rand*pc(end);
                [~,order2]=sort([pc;r2]);
                s2=find(order2==(ceil(n2*p_p)+2))-1;
            elseif Selection_function==2%Tounament
                r1=randperm(n1,Tournament_size);
                r2=randperm(n2,Tournament_size);
                s1=min(r1);
                s2=min(r2);
            end
            %Crossover
            if Crossover_function==0%half-half
                col1=[ones(ceil(m/2),1),zeros(floor(m/2),1)];
            elseif Crossover_function==1%one point
                cut=randi(m-1);%uniform distributed
                col1=[ones(cut,1),zeros(m-cut,1)];
            elseif Crossover_function==2%two points
                cut=randperm(m,2);
                col1=[ones(min(cut),1),zeros(max(cut)-min(cut),1),ones(m-max(cut),1)];
            elseif Crossover_function==3%scatters
                col1=randi(2,m,1)-1;
            end
            child=logical(P2(:,s1).*col1+P2(:,s2).*~col1); % generate children
            P2_new(:,floor(n2*p_E)+i)=child; % Keep the elite and save the children
        end
        % % mutation
        for i=1:(n2-floor(n2*p_E)-ceil(n2*c_f))
            r1=randi(ceil(n2*c_f)+floor(n2*p_E));
            Mutant=P2_new(:,r1);
            r2=randi(ceil(m*p_m));
            pr=randperm(m,r2);
            Mutant(pr)=~Mutant(pr);
            P2_new(:,ceil(n2*c_f)+floor(n2*p_E)+i)=Mutant; % Save the mutant
        end
        P2=P2_new; %Save the new generation
    end
    
    A1=[A1,F1(1)];% Save the publisher's payoff (Choose the 1st strategy from population 1)
    A2=[A2,sum(P1(:,1))];% Save the # released SNPs
    AA=[AA,P1(:,1)];% Save the publisher's strategy
    if sum(P1(:,1))==0% if sharer releases nothing
        A4=[A4,0];% save the true positive count
        A5=[A5,0];% save the false positive count
        A8=[A8,0];% save the adversary's payoff
        A9=[A9,0];% save the adversary's estimated payoff
        A10=[A10,0];% save the publisher's benefit
        A11=[A11,0];% save the publisher's cost
    else
        % % compute the adversary's payoff given the publisher's strategy
        benefit=U*sum(utility(P1(:,1)))/sum_utility; % calculate the benefit
        Selection= logical(P1(:,1)); % Make a copy of the publisher's strategy
        if relax2==1 %if the adversary only uses released independent SNPs
            for iii=1:(m-1)
                if Selection(iii)==false
                    continue
                end
                Selection(iii:m)=Selection(iii:m) & not(BB(iii:m,iii));
            end
        end
        %For targets in the pool
        sllr=lralternate(Selection,:); %vector of log-likelihood ratios of SNPs to be shared
        lr=exp(sum(sllr,1)'); %calculate the likelihood ratio for the sharing strategy
        if Homer_Attack==1 %if the adversary is using the Homer's attack
            est_dn_tp_i=1./(1+(1-sample_rate)/(sample_rate)./lr);
        else %if the adversary is using the Jordan's attack
            est_dn_tp_i=lr*sample_rate;
            est_dn_tp_i(est_dn_tp_i>1)=1; %the true positive count vecotr estimated by the adversary
        end
        est_dn_fp_i=1-est_dn_tp_i; %the false positive count vecotr estimated by the adversary
        if Attacker_Baseline ==0 % if the adversary is acting rationally (economically motivated); 
            attack=(G*est_dn_tp_i-(Cp+Cm))>0; %the attack decision
        else % if the adversary always attacks all individuals;
            attack= true(n,1); %the attack decision
        end
        est_dn_fp=attack.*est_dn_fp_i; %the false positive count estimated by the adversary
        est_dn_tn=(~attack).*est_dn_fp_i; %the true negative count estimated by the adversary
        
        %For targets in the reference
        sllr_ref=lralternate_ref(Selection,:); %vector of log-likelihood ratios of SNPs to be shared
        lr_ref=exp(sum(sllr_ref,1)'); %calculate the likelihood ratio for the sharing strategy
            if Homer_Attack==1 %if the adversary is using the Homer's attack
                est_dn_tp_i_ref=1./(1+(1-sample_rate)/(sample_rate)./lr_ref);
            else %if the adversary is using the Jordan's attack
                est_dn_tp_i_ref=lr_ref*sample_rate;
                est_dn_tp_i_ref(est_dn_tp_i_ref>1)=1; %the true positive count vecotr estimated by the adversary
            end
        est_dn_fp_i_ref=1-est_dn_tp_i_ref; %the false positive count vecotr estimated by the adversary
        if Attacker_Baseline ==0 % if the adversary is acting rationally (economically motivated); 
            attack_ref=(G*est_dn_tp_i_ref-(Cp+Cm))>0; %the attack decision
        else % if the adversary always attacks all individuals;
            attack_ref=true(n_ref,1); %the attack decision
        end
        est_dn_fp_ref=attack_ref.*est_dn_fp_i_ref; %the false positive count estimated by the adversary
        est_dn_tn_ref=(~attack_ref).*est_dn_fp_i_ref; %the true negative count estimated by the adversary
        
        A4=[A4,sum(attack)*select_rate];% save expected true positive count (Bayes' estimation)
        A5=[A5,sum(attack_ref)*select_rate_ref]; % save expected false positive count (Bayes' estimation)
        A8=[A8,(G-(Cp+Cm))*sum(attack)*select_rate-(Cp+Cm)*sum(attack_ref)*select_rate_ref];% save the adversary's expected payoff (Bayes' estimation)
        A9=[A9,sum(attack.*(G*est_dn_tp_i-(Cp+Cm)))*select_rate+sum(attack_ref.*(G*est_dn_tp_i_ref-(Cp+Cm)))*select_rate_ref];% save teh adversary's expected estimated payoff (Bayes' estimation)
        A10=[A10,benefit]; % save the publisher's benefit
        A11=[A11,L*sum(attack)*select_rate]; % save the publisher's expected cost (Bayes' estimation)
    end
end
toc % read the stopwatch timer

% % Plot figures of results
xx=(0:K)/K*G; % x-axis of measured point

figure
plot(xx,A2,'-*b')
xlabel('Expected cost of penalty to the recipient')
ylabel('Optimal # of Released SNPs')

figure
plot(xx,A2/m,'-*b')
xlabel('Expected cost of penalty to the recipient')
ylabel('Optimal portion of released SNPs')
axis([0 xx 0 1])

figure
plot(xx,A1,'-*b',xx,A8,'-*r',xx,A9,'-*m')
xlabel('Expected cost of penalty')
legend('Publisher''s maximal payoff','Adversary''s payoff','Adversary''s estimated payoff','Location','NorthWest')

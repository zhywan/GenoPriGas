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
close all
clear
format compact
rng(0);
tic

%--Sensitivity analysis on the expected cost of penalty to the adversary--%


% % Game Configuration
Homer_Attack=0; % 1=Homer's attack; 0=Jordan's attack (default)
relax2=0; % 1=Adversary only use released independent SNPs; 0=Adversry use all released SNPs (default)
Demo=0; % 1=Show the GA running process; 0=Hide the GA running process (default)
N_Class=2; % 1=1-class classification (LR test: in pool or in population); 
           % 2=2-classes classification (LR test: in pool or in reference) (default)
Attacker_Baseline=0;% 1=adversary will attack all individuals; 
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

% % Handling missing values
MissingSNP=sum(Pool==-1);%A vector, in which each element is # of individuals in pool with missing values for a SNP
RetainedSNP=MissingSNP<=n_pool*theta_miss;%A binary vector, in which each element indicates a SNP will be retained or not
RetainedID=1:size(Pool,2);%Create a ID for each individual in pool
Pool=Pool(:,RetainedSNP);%Discard all the SNPs that should not be retained from the pool
Ref=Ref(:,RetainedSNP);%Discard all the SNPs that should not be retained from the reference
RetainedID=RetainedID(RetainedSNP);%Record the IDs of those retained SNPs

% % Compute MAF
MissPool=(Pool==-1);%A binary matrix, in which each element indicates a SNP for an individual in pool is missed or not
MissRef=(Ref==-1);%A binary matrix, in which each element indicates a SNP for an individual in reference is missed or not
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
m=size(RetainedPool,2);% the # of SNPs in the new pool
RetainedScore=Score(~Reject);% Discard SNPs with too small MAF in the utilty weight vector
RetainedRef=Ref(:,~Reject);% Discard SNPs with too small MAF freom the reference
RetainedPoolfreq=Poolfreq(:,~Reject);% Discard SNPs with too small MAF from the poof MAF vector
RetainedReffreq=Reffreq(:,~Reject);% Discard SNPs with too small MAF from the reference MAF vector
RetainedID=RetainedID(~Reject);%Record the IDs of remaining SNPs

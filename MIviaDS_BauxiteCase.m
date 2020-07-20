%----------------------------------------------------------------------------------------------------------------
% Authors: Hamed Mohammadi, Sajjad Talesh Hosseini, Omid Asghari
% Paper: "A new approach for multivariate imputation of unequally sampled compositional variables 
%         in the 3D datasets via a direct sampling-based algorithm (with three case studies)"
% Journal: Computer & Geosciences
%----------------------------------------------------------------------------------------------------------------

clc
clear

%%%--------------------------
% Simulation parameters
%%%--------------------------
Validation=1;                           % 1= Use of code to validation with a refrence data file
n_cat = 1;                              % number of categorical variables                             
n = 12;                                 % maximum number of points in the data event
f = 1;                                  % maximum fraction of the training data to be scanned
t1 = 0.2;                               % acceptance treshold for multivariate distance between central nodes (between 0 and 1)
t2 = 0.02;                              % acceptance threshold for multivariate distance between Data events
nb_of_real=1;                           % number of realizations
expo_w=1;                               % exponent of the inverse distance weights
search_radius_z=20;                     % Max search radius in vertical direction (in direction of drill-hole data),compute in data coordinate units (meter)
search_radius_xy=20;                    % Max search radius in horizontal direction, (meter)
z_comp_length=1;                        % Sample interval length (meter)

% loading the unequally sample data file:
% In this file missed valued must be set to -999
% first column are the id of each sample, 2nd,3rd,4th Columns are x,y,z coordinate of each samles
% 5th to end column filles by continuous variables followed by the categorical variables

Data=load ('Missing_Data.txt'); Data=sortrows(Data);
nb_of_var=size(Data,2)-4; %Number of variables (contionuous + categorical)
Data(Data==-999)=NaN;         
dmax=max(Data(:,5:end))-min(Data(:,5:end)); %compute range for each variable

% For Validation
% If the number of variables was diffren from this example(Bauxite Case), these section must be modified by user
if Validation==1
id_of_removed_6=find(isnan(Data(:,6)));
id_of_removed_7=find(isnan(Data(:,7)));
id_of_removed_8=find(isnan(Data(:,8)));
end

% Divide input sample files (Data) into Sim(Unequally sampled) and TD(Equally Sampled as Traning Data)
A=sum(isnan(Data),2);
raw_Sim=find(A>0);raw_TD=find(A==0);
Sim=Data(raw_Sim,:);TD=Data(raw_TD,:);
% 3D Plot of Sim(Unequally sampled) and TD(Equally Sampled)
figure(1)
scatter3(TD(:,2),TD(:,3),TD(:,4),5,'b','filled'); hold on
scatter3(Sim(:,2),Sim(:,3),Sim(:,4),5,'r','filled');
axis equal tight;
view(40,35)

%%

tic
% Simulation loop
for i0=1:nb_of_real
    Neighbors=TD; %Neighbors is a file contain conditional samples, at first it is equal to TD and each sample that its simulation was completed added to it
    
    %define a random path of Sim Samples
    path_sim = randperm(size(Sim,1));
    for i1=1:size(Sim,1)
    %Coordinate of Simulation Sample
    Sample=Sim(path_sim(i1),:);
    %Find spatial distance with all sample in TD file dxy=Lataral Distance,
    %dz=Vertical Distance
    dxy = sqrt(((Sample(3)-Neighbors(:,3)).^2)+((Sample(2)-Neighbors(:,2)).^2));
    dz=abs((Sample(4)-Neighbors(:,4)));
    dz1=(Sample(4)-Neighbors(:,4));
    %sorting
    dz((find(dxy>=search_radius_xy)))=inf;
    [dz,s] = sort(dz);
    
    %checking how many points are within outer search radius
    nb_in_search_radius = sum(dz<=search_radius_z);
    
    if nb_in_search_radius < n
       informed_nodes_counter = nb_in_search_radius;
    else
        informed_nodes_counter = n;
    end
    
    % Agar hich hamsayegi dar r teaain shode peida nashod boro ye nemone
    % dige peydakon
    
    if nb_in_search_radius==0
       Neighbors=vertcat(Neighbors,Sample); % This sample wasnt simulated!
        continue
    end
    
    %Compute the z vector to informed nodes by a tolerance of half
    %composite length
    vect_to_Z=dz1(s(1:informed_nodes_counter));
    wd=abs(vect_to_Z);
    vect_to_Z_plus1=vect_to_Z+(z_comp_length/2);
    vect_to_Z_minus1=vect_to_Z-(z_comp_length/2);
    
    Sim_data_event=zeros(informed_nodes_counter,nb_of_var);
    TD_data_event=zeros(informed_nodes_counter,nb_of_var);
    
    for i3=1:informed_nodes_counter
        Sim_data_event(i3,:)=Neighbors(s(i3),5:end);  
    end
    
    missed_var_insample=find(isnan(Sample)); % Which variable are missed in Simulation Sample?
    
    %Simulation of each missed variable in a random path
    path_sim_var = randperm(numel(missed_var_insample)); %Defininton a random path of missed variable in simulation path
    for i2=1:numel(missed_var_insample)
        
        Central_SimE=Sample(5:end);
        Central_missed_var=find(isnan(Central_SimE));
        Central_SimE(Central_missed_var)=[];
        
        %
        distances=inf(size(TD,1),1);
        TD_central_values=nan(size(TD,1),1);
        
        % Scanning TD in a Random Path
        path_td = randperm(size(TD,1));
        
        for i4=1:f*(size(TD,1))
            
            xTD=TD(path_td(i4),2);
            yTD=TD(path_td(i4),3);
            zTD=TD(path_td(i4),4);
            Central_TDE=TD(path_td(i4),5:end);
            Central_TDE(Central_missed_var)=[];
            
 %compute missmatch between variable in simulation sample and condidated TD Sample

            D1_Cat=Central_TDE(end-n_cat+1:end)==Central_SimE(end-n_cat+1:end);
            D1_Con=(Central_TDE(:,1:end-n_cat)-Central_SimE(:,1:end-n_cat)).^2;
            dmax1=dmax;dmax1(Central_missed_var)=[];
            dmax1=repmat(dmax1(:,1:end-n_cat).^2,size(D1_Con,1),1);
            D1_Con=D1_Con./dmax1;
            D1_Con=sum(D1_Con);
            D1_Con=sqrt(D1_Con);           
% If, mismatch between Rock type (D1_Cat) of these samples or Disatance
% between contiuous variable is greter than t1 go to next sample of TD
            if D1_Cat==0 || D1_Con>t1
                continue
            end
             
            
            TD_dxy=sqrt(((xTD-TD(:,2)).^2)+((yTD-TD(:,3)).^2));
            %Taeine TD_data_event
            Same_bh_nodes=TD(TD_dxy<=search_radius_xy,:);
            ZsTD_Up=zTD+vect_to_Z_plus1;
            ZsTD_blw=zTD+vect_to_Z_minus1;
            
            for i5=1:size(Sim_data_event,1)
                row=find((Same_bh_nodes(:,4)>=ZsTD_blw(i5))&(Same_bh_nodes(:,4)<=ZsTD_Up(i5)));
                if isempty(row)
                    TD_data_event(i5,:)=nan;
                else
                    TD_data_event(i5,:)=Same_bh_nodes(row(1),5:end);
                end
            end
            
            
            if sum(~isnan(TD_data_event(:,1)))>=2
                nan_nodes=find(isnan(TD_data_event(:,1)));
                TD_data_event(nan_nodes,:)=[]; Sim_data_event(nan_nodes,:)=[];
                wd(nan_nodes)=[];
            else
                continue
            end
            
            %Compute Multivariate Euqlidisian Distance Between Data Events
            D2_Cat=(TD_data_event(end-n_cat+1:end)~=Sim_data_event(end-n_cat+1:end)); D2_Cat=mean(mean(D2_Cat));
            
            D2_Con=(TD_data_event(:,1:end-n_cat)-Sim_data_event(:,1:end-n_cat)).^2;
            dmax1=repmat(dmax(:,1:end-n_cat).^2,size(D2_Con,1),1);
            D2_Con=D2_Con./dmax1;
            wd=(wd.^expo_w)./sum(wd.^expo_w);
            D2_Con=D2_Con .* wd;
            D2_Con=sum(D2_Con);
            D2_Con=sqrt(D2_Con);
            D2_Con=mean(D2_Con);
            Distance1=(D2_Cat+D2_Con)./2;
            
            % Save distance and the condidate TD sample values, If find a
            % value less than t2, value of TD sample directly pasted to
            % Sample(i2)
            %if no distance less than t2, smallest distance was chosen
            distances(path_td(i4))=Distance1;
            TD_central_values(path_td(i4))=TD(path_td(i4),missed_var_insample(path_sim_var(i2)));
            
            if Distance1 <= t2
                break
            end
            
        end
        
        [distances,s1]=sort(distances);
        Sample(missed_var_insample(path_sim_var(i2)))=TD_central_values(s1(1));
    end
    %Simulation of sample was completed and assigned it to Neighbors
    Neighbors=vertcat(Neighbors,Sample);
    
    end
% Complete Realization (Reconstructed Data and Save it
Reconstructed_Data=sortrows(Neighbors);
Name_of_output=i0;
xlswrite(num2str(i0),Reconstructed_Data);

not_sim=sum(sum(isnan(Reconstructed_Data)));
disp(not_sim);%display number of samples that that remain Unequal in each realization
no_sim(i0)=not_sim; %Save number of samples that that remain Unequal in each realization

if Validation==1
    Imputed_6=Reconstructed_Data(:,6);Imputed_6=Imputed_6(id_of_removed_6);
    Imputed_66(:,i0)=Imputed_6;
    
    Imputed_7=Reconstructed_Data(:,7);Imputed_7=Imputed_7(id_of_removed_7);
    Imputed_77(:,i0)=Imputed_7;
    
    Imputed_8=Reconstructed_Data(:,8);Imputed_8=Imputed_8(id_of_removed_8);
    Imputed_88(:,i0)=Imputed_8;
end

end
toc
%%
% Validation Plots 
if Validation==1

Data=sortrows(Data);
id_of_removed_6=find(isnan(Data(:,6)));
id_of_removed_7=find(isnan(Data(:,7)));
id_of_removed_8=find(isnan(Data(:,8)));

Org_Data=load('Refrence_Data.txt');
Org_Data=sortrows(Org_Data);
True_6=Org_Data(:,6);True_6=True_6(id_of_removed_6);
True_7=Org_Data(:,7);True_7=True_7(id_of_removed_7);
True_8=Org_Data(:,8);True_8=True_8(id_of_removed_8);

e_type_6=mean(Imputed_66,2);
e_type_7=mean(Imputed_77,2);
e_type_8=mean(Imputed_88,2);

figure(2); hold on
subplot(1,5,2)
scatter(e_type_6,True_6,5); axis equal tight
refline(1,0)
subplot(1,5,3)
scatter(e_type_7,True_7,5); axis equal tight
refline(1,0)
subplot(1,5,4)
scatter(e_type_8,True_8,5); axis equal tight
refline(1,0)

figure(4)
subplot(1,5,2);hold on
for i=1:nb_of_real
ecdf(Imputed_66(:,i))
end
ecdf(True_6)
subplot(1,5,3);hold on
for i=1:nb_of_real
ecdf(Imputed_77(:,i))
end
ecdf(True_7)
subplot(1,5,4);hold on
for i=1:nb_of_real
ecdf(Imputed_88(:,i))
end
ecdf(True_8)

figure(5)
subplot(1,5,2);hold on
x6=Org_Data(:,2);x6=x6(id_of_removed_6);
y6=Org_Data(:,3);y6=y6(id_of_removed_6);
z6=Org_Data(:,4);z6=z6(id_of_removed_6);
for i=1:nb_of_real
d = variogram([x6 y6 z6],Imputed_66(:,i),'plot',true,'nrbins',30);
end
d = variogram([x6 y6 z6],True_6,'plot',true,'nrbins',30);

subplot(1,5,3);hold on
x7=Org_Data(:,2);x7=x7(id_of_removed_7);
y7=Org_Data(:,3);y7=y7(id_of_removed_7);
z7=Org_Data(:,4);z7=z7(id_of_removed_7);
for i=1:nb_of_real
d = variogram([x7 y7 z7],Imputed_77(:,i),'plot',true,'nrbins',30);
end
d = variogram([x7 y7 z7],True_7,'plot',true,'nrbins',30);

x8=Org_Data(:,2);x8=x8(id_of_removed_8);
y8=Org_Data(:,3);y8=y8(id_of_removed_8);
z8=Org_Data(:,4);z8=z8(id_of_removed_8);
subplot(1,5,4);hold on
for i=1:nb_of_real
d = variogram([x8 y8 z8],Imputed_88(:,i),'plot',true,'nrbins',30);
end
d = variogram([x8 y8 z8],True_8,'plot',true,'nrbins',30);
end

[r,c]=find(isnan(sum(Data,2)));
T6=Org_Data(r,6);
T7=Org_Data(r,7);
T8=Org_Data(r,8);
T9=Org_Data(r,9);
figure(6)
subplot(2,3,1)
scatter(T6,T7,5,T9);
title('SiO2(X)-Fe2O3(Y)')
subplot(2,3,2)
scatter(T6,T8,5,T9);
title('SiO2(X)-TiO2(Y)')
subplot(2,3,3)
scatter(T7,T8,5,T9);
title('Fe2O3(X)-TiO2(Y)')

Reconstructed_Etype=Org_Data;
Reconstructed_Etype(id_of_removed_6,6)=e_type_6;
Reconstructed_Etype(id_of_removed_7,7)=e_type_7;
Reconstructed_Etype(id_of_removed_8,8)=e_type_8;
I6=Reconstructed_Etype(r,6);
I7=Reconstructed_Etype(r,7);
I8=Reconstructed_Etype(r,8);
I9=Reconstructed_Etype(r,9);
subplot(2,3,4)
scatter(I6,I7,5,I9);
subplot(2,3,5)
scatter(I6,I8,5,I9);
subplot(2,3,6)
scatter(I7,I8,5,I9);

clc
clear all
close all
n_ext=6;
set=1:n_ext;
%% 2-neighborhood
mat1=eye(n_ext,n_ext);
mat2=-eye(n_ext,n_ext);
N_2=[mat1;mat2];
num_vecinos_N_2=2*n_ext;
%% Lflat neighborhood 
N_Lflat=[(2^((n_ext)+1))-2,n_ext];

k=1;
 for i=0:n_ext-1
    sub=nchoosek(set,n_ext-i);
    [f,~] = size(sub);
    for j=1:f 
       N_Lflat(k,1:n_ext)=ismember(set,sub(j,:));
       k=k+1;
       N_Lflat(k,1:n_ext)=-ismember(set,sub(j,:)); 
       k=k+1; 
    end
end
N_Lflat;
num_vecinos_N_Lflat=(2^((n_ext)+1))-2;

%% Mflat neighborhood
mat1=[mat1;zeros(1,n_ext)];
[f,~] = size(mat1);
k=1;
for i=1:f
   for j=1:f
       if i~=j
          N_Mflat(k,1:n_ext)=mat1(i,:)-mat1(j,:); 
          k=k+1;      
       end
   end
   
end
N_Mflat;
num_vecinos_N_Mflat=n_ext*(n_ext+1);

%% Infinity neighborhood
N_Infinity=npermutek([-1 0 1],n_ext);
N_Infinity = N_Infinity(any(N_Infinity,2),:);
num_vecinos_N_Infty=(3^n_ext)-1;

%% Export to txt file example
vec=1:num_vecinos_N_Infty';
matt=N_Infinity;
fileID=fopen(strcat(num2str(n_ext),'_','Infinity.csv'),'w');
for i=1:length(vec)
fprintf(fileID,strcat('d',num2str(vec(i)),',',regexprep(num2str(matt(i,:)),'\s+',','),'\n'));
end
fclose(fileID);

vec=1:num_vecinos_N_2';
matt=N_2;
fileID=fopen(strcat(num2str(n_ext),'_','Separable.csv'),'w');
for i=1:length(vec)
fprintf(fileID,strcat('d',num2str(vec(i)),',',regexprep(num2str(matt(i,:)),'\s+',','),'\n'));
end
fclose(fileID);

vec=1:num_vecinos_N_Lflat';
matt=N_Lflat;
fileID=fopen(strcat(num2str(n_ext),'_','Lflat.csv'),'w');
for i=1:length(vec)
fprintf(fileID,strcat('d',num2str(vec(i)),',',regexprep(num2str(matt(i,:)),'\s+',','),'\n'));
end
fclose(fileID);

vec=1:num_vecinos_N_Mflat';
matt=N_Mflat;
fileID=fopen(strcat(num2str(n_ext),'_','Mflat.csv'),'w');
for i=1:length(vec)
fprintf(fileID,strcat('d',num2str(vec(i)),',',regexprep(num2str(matt(i,:)),'\s+',','),'\n'));
end
fclose(fileID);


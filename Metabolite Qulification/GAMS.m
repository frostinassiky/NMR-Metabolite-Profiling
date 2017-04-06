%% WriteGAMS_data
fileID = fopen('GAMS_data.txt','w');

% empty line
fprintf(fileID,'\n');

% i, j, m
fprintf(fileID,'set i /1*%d/;\nset j /1*%d/;\nset m /1*%d/;\n',...
    length(peak_part_th),length(peak_part_s),settings.num);
fprintf(fileID,'\n');

% S
fprintf(fileID,'set S(i,j);  S(i,j)=no;\n');
[row,col,v] = find(pstrength);
for k = 1:length(v)
    fprintf(fileID,'S(''%d'',''%d'')=yes;\n',row(k),col(k));    
end
fprintf(fileID,'\n');

% K
fprintf(fileID,'set K(m,i);  K(m,i)=no;\n');
for k = 1:length(peak_part_th)
    fprintf(fileID,'K(''%d'',''%d'')=yes;\n',peak_index(peak_part_th(k)),k);    
end
fprintf(fileID,'\n');

% parameter c2
fprintf(fileID,'parameter c2(i,j);\n');
for k = 1:length(v)
    loc1 = locs_th( peak_part_th(row(k)) );
    loc2 = locs_s( peak_part_s(col(k)) );
    c2 =  gaussmf( loc1,[gauss_p1, loc2]);
    fprintf(fileID,'c2(''%d'',''%d'')=%e;\n',row(k),col(k),c2);    
end
fprintf(fileID,'\n');

% parameter pks_sample(j);
fprintf(fileID,'parameter pks_sample(j);\n');
for k = 1:length(peak_part_s)
    fprintf(fileID,'pks_sample(''%d'')= %e;\n',k,pks_s(peak_part_s(k)));    
end
fprintf(fileID,'\n');

% parameter pks_theory(i);
fprintf(fileID,'parameter pks_theory(i);\n');
for k = 1:length(peak_part_th)
    fprintf(fileID,'pks_theory(''%d'')= %e;\n',k,pks_th(peak_part_th(k)));    
end
fprintf(fileID,'\n');

% parameters 
fprintf(fileID,'parameter sigma;\n');
fprintf(fileID,'sigma = %e;\n',gauss_p2 );
fprintf(fileID,'\n');
fprintf(fileID,'scalar maxConcentration /9.0/;\n');
fprintf(fileID,'\n');

fclose(fileID);

%% Run GAMS
gamsdir = '"C:\GAMS\win64\23.9\gams"';  
jobdir= pwd ;
roFileName = 'find_match.gms';
r = system(['cd ' jobdir ' & ' gamsdir ' ' roFileName ' lo=3 > roLog.txt']);
if r ~= 0 && r ~= 112
    msgbox('error! PLS check files.');
    error(['error: ' , num2str(r)])
end

%% ReadGAMS_output
fileID = fopen('GAMS_output.txt');
GAMS_output = struct('i',[],'j',[],'m',[],'c',[]);
file = textscan(fileID,'%s','delimiter','#'); 
k1 = 1; k2 = 1; index=0;
for k = 1 : length(file{1})
    if isempty(file{1}{k})
        index = index+1;
        continue
    end
    switch index
        case 1
            ites = textscan(file{1}{k},'%d,%d');
            GAMS_output.i(k1)=ites{1};
            GAMS_output.j(k1)=ites{2};
            k1 = k1 + 1;
        case 2
            ites = textscan(file{1}{k},'%d,%f');
            GAMS_output.m(k2)=ites{1};
            GAMS_output.c(k2)=ites{2};
            k2 = k2 +1;
    end
end
fclose(fileID);     


%データセーブ

mkdir("design_data")
cd("design_data")
mkdir(pj_name)
cd(pj_name)

pj_variables=strcat(pj_name,'_variables.txt')
save(pj_variables)

cd(now_work)
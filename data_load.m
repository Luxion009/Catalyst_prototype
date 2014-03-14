%データロード

now_work=pwd;
cont=input("データを読み出さず設計を開始しますか？ yes:1 no:0 ---");

if(cont!=1)
	clear all
	new=input('機体データを新規作成しますか? yes:1 no:0 ---');
	if(new==1)
		now_work=pwd;
		pj_name=input("機体名を入力してください ---",['s']);
		cd("design_data")
		mkdir(pj_name);
		cd(pj_name)

		pj_variables=strcat(pj_name,'_variables.txt')
		save(pj_variables)

		cd(now_work)

	else
		pj_name=input("ロードする機体名を入力してください ---",['s']);
		cd("design_data")
		cd(pj_name)

		l_data=strcat(pj_name,'_variables.txt');
		load(l_data)

		cd(now_work)

	endif
endif
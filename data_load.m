%�f�[�^���[�h

now_work=pwd;
cont=input("�f�[�^��ǂݏo�����݌v���J�n���܂����H yes:1 no:0 ---");

if(cont!=1)
	clear all
	new=input('�@�̃f�[�^��V�K�쐬���܂���? yes:1 no:0 ---');
	if(new==1)
		now_work=pwd;
		pj_name=input("�@�̖�����͂��Ă������� ---",['s']);
		cd("design_data")
		mkdir(pj_name);
		cd(pj_name)

		pj_variables=strcat(pj_name,'_variables.txt')
		save(pj_variables)

		cd(now_work)

	else
		pj_name=input("���[�h����@�̖�����͂��Ă������� ---",['s']);
		cd("design_data")
		cd(pj_name)

		l_data=strcat(pj_name,'_variables.txt');
		load(l_data)

		cd(now_work)

	endif
endif
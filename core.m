%�������Ђ�
%WASA14�N�݌v��C�@�Ȃہ[
%twitter:Luxion009

printf("\n");
disp("�l�͔�s�@�����݌v�v���O����")
disp("---Catalyst prototype---");
printf("\n");
disp("");

data_load

do
	printf("[1]��͉��\n[2]���^�쐬\n[3]���ʌ^�쐬\n[4]��͉��\n[5]�f�[�^�ۊ�\n[6]���݌v\n[6.5]���Ē�`\n[7]�������[�����g��`\n[7.5]������W���v�Z\n[8]�^�����\n[9]�f�[�^�o��\n[10]�I��\n");
	mode=input("���[�h�Z���N�g ---");

	if(mode==1)
		%�z�����߂�
		aerodesign

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=2;
		endif

	endif

	if(mode==2)

		%���^�쐬��⏕����
		defprofilepal

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif

		% %�쐬�����������f�[�^����e�����悭�[���^�𐶐�����
		% defprofile]

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=3;
		endif

	endif

	if(mode==3)

		%���ʌ`�����
		defplane2

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=4;
		endif
		
	endif

	if(mode==4)
		%��͉��
		wing_analysis

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif
		
		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=5;
		endif

	endif

	if(mode==5)
		% %���W�_�f�[�^�擾
		% foil_p_input
		% %������`
		% defthn
		%�}���h�����݌v
		defmand
		%�f�[�^���
		filldata

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif
		
		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=6;
		endif
		
	endif

	if(mode==6)
		
		%�ޗ��v�Z
		defspar

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif
		
		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=7;
		endif
		
	endif

	if(mode==6.5)

		%�����I�ޗ��v�Z
		plyrefine
		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif
		
		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=7;
		endif
	endif

	if(mode==7)

		%�������[�����g��`
		definertia

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=7.5;
		endif
	endif

	if(mode==7.5)
		%�^����͂ɕK�v�Ȕ��W���̐��Z
		%�㔼�p����
		defclbeta
		%���[���p���x
		def_p
		%���[�p���x
		def_r

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=8;
		endif
	endif


	if(mode==8)
		
		%�e�[����`
		deftail

		%����������`
		deftwing

		%����������`
		defvwing

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif
		
	endif

	if(mode==9)

	

	endif

until(mode==10)
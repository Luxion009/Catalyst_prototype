%�^����͂��s������������


%�A�X�y�N�g������߂�
aspect_ratio=(span+wl_span*2)^2/mw_s;

%�嗃��͕��ϗ��������߂�
cbar=2/mw_s*sum(f_data(:,2).^2.*line_e_d(:,1));

%��͕��ϗ����̈ʒu�i�嗃�O�������_�Ƃ�����W�n�j�����߂�
cbar_pos=interp1(f_data(51:250,2),f_data(51:250,1),cbar,'linear');%��Ԋ֐��͋��`�P�������łȂ��Ă͂Ȃ�Ȃ�
%��͕��ϗ����̌��ʒu(�d�S�ʒu)�����߂�
cbar_sparpos=interp1(f_data(:,1),spar_pos(:,1),cbar_pos,'linear');
cog=cbar_sparpos;

%��͕��ϗ����ɂ����錅�ʒu�����͒��S(0.25)�܂ł̋���(m)�����߂�
cbar_ac=cbar*(cbar_sparpos-0.25);

%�嗃�g�͌X�΂����߂�(0.5�x�ŉ��)
%L_cla=clalpha_analysis_m(f_data,Q_ij,mt_ang,line_e_d,lst_re,lst_i_ang,lst_comp_vel,d_vel,wld_c,we_div,a_density,adata)(1,1);

cla_analysis

alphaw=(cla_L_lst*9.8-L_lst*9.8)*2/a_density/d_vel^2/mw_s*180/pi/0.5;

%�����ʒu�ł̐������낵�X�΂����߂�(�ȉ~���Ƌߎ�����)
% dep_dal=(mean(lst_i_ang(1:we_div/5)).*2.-mean(cla_lst_i_ang(1:we_div/5)).*2)*180/pi/0.5;
dep_dal=2*alphaw/pi/aspect_ratio;

% %����������̓f�[�^�ǂݍ���
% cd("NACA0012")

% %-----���^�̓Ǎ�
% fp=fopen('hsfoil.dat');
% hsfoil=fgetl(fp);
% hsfoil=(fscanf(fp,'%f',[2,250]))';
% fclose(fp);
% fshs=0;


% %-----������͉�̓f�[�^��
% %-----���C�m���Y���X�g�ǂݍ���
% Relist={50000,'0.050';100000,'0.100'; 150000,'0.150'; 200000,'0.200'; 250000,'0.250'; 300000,'0.300';350000,'0.350'; 400000,'0.400';450000,'0.450'; 500000,'0.500';550000,'0.550';600000,'0.600';};;

% for m=1:12
% %m�͏���50000,100000,150000,200000,250000,300000,350000,400000,450000,500000,550000,600000

% 	foilname='hsfoil';
% 	datafile_a='_T1_Re';
% 	Re=Relist{12+m};
% 	datafile_b='_M0.00_N9.0.txt';
% 	dataname=strcat(foilname,datafile_a,Re,datafile_b);

% 	fpr=fopen(dataname(1,:));
% 	do
% 		buffer=fscanf(fpr,'%c',[1,1]);
% 	until(buffer(1,1)=='-')
% 	buffer=fgetl(fpr);
% 	ahdata{1,m}=(fscanf(fpr,'%f',[10,100]))';
% 	fclose(fpr);
% endfor

% cd(now_work)

% %-----�����������^�̓Ǎ�
% %-----�g�p���^��Ǎ�
% cd("NACA0010")

% %-----���^�̓Ǎ�
% fp=fopen('vsfoil.dat');
% vsfoil=fgetl(fp);
% vsfoil=(fscanf(fp,'%f',[2,250]))';
% fclose(fp);
% fsvs=0;


% %-----������͉�̓f�[�^��
% %-----���C�m���Y���X�g�ǂݍ���
% Relist={50000,'0.050';100000,'0.100'; 150000,'0.150'; 200000,'0.200'; 250000,'0.250'; 300000,'0.300';350000,'0.350'; 400000,'0.400';450000,'0.450'; 500000,'0.500';550000,'0.550';600000,'0.600';};;

% for m=1:12
% %m�͏���50000,100000,150000,200000,250000,300000,350000,400000,450000,500000,550000,600000

% 	foilname='vsfoil';
% 	datafile_a='_T1_Re';
% 	Re=Relist{12+m};
% 	datafile_b='_M0.00_N9.0.txt';
% 	dataname=strcat(foilname,datafile_a,Re,datafile_b);

% 	fpr=fopen(dataname(1,:));
% 	do
% 		buffer=fscanf(fpr,'%c',[1,1]);
% 	until(buffer(1,1)=='-')
% 	buffer=fgetl(fpr);
% 	avdata{1,m}=(fscanf(fpr,'%f',[10,100]))';
% 	fclose(fpr);
% endfor

% cd(now_work)


	% %�e�[���p�C�v�쐬
	% disp("�e�[���p�C�v�\���i�S�蓮���́j");
	% t_manddata(1,1)=input("�e�[����(mm)=");
	% t_manddata(1,2)=input("�e�[���a(mm)=");
	% t_manddata(1,3)=input("�ڍ�������(mm)=");
	% t_manddata(1,4)=t_manddata(1,2);%�}���h�����a�A�e�[���a�Ɠ���
	% t_manddata(1,5)=input("�e�[���}���h��������(mm)=");

	% %�ϑw�\������

	% %�����\������
	% t_plydata(1:4,2)=24;
	% %24t
	% t_plydata(1:4,3)=90;
	% %�ϑw�p�x
	% t_plydata(1:3:4,4)=90;
	% t_plydata(2:3,4)=45;
	% %�ϑw����
	% t_plydata(1:4,5)=0;
	% t_plydata(1:4,6)=t_manddata(1,1);

	% t_plydata(1,1)=t_manddata(cw,4)+C_T(lookup(C_T(:,1),t_plydata(1,2)),2)*2;
	% for j=2:4
	% 		t_plydata(j,1)=t_plydata(j-1,1)+C_T(lookup(C_T(:,1),t_plydata(j,2)),2)*2;
	% endfor
	% %�O�`�v�Z

	% %�蓮����
	% do
	% 	ts=input("�P�F�����ϑw�����@�Q�F����ϑw�����@�R�F�ϑw���� ---");

	% 	if ts==1
	% 		ra=input("�ϑw�p�x=");
	% 		t_plydata=ply_ex(t_plydata,1,ra,C_T,t_manddata(1,1))

	% 	elseif ts==2
	% 		rs=input("�w�ԍ�=");
	% 		t_plydata=ply_ins(t_plydata,rs,C_T,t_manddata(1,1));

	% 	elseif ts==3
	% 		rs=input("�w�ԍ�=");
	% 		t_plydata(rs,:)=[];
	% 		for j=2:rows(t_plydata)
	% 			t_plydata(j,1)=t_plydata(j-1,1)+C_T(lookup(C_T(:,1),t_plydata(j,2)),2)*2;
	% 		end
	% 	endif

	% 	%�ϑw�\����\��
	% 	ply_disp(t_plydata,000);
		
	% 	tsign=input("ok:1 return:0 ---")
			
	% until(tsign==0)

buf=input("deftail end")
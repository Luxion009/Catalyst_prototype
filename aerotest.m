%�܂��ȉ~���Ƃ��ċȂ����[�����g�����߂�
disp("��X�p����0.8~1.4�ŒT��")
span=input("��X�p��(m)=");
h_span=span/2;
d_s=h_span/100/2;
%��S����

%����ݗʂ����肷��
do
	disp("����ݗʂ͌Œ�ŒT��")
	max_tawami=input("�ő傽���(m)=");
	%�X�p��2m������0.154m���炢�ł����񂶂�Ȃ�����
	alpha=max_tawami/(h_span*h_span);
	a=0;
	d_count=1;
	f_v=0;

	%���f�̒[�_��y���W�����߂�i���_�Ȃ��j
	do
		line_e(d_count,1)=newton_met(a,alpha,d_s,f_v);
		a=line_e(d_count,1);
		f_v=a;
		d_count=d_count+1;
	
	until(d_count>100)

	%test=rows(line_e)

	d_count=1;

	%line_e(1,2)=alpha*((line_e(1,1))^2);

	%disp(line_e(1,:));
	
	%���f�̒[�_�̂����W�����߂�i���_�Ȃ��j
	do
		line_e(d_count,2)=alpha*((line_e(d_count,1))^2);
		d_count=d_count+1;
	
	until(d_count>100)

	%disp(line_e);
	
	%���f�̂��ꂼ��̒��������߂�i�P�O�O�����Ȃ̂łP�O�O�j
	line_e_d(1,1)=sqrt((line_e(1,1))^2+(line_e(1,2))^2);

	d_count=2;

	do
		line_e_d(d_count,1)=sqrt((line_e(d_count,1)-line_e(d_count-1,1))^2+(line_e(d_count,2)-line_e(d_count-1,2))^2);
		d_count=d_count+1;
	
	until(d_count>100)

	%disp(line_e_d);

	%test2=rows(line_e_d)
	
	printf("���K�ő傽���(m)=%f\n",line_e(100,2));
	sign=input("ok:1 return:0 ---");
	
until(sign==1)

%�R���g���[���|�C���g�i�㔼�p���l���������j�̍��W�����߂�
line_e_cp(1,1)=line_e(1,1)/2;
line_e_cp(1,2)=line_e(1,2)/2;
d_count=2;
do
	line_e_cp(d_count,1)=(line_e(d_count,1)-line_e(d_count-1,1))+line_e(d_count-1,1);
	line_e_cp(d_count,2)=(line_e(d_count,2)-line_e(d_count-1,2))+line_e(d_count-1,2);
	d_count=d_count+1;
	
until(d_count==101)

%�e���f�ɂ�����Ǐ��㔼�p�i���W�A���j�����߂�
line_e_ang(1,1)=atan(line_e(1,2)/line_e(1,1));
d_count=2;
do
	line_e_ang(d_count,1)=atan((line_e(d_count,2)-line_e(d_count-1,2))/(line_e(d_count,1)-line_e(d_count-1,1)));
	d_count=d_count+1;

until(d_count==101)

%�㔼�p���l�����Ȃ���������CP�i��y���W�̂݁j�����߂�
d_count=1;
do
	nang_cp(d_count,1)=d_s+(d_count-1)*2*d_s;
	d_count=d_count+1;

until(d_count==101)

%�����̊�Ȃ����[�����g��ȉ~����p���ĎZ�o����
do
	weight=input("�S���d��(kg)=");
	d_vel=input("�݌v�@��(m/s)=");
	a_density=1.164;
	nu=1.604*10^(-5);

	%�ȉ~�z�Ɖ��肵�ďz���z�����߂�
	root_gamma=4*weight*9.81/pi/a_density/d_vel/span;
	d_count=1;
	do
		d_gamma(d_count,1)=sqrt(root_gamma^2-(nang_cp(d_count,1)*root_gamma/span*2)^2);
		d_count=d_count+1;
		
	until(d_count==101)

	%�z���z���㔼�p���l����������CP�ɓK�����ėg�͂ƃ��[�����g�����Ƃ߂�
	%�܂��Ǐ��g�͂����߂�
	d_count=1;
	do
		d_l(d_count,1)=a_density*d_vel*d_gamma(d_count,1)*cos(line_e_ang(d_count,1))*line_e_d(d_count,1);
		d_count=d_count+1;

	until(d_count==101)
	%line_e_d��d_s�ł����Ȃ��H
	e_l=2*sum(d_l(:,1));
	%�����������Z�o���Ă��Ȃ��̂œ�{
	printf("�㔼�p���l���������̗g�́iN�j=%f\n",e_l);
	printf("�d�ʊ��Z��(kg)=%f\n",e_l/9.81);
	
	sign=input("ok:1 return:0 ---");

until(sign==1)

%�Ǐ����[�����g�����߂�i�ޗ͂ɂ����郂�[�����g�ł͂Ȃ��j
d_count=1;
do
	d_b(d_count,1)=a_density*d_vel*d_gamma(d_count,1)*(line_e_cp(d_count,1)*cos(line_e_ang(d_count,1))+line_e_cp(d_count,2)*sin(line_e_ang(d_count,1)))*line_e_d(d_count,1);
	d_count=d_count+1;

until(d_count==101)

bend_m=sum(d_b(:,1));
bend_m=bend_m*3*pi/2;
%�����܂őȉ~���Ɖ��肵���v�Z

%�U���R�͉�́i�E�C���O���b�g�͊܂܂Ȃ��j
div_count=input("��͕����� ---");
if(div_count!=0)
	for i=1:div_count
		s_span(i,1)=span*(0.8+0.6/(div_count-1)*(i-1));
	endfor

	eq_i=(div_count-1)/3+1;


	%�w��͈͂ŗU���R�͂̒l���Z�o
	for s_count=1:div_count
		printf("count %d/%d \n",s_count,div_count);
		%�ēx�X�p���Ƃ���݂��`���Ė{�v�Z�̏������s���i���܂ł̑ȉ~���̎��̃f�[�^�͋Ȃ����[�����g������������j
		span=s_span(s_count,1);
		h_span=span/2;
		d_s=h_span/100/2;
		%��S����
		
		%����ݗʂ͑ȉ~������Œ�
		alpha=max_tawami/(h_span*h_span);
		a=0;
		d_count=1;
		f_v=0;

		%���f�̒[�_��y���W�����߂�i���_�Ȃ��j
		do
			line_e(d_count,1)=newton_met(a,alpha,d_s,f_v);
			a=line_e(d_count,1);
			f_v=a;
			d_count=d_count+1;
		
		until(d_count>100)

		%test=rows(line_e)

		d_count=1;

		%line_e(1,2)=alpha*((line_e(1,1))^2);

		%disp(line_e(1,:));
		
		%���f�̒[�_�̂����W�����߂�i���_�Ȃ��j
		do
			line_e(d_count,2)=alpha*((line_e(d_count,1))^2);
			d_count=d_count+1;
		
		until(d_count>100)

		%disp(line_e);
		
		%���f�̂��ꂼ��̒��������߂�i�P�O�O�����Ȃ̂łP�O�O�j
		line_e_d(1,1)=sqrt((line_e(1,1))^2+(line_e(1,2))^2);

		d_count=2;

		do
			line_e_d(d_count,1)=sqrt((line_e(d_count,1)-line_e(d_count-1,1))^2+(line_e(d_count,2)-line_e(d_count-1,2))^2);
			d_count=d_count+1;
		
		until(d_count>100)

		%disp(line_e_d);

		%test2=rows(line_e_d)

		%�R���g���[���|�C���g�i�㔼�p���l���������j�̍��W�����߂�
		line_e_cp(1,1)=line_e(1,1)/2;
		line_e_cp(1,2)=line_e(1,2)/2;
		d_count=2;
		do
			line_e_cp(d_count,1)=(line_e(d_count,1)-line_e(d_count-1,1))/2+line_e(d_count-1,1);
			line_e_cp(d_count,2)=(line_e(d_count,2)-line_e(d_count-1,2))/2+line_e(d_count-1,2);
			d_count=d_count+1;
			
		until(d_count==101)
		%dlmwrite( "y.txt",line_e_cp(:,1),"delimiter","***","newline","\n");
		%dlmwrite( "z.txt",line_e_cp(:,2),"delimiter","***","newline","\n");

		%�e���f�ɂ�����Ǐ��㔼�p�i���W�A���j�����߂�
		line_e_ang(1,1)=atan(line_e(1,2)/line_e(1,1));
		d_count=2;
		do
			line_e_ang(d_count,1)=atan((line_e(d_count,2)-line_e(d_count-1,2))/(line_e(d_count,1)-line_e(d_count-1,1)));
			d_count=d_count+1;

		until(d_count==101)
		%dlmwrite( "ang.txt",line_e_ang(:,1),"delimiter","***","newline","\n");

		%�㔼�p���l�����Ȃ���������CP�i��y���W�̂݁j�����߂�
		d_count=1;
		do
			nang_cp(d_count,1)=d_s+(d_count-1)*2*d_s;
			d_count=d_count+1;

		until(d_count==101)

		%Ci�����߂�
		disp("calculating Ci.");
		d_count=1;
		do
			c_i(d_count,1)=2*cos(line_e_ang(d_count,1))*d_s;
			d_count=d_count+1;

		until(d_count==101)

		%Bi�����߂�
		disp("calculating Bi.");
		d_count=1;
		do
			b_i(d_count,1)=3*pi/2*(line_e_cp(d_count,1)*cos(line_e_ang(d_count,1))+line_e_cp(d_count,2)*sin(line_e_ang(d_count,1)))*d_s;
			d_count=d_count+1;
			
		until(d_count==101)

		%yd_ij�Ƃ������߂�
		disp("calculating yd_ij.");
		for i=1:100
			for j=1:100
				yd_ij(i,j)=(line_e_cp(i,1)-line_e_cp(j,1))*cos(line_e_ang(j,1))+(line_e_cp(i,2)-line_e_cp(j,2))*sin(line_e_ang(j,1));
				zd_ij(i,j)=-1*(line_e_cp(i,1)-line_e_cp(j,1))*sin(line_e_ang(j,1))+(line_e_cp(i,2)-line_e_cp(j,2))*cos(line_e_ang(j,1));
				ydd_ij(i,j)=(line_e_cp(i,1)+line_e_cp(j,1))*cos(line_e_ang(j,1))-(line_e_cp(i,2)-line_e_cp(j,2))*sin(line_e_ang(j,1));
				zdd_ij(i,j)=(line_e_cp(i,1)+line_e_cp(j,1))*sin(line_e_ang(j,1))-(line_e_cp(i,2)-line_e_cp(j,2))*cos(line_e_ang(j,1));
			endfor
			%printf("%d/100\n",i);
		endfor
		%dlmwrite( "ydij.txt",yd_ij(:,1),"delimiter","***","newline","\n");
		%dlmwrite( "zdij.txt",zd_ij(:,1),"delimiter","***","newline","\n");

		%R2_ij�Ƃ������߂�
		disp("calculating R2_ij.");
		for i=1:100
			for j=1:100
				R2_pij(i,j)=(yd_ij(i,j)-d_s)^2+(zd_ij(i,j))^2;
				R2_mij(i,j)=(yd_ij(i,j)+d_s)^2+(zd_ij(i,j))^2;
				Rd2_pij(i,j)=(ydd_ij(i,j)+d_s)^2+(zdd_ij(i,j))^2;
				Rd2_mij(i,j)=(ydd_ij(i,j)-d_s)^2+(zdd_ij(i,j))^2;
			endfor
			%printf("%d/100\n",i)
		endfor
		%dlmwrite( "r2pij.txt",R2_pij(:,1),"delimiter","***","newline","\n");

		%Q_ij�����߂�
		disp("calculating Q_ij.");
		for i=1:100
			for j=1:100
				Q_ij(i,j)=-1/2/pi*(((yd_ij(i,j)-line_e_d(j,1)./2)./R2_pij(i,j)-1*(yd_ij(i,j)+line_e_d(j,1)./2)./R2_mij(i,j))*cos(line_e_ang(i,1)-line_e_ang(j,1))+(zd_ij(i,j)./R2_pij(i,j)-zd_ij(i,j)./R2_mij(i,j))*sin(line_e_ang(i,1)-line_e_ang(j,1))+((ydd_ij(i,j)-line_e_d(j,1)./2)./Rd2_mij(i,j)-1*(ydd_ij(i,j)+line_e_d(j,1)./2)./Rd2_pij(i,j))*cos(line_e_ang(i,1)+line_e_ang(j,1))+(zdd_ij(i,j)./Rd2_mij(i,j)-zdd_ij(i,j)./Rd2_pij(i,j))*sin(line_e_ang(i,1)+line_e_ang(j,1)));
			endfor
		endfor
		%dlmwrite( "qij.txt",Q_ij(:,1),"delimiter","***","newline","\n");

		%A_ij�����߂�
		disp("calculating A_ij.");
		for i=1:100
			for j=1:100
				A_ij(i,j)=pi*Q_ij(i,j)*line_e_d(i,1)/2;
			endfor
		endfor

		%dlmwrite( "aij.txt",A_ij(:,1),"delimiter","***","newline","\n");

		%�œK���s������
		%�܂�A_ij����Ȃ�ꕔ�������
		for i=1:100
			for j=1:100
				Op_p1(i,j)=A_ij(i,j)+A_ij(j,i);
			endfor
		endfor

		%�T�C�h�̈ꕔ�������
		Op_p2=[-c_i -b_i];

		%���̈ꕔ�������
		Op_p3=[-(c_i)' 0 0;-(b_i)' 0 0];

		%�S�����킹��
		Op_mat=[Op_p1 Op_p2;Op_p3];

		%0000-1beta�̍s������
		for i=1:100
			Op_zero(i,1)=0;
		endfor
		%�݌v�g��
		del=weight*9.81;
		Op_r=[Op_zero;-del;-bend_m];
		a=[1.0 2.0;3.0 4.0];

		%�m�F�p�ɍs����o��
		%dlmwrite( "Opr.txt",Op_mat(:,1:10),"delimiter","***","newline","\n");

		%�{�v�Zgi�����߂�
		disp("calculating g_i.");
		g_i=Op_mat\Op_r;

		%gi�����ɏz���z�����߂�
		for i=1:100
			Op_gamma(i,1)=g_i(i,1)/2/a_density/d_vel;
		endfor

		%���������U�����x
		disp("calculating Vn_i.");
		for i=1:100
			for j=1:100
				Vn_i_d(j,i)=Q_ij(i,j)*Op_gamma(j,1)./2;
			endfor
			Vn_i(i,1)=sum(Vn_i_d(:,i));
			%printf("%d/100\n",i);
		endfor

		%�U���R��
		disp("calculating di_d.");
		for i=1:100
			di_d(i,1)=Op_gamma(i,1)*a_density*Vn_i(i,1)*line_e_d(i,1);
		endfor
		%�X�p�����ƂɗU���R�͂��i�[
		di(s_count,1)=2*sum(di_d(:,1));
		
	endfor

	%�w��͈̗͂U���R�͂̒l���o��
	figure(1);
		subplot(2,1,1);
		plot([0.8:0.6/(div_count-1):1.4],di(:,1)/di(eq_i,1));
		xlabel('span retio');
		ylabel('Di retio');
		xlim([0.8,1.4]);
		grid on;
		
		subplot(2,1,2);
		plot(s_span(:,1),di(:,1));
		xlabel('span');
		ylabel('Di');
		xlim([s_span(1,1),s_span(div_count,1)]);
		grid on;
		
		%-r100�͉𑜓x��100dpi�ɂ���D
	cd("design_data")
	cd(pj_name)
	print("Di_span.png",'-dpng','-r100')
	cd(now_work)
	
endif

%�O���t����X�p����I�ѐ����Ɍv�Z
%�ēx�X�p���Ƃ���݂��`���Ė{�v�Z�̏������s���i���܂ł̑ȉ~���̎��̃f�[�^�͋Ȃ����[�����g������������j
%���f�������͉�
do

	dw_divn=50;%������������
	disp("������{�v�Z�ł�");
	span=input("�{�X�p��(m)=");
	w_div_c=input("�嗃������=");
	we_div=w_div_c*dw_divn;%�嗃���ʕ�������

	for i=1:w_div_c
		disp(i);
		divw_span(i,1)=input("����������(m)=");
	endfor

	h_span=span/2;
	% d_s=h_span/we_div/2
	%�g���̂�߂���
	d_s(:,1)=divw_span(:,1)./dw_divn;

	%����ݗʂ����肷��
	do
		max_tawami=input("�ő傽���(m)=");
		%�X�p��2m������0.154m���炢�ł����񂶂�Ȃ�����
		alpha=max_tawami/(h_span*h_span);
		a=0;
		d_count=1;
		f_v=0;

		%���f�̒[�_��y���W�����߂�i���_�Ȃ��j
		for i=1:w_div_c
			for j=1:dw_divn
				line_e(dw_divn*(i-1)+j,1)=newton_met(a,alpha,d_s,f_v);
				a=line_e(d_count,1);
				f_v=a;
				d_count=d_count+1;
			
		until(d_count>we_div)

		%test=rows(line_e)

		d_count=1;

		%line_e(1,2)=alpha*((line_e(1,1))^2);

		%disp(line_e(1,:));
		
		%���f�̒[�_�̂����W�����߂�i���_�Ȃ��j
		do
			line_e(d_count,2)=alpha*((line_e(d_count,1))^2);
			d_count=d_count+1;
		
		until(d_count>we_div)

		%disp(line_e);
		
		%���f�̂��ꂼ��̒��������߂�i�P�O�O�����Ȃ̂łP�O�O�j
		line_e_d(1,1)=sqrt((line_e(1,1))^2+(line_e(1,2))^2);

		d_count=2;

		do
			line_e_d(d_count,1)=sqrt((line_e(d_count,1)-line_e(d_count-1,1))^2+(line_e(d_count,2)-line_e(d_count-1,2))^2);
			d_count=d_count+1;
		
		until(d_count>we_div)

		%disp(line_e_d);

		%test2=rows(line_e_d)
		
		printf("���K�ő傽���(m)=%f\n",line_e(we_div,2));
		sign=input("ok:1 return:0 ---");
		
	until(sign==1)

	%�R���g���[���|�C���g�i�㔼�p���l���������j�̍��W�����߂�
	line_e_cp(1,1)=line_e(1,1)/2;
	line_e_cp(1,2)=line_e(1,2)/2;
	d_count=2;
	do
		line_e_cp(d_count,1)=(line_e(d_count,1)-line_e(d_count-1,1))/2+line_e(d_count-1,1);
		line_e_cp(d_count,2)=(line_e(d_count,2)-line_e(d_count-1,2))/2+line_e(d_count-1,2);
		d_count=d_count+1;
		
	until(d_count==we_div+1)
	%dlmwrite( "y.txt",line_e_cp(:,1),"delimiter","***","newline","\n");
	%dlmwrite( "z.txt",line_e_cp(:,2),"delimiter","***","newline","\n");

	%�e���f�ɂ�����Ǐ��㔼�p�i���W�A���j�����߂�
	line_e_ang(1,1)=atan(line_e(1,2)/line_e(1,1));
	d_count=2;
	do
		line_e_ang(d_count,1)=atan((line_e(d_count,2)-line_e(d_count-1,2))/(line_e(d_count,1)-line_e(d_count-1,1)));
		d_count=d_count+1;

	until(d_count==we_div+1)
	%dlmwrite( "ang.txt",line_e_ang(:,1),"delimiter","***","newline","\n");
	
	%��������E�C���O���b�g�쐬��
	%�E�C���O���b�g���쐬����
	wl_ang=input("�E�C���O���b�g�̎��t���p�x�i����݂Ȃ����A����v���X�j=");
	wl_span=input("�E�C���O���b�g�̒����im�j=");
	
	if(wl_span==0)
	wld_c=0;
	endif
	
	if(wl_span>0)
		%�E�C���O���b�g�̕����������ɗ��f�̑傫�������߂�A�����ł͂Q�O����
		wld_c=20;
		wl_d=wl_span/wld_c;
		%�㔼�p���l�������Ƃ��̐�Ύ��t���p���Z�o����
		wl_ang_a=line_e_ang(we_div,1)+wl_ang/180*pi;
		%���f�P������̂��C�����W�ω������߂�
		wl_d_y=wl_d*cos(wl_ang_a);
		wl_d_z=wl_d*sin(wl_ang_a);
		%�P�F�Q�O�܂ł̗��f�̒[�_���W�����߂�
		for i=1:wld_c
			wl_e(i,1)=line_e(we_div,1)+wl_d_y*i;
			wl_e(i,2)=line_e(we_div,2)+wl_d_z*i;
		endfor
		%1:20�܂ł�CP���W�����߂�
		wl_e_cp(1,1)=line_e(we_div,1)+wl_d_y/2;
		wl_e_cp(1,2)=line_e(we_div,2)+wl_d_z/2;
		
		for i=2:wld_c
			wl_e_cp(i,1)=wl_e_cp(i-1,1)+wl_d_y;
			wl_e_cp(i,2)=wl_e_cp(i-1,2)+wl_d_z;
		endfor
		
		%���f�̃f�[�^���i�[���Ă���ϐ��ɃE�B���O���b�g�̃f�[�^���i�[����
		line_e((we_div+1):(wld_c+we_div),1)=wl_e(1:wld_c,1);
		line_e((we_div+1):(wld_c+we_div),2)=wl_e(1:wld_c,2);
		line_e_cp((we_div+1):(wld_c+we_div),1)=wl_e_cp(1:wld_c,1);
		line_e_cp((we_div+1):(wld_c+we_div),2)=wl_e_cp(1:wld_c,2);
		line_e_ang((we_div+1):(wld_c+we_div),1)=wl_ang_a;
		line_e_d((we_div+1):(wld_c+we_div),1)=wl_d;
		
	endif
	
	%�㔼�p���l�����Ȃ���������CP�i��y���W�̂݁j�����߂�
	d_count=1;
	do
		nang_cp(d_count,1)=line_e_d(d_count,1)/2+(d_count-1)*line_e_d(d_count,1);
		d_count=d_count+1;

	until(d_count==we_div+1)
	
	if(wl_span>0)
		nang_cp(we_div+1,1)=nang_cp(we_div,1)+line_e_d(we_div,1)/2+wl_d/2;
		for i=(we_div+2):(we_div+wld_c)
			nang_cp(i,1)=nang_cp(i-1,1)+wl_d;
		endfor

	endif
	
	dlmwrite( "log/nang_cp.txt",nang_cp(:,1),"delimiter","***","newline","\n");

	%Ci�����߂�
	disp("calculating Ci.");
	d_count=1;
	do
		c_i(d_count,1)=2*cos(line_e_ang(d_count,1))*line_e_d(d_count,1)/2;
		d_count=d_count+1;

	until(d_count==wld_c+we_div+1)

	%Bi�����߂�
	disp("calculating Bi.");
	d_count=1;
	do
		b_i(d_count,1)=3*pi/2*(line_e_cp(d_count,1)*cos(line_e_ang(d_count,1))+line_e_cp(d_count,2)*sin(line_e_ang(d_count,1)))*line_e_d(d_count,1)/2;
		d_count=d_count+1;
		
	until(d_count==wld_c+we_div+1)

	%yd_ij�Ƃ������߂�
	disp("calculating yd_ij.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			yd_ij(i,j)=(line_e_cp(i,1)-line_e_cp(j,1))*cos(line_e_ang(j,1))+(line_e_cp(i,2)-line_e_cp(j,2))*sin(line_e_ang(j,1));
			zd_ij(i,j)=-1*(line_e_cp(i,1)-line_e_cp(j,1))*sin(line_e_ang(j,1))+(line_e_cp(i,2)-line_e_cp(j,2))*cos(line_e_ang(j,1));
			ydd_ij(i,j)=(line_e_cp(i,1)+line_e_cp(j,1))*cos(line_e_ang(j,1))-(line_e_cp(i,2)-line_e_cp(j,2))*sin(line_e_ang(j,1));
			zdd_ij(i,j)=(line_e_cp(i,1)+line_e_cp(j,1))*sin(line_e_ang(j,1))-(line_e_cp(i,2)-line_e_cp(j,2))*cos(line_e_ang(j,1));
		endfor
		%printf("%d/100\n",i);
	endfor
	%dlmwrite( "ydij.txt",yd_ij(:,1),"delimiter","***","newline","\n");
	%dlmwrite( "zdij.txt",zd_ij(:,1),"delimiter","***","newline","\n");

	%R2_ij�Ƃ������߂�
	disp("calculating R2_ij.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			R2_pij(i,j)=(yd_ij(i,j)-line_e_d(j,1)/2)^2+(zd_ij(i,j))^2;
			R2_mij(i,j)=(yd_ij(i,j)+line_e_d(j,1)/2)^2+(zd_ij(i,j))^2;
			Rd2_pij(i,j)=(ydd_ij(i,j)+line_e_d(j,1)/2)^2+(zdd_ij(i,j))^2;
			Rd2_mij(i,j)=(ydd_ij(i,j)-line_e_d(j,1)/2)^2+(zdd_ij(i,j))^2;
		endfor
		%printf("%d/100\n",i)
	endfor
	%dlmwrite( "r2pij.txt",R2_pij(:,1),"delimiter","***","newline","\n");

	%Q_ij�����߂�
	disp("calculating Q_ij.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Q_ij(i,j)=-1/2/pi*(((yd_ij(i,j)-line_e_d(j,1)./2)./R2_pij(i,j)-1*(yd_ij(i,j)+line_e_d(j,1)./2)./R2_mij(i,j))*cos(line_e_ang(i,1)-line_e_ang(j,1))+(zd_ij(i,j)./R2_pij(i,j)-zd_ij(i,j)./R2_mij(i,j))*sin(line_e_ang(i,1)-line_e_ang(j,1))+((ydd_ij(i,j)-line_e_d(j,1)./2)./Rd2_mij(i,j)-1*(ydd_ij(i,j)+line_e_d(j,1)./2)./Rd2_pij(i,j))*cos(line_e_ang(i,1)+line_e_ang(j,1))+(zdd_ij(i,j)./Rd2_mij(i,j)-zdd_ij(i,j)./Rd2_pij(i,j))*sin(line_e_ang(i,1)+line_e_ang(j,1)));
		endfor
	endfor
	%dlmwrite( "qij.txt",Q_ij(:,1),"delimiter","***","newline","\n");

	%A_ij�����߂�
	disp("calculating A_ij.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			A_ij(i,j)=pi*Q_ij(i,j)*line_e_d(i,1)/2;
		endfor
	endfor

	%dlmwrite( "aij.txt",A_ij(:,1),"delimiter","***","newline","\n");

	%�œK���s������
	%�܂�A_ij����Ȃ�ꕔ�������
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Op_p1(i,j)=A_ij(i,j)+A_ij(j,i);
		endfor
	endfor

	%�T�C�h�̈ꕔ�������
	Op_p2=[-c_i -b_i];

	%���̈ꕔ�������
	Op_p3=[-(c_i)' 0 0;-(b_i)' 0 0];

	%�S�����킹��
	Op_mat=[Op_p1 Op_p2;Op_p3];

	%0000-1beta�̍s������
	for i=1:(wld_c+we_div)
		Op_zero(i,1)=0;
	endfor
	%�݌v�g��
	del=weight*9.81;
	Op_r=[Op_zero;-del;-bend_m];
	a=[1.0 2.0;3.0 4.0];

	%�m�F�p�ɍs����o��
	%dlmwrite( "Opr.txt",Op_mat(:,1:10),"delimiter","***","newline","\n");

	%�{�v�Zgi�����߂�
	disp("calculating g_i.");
	g_i=Op_mat\Op_r;

	%gi�����ɏz���z�����߂�
	for i=1:(wld_c+we_div)
		Op_gamma(i,1)=g_i(i,1)/2/a_density/d_vel;
	endfor
	dlmwrite( "log/Opgamma.txt",Op_gamma(:,1),"delimiter","***","newline","\n");

	%���������U�����x
	disp("calculating Vn_i.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Vn_i_d(j,i)=Q_ij(i,j)*Op_gamma(j,1)./2;
		endfor
		Vn_i(i,1)=sum(Vn_i_d(:,i));
		%printf("%d/100\n",i);
	endfor

	%�U���p�x
	disp("calculating i_ang.");
	for i=1:(wld_c+we_div)
		i_ang(i,1)=atan(Vn_i(i,1)/d_vel)*180/pi;
	endfor
	
	%�������x�v�Z
	disp("calculating comp_vel");
	for i=1:(wld_c+we_div)
		comp_vel(i,1)=sqrt(d_vel^2+Vn_i(i,1)^2);
	endfor

	%�U���R��
	disp("calculating di_d.");
	for i=1:(wld_c+we_div)
		di_d(i,1)=Op_gamma(i,1)*Vn_i(i,1)*line_e_d(i,1);
	endfor
	di=2*a_density*sum(di_d(:,1));

	%�Ǐ��g�͂����߂�(�ޗ͗p,�㔼�Ȃ�)
	disp("calculating d_l.");
	for i=1:(wld_c+we_div)
		d_l(i,1)=Op_gamma(i,1)*d_vel*a_density*line_e_d(i,1);
	endfor

	%����f�͂����߂�(�E�C���O���b�g�Ȃ�)
	disp("calculating Q.");
	Q(we_div,1)=d_l(we_div,1);
	for i=(we_div-1):-1:1
		Q(i,1)=Q(i+1,1)+d_l(i,1);
	endfor
	%dlmwrite( "Q.txt",Q(:,1),"delimiter","***","newline","\n");

	if(wl_span>0)
		%�E�C���O���b�g���̂���f�͂����߂�
		Q(wld_c+we_div,1)=d_l(wld_c+we_div,1);
		for i=(we_div-1+wld_c):-1:we_div+1
			Q(i,1)=Q(i+1,1)+d_l(i,1);
		endfor
	endif
	
	%�Ǐ����[�����g�����߂�(mm�ɕϊ�,�E�C���O���b�g�Ȃ�)
	disp("calculating d_m.");
	d_m(we_div,1)=Q(we_div,1)*line_e_d(we_div,1)/2*1000;
	for i=(we_div-1):-1:1
		d_m(i,1)=(Q(i,1)+Q(i+1))*line_e_d(i,1)/2*1000;
	endfor
	%dlmwrite( "dM.txt",d_m(:,1),"delimiter","***","newline","\n");
	
	if(wl_span>0)
		%�E�C���O���b�g���̋Ǐ����[�����g�����߂�
		d_m(wld_c+we_div,1)=Q(wld_c+we_div,1)*line_e_d(wld_c+we_div,1)/2*1000;
		for i=(we_div-1+wld_c):-1:(we_div+1)
			d_m(i,1)=(Q(i,1)+Q(i+1))*line_e_d(i,1)/2*1000;
		endfor
		
		%�E�C���O���b�g�̃��[�����g�����߂ė��[�Ƀ}�C�i�X�������[�����g�׏d�Ƃ��Čv�Z����
		M(wld_c+we_div,1)=d_m(wld_c+we_div,1);
		for i=(we_div-1+wld_c):-1:we_div+1
			M(i,1)=M(i+1,1)+d_m(i,1);
		endfor
		%M(we_div+1,1)�����[�ɂ�����ƍl�����郂�[�����g�׏d
	endif

	if(wl_span==0)	
		M(we_div+1,1)=0;
	endif
		
	%���[�����g�����߂�(N*mm)
	disp("calculating M.");
	M(we_div,1)=d_m(we_div,1)-M(we_div+1,1);
	for i=(we_div-1):-1:1
		M(i,1)=M(i+1,1)+d_m(i,1);
	endfor
	%M(1,1)���������Ȃ����[�����g
	%dlmwrite( "M.txt",M(:,1),"delimiter","***","newline","\n");

	%�g�͌v�Z
	disp("calculating L.");
	for i=1:(wld_c+we_div)
		dd_l(i,1)=Op_gamma(i,1)*cos(line_e_ang(i,1))*line_e_d(i,1)/2;
	endfor
	re_def_l=4*a_density*d_vel*sum(dd_l(:,1));

	%�ȉ~�z�Ɖ��肵�ďz���z�����߂�(�z���z���m�F���邽��)
	root_gamma=4*weight*9.81/pi/a_density/d_vel/(span+wl_span*2);
	d_count=1;
	do
		d_gamma(d_count,1)=sqrt(root_gamma^2-(nang_cp(d_count,1)*root_gamma/(span+wl_span*2)*2)^2);
		d_count=d_count+1;
		
	until(d_count==wld_c+we_div+1)


	%%%�����̂��ߏz���z�𑽊p�`�ߎ�����
	for i=1:5
		gamma_d(i,1)=(Op_gamma(i*40+1,1)-Op_gamma((i-1)*40+1,1))/40;
		
		for j=1:40
			Hex_gamma((i-1)*40+j,1)=Op_gamma((i-1)*40+1,1)+gamma_d(i,1)*(j-1);
		endfor
	endfor

	Hex_gamma(201:220,1)=Op_gamma(201:220,1);


	%���������U�����x
	disp("calculating Hex_Vn_i.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Hex_Vn_i_d(j,i)=Q_ij(i,j)*Hex_gamma(j,1)./2;
		endfor
		Hex_Vn_i(i,1)=sum(Hex_Vn_i_d(:,i));
		%printf("%d/100\n",i);
	endfor

	%�U���p�x
	disp("calculating Hex_i_ang.");
	for i=1:(wld_c+we_div)
		Hex_i_ang(i,1)=atan(Hex_Vn_i(i,1)/d_vel)*180/pi;
	endfor
	
	%�������x�v�Z
	disp("calculating Hex_comp_vel");
	for i=1:(wld_c+we_div)
		Hex_comp_vel(i,1)=sqrt(d_vel^2+Hex_Vn_i(i,1)^2);
	endfor

	%�U���R��
	disp("calculating Hex_di_d.");
	for i=1:(wld_c+we_div)
		Hex_di_d(i,1)=Hex_gamma(i,1)*Hex_Vn_i(i,1)*line_e_d(i,1);
	endfor
	Hex_di=2*a_density*sum(Hex_di_d(:,1));



	%�v�Z�l���o��
	disp("�v�Z����");
	printf("�@��(m/s)=%f\n",d_vel);
	printf("�X�p��(m)=%f\n",span);
	printf("�U���R��(N)=%f\n",di);
	printf("���p�`�U���R��(N)=%f\n",Hex_di);

	%�o��
	figure(2);

		subplot(2,1,1);
			plot(nang_cp(:,1),d_gamma(:,1),nang_cp(:,1),Op_gamma(:,1),"-",nang_cp(:,1),Hex_gamma(:,1),"-");
			xlabel('y[m]');
			ylabel('gamma');
			xlim([0,span/2+wl_span+1]);
			grid on;
		
		subplot(2,1,2);
			plot(nang_cp(:,1),di_d(:,1),"-");
			xlabel('y[m]');
			ylabel('Di');
			xlim([0,span/2+wl_span+1]);
			grid on;
	cd("design_data")
	cd(pj_name)
	print("Di_gamma.png","-dpng","-r100")
	cd(now_work)

	figure(5);
		subplot(2,1,1);
			plot(nang_cp(:,1),Vn_i(:,1),"-");
			xlabel('y[m]');
			ylabel('Vn_i[m/s]');
			xlim([0,span/2+wl_span+1]);
			grid on;
			
		subplot(2,1,2);
			plot(nang_cp(:,1),i_ang(:,1),"-");
			xlabel('y[m]');
			ylabel('i_angle');
			xlim([0,span/2+wl_span+1]);
			grid on;
			
		%-r100�͉𑜓x��100dpi�ɂ���D
	cd("design_data")
	cd(pj_name)
	print("Vn_iang.png","-dpng","-r100")
	cd(now_work)


	figure(3);
		plot(line_e(:,1),line_e(:,2),"-;wing curve;");
		xlabel('y[m]');
		ylabel('z[m]');
		ylim([-1,5]);
		axis("equal");
		grid on;

	cd("design_data")
	cd(pj_name)
	print("wingcurve.png",'-dpng','-r100')
	cd(now_work)
	
	figure(4);
		plot(1:we_div,M(1:we_div,1),"-;M;");
		xlabel('CP number');
		ylabel('M[N*mm]');
		grid on;

	cd("design_data")
	cd(pj_name)
	print("CP_M.png",'-dpng','-r100')
	cd(now_work)

	sign=input("ok:1 return:0 ---");
	
until(sign==1)
%���݌v
%�������͂���͂�����brazier�����Ƃ�ply�������肷��

%-----�J�[�{������(2���24t,3���40t)kgfmm
	E=[0,24,40;0,13000,22000;45,1900,1900;90,900,800];
%-----�v���v���O���x(g/mm^3)
	C_D=[24,0.001496;40,0.001559];
%-----�`����(mm)
	C_T=[24,0.125;40,0.111];
%-----�񎟕��ޖ��x(kg/m2)
	Snd_rho=0.4;
%-----�|�A�\����
	Pois=0.3;
%%skyscraper���玝���Ă���

%���̉�͎��̕����_
M_divp=[25 50 75 100];

%�����ϑw�p�x
Rply_ang=[45 40 35 30 25];

%�����ϑw�̌��E��
Cply_lim=5;

%�v�����N��
p_thn=2.5;

%���E���v�����N�ԗ]�T
lim_gap=2;

do
	clear ply_data;

	% LBS=input("��������(Mpa)=");
	XLS=input("X�������E�Ȃ�����(Mpa)=");
	E_xy=22000*9.8;
	Poison_r=0.3;
	%D_t_stand=1/LBS*2*sqrt(2)/9*E_xy/(1-Poison_r^2);
	%d/t�𓱏o

	%�S���ɂ킽��90/45/45/90��ϑw
	for cw=1:w_div_c
		
		ply_data{cw}(1:4,2)=24;
		%24t
		ply_data{cw}(1:4,3)=90;
		%�ϑw�p�x
		ply_data{cw}(1:3:4,4)=90;
		ply_data{cw}(2:3,4)=45;
		%�ϑw����
		ply_data{cw}(1:4,5)=0;
		ply_data{cw}(1:4,6)=divw_span(cw,1)*1000;
		
		ply_data{cw}(1,1)=mand_data(cw,4)+C_T(lookup(C_T(:,1),ply_data{cw}(1,2)),2)*2;
		for j=2:4
				ply_data{cw}(j,1)=ply_data{cw}(j-1,1)+C_T(lookup(C_T(:,1),ply_data{cw}(j,2)),2)*2;
		endfor
		%�O�`�v�Z
		
	endfor
	disp("generated start ply");

	%X�����ő剞�͊�ň���ϑw������

	for xpc=1:w_div_c %X��v���C�̃J�E���g
		xplyc(xpc,1:length(M_divp))=0;

	endfor

	for ci=1:w_div_c

		%�ꗃ���Ƃɐϑw������
		%�����l�������č\������
		
		printf("generating x ply [%d wing] \n",ci);		
		x_bend_disp
		
		% figure(18);
		% 	plot(S_mat(:,1),XLS,"-",S_mat(:,1),S_max_X(:,1),"-");
		% 	xlabel('y[m]');
		% 	ylabel('XS[Mpa]');
		% 	grid on;
		% print("output/defspar.png","-dpng","-r100")

		% secend=input("PRESS ENTER---");

		for cv=1:length(M_divp)%�Z�N�V�����ԍ�1-4
			
			printf("generating x ply [%d section] \n",cv);
			
			if(S_max_X(dp_num*(ci-1)+M_divp(1,cv)-1,1)>XLS)%X�����ő剞�͂��傫����΃v���C�𑝂�

				disp("plying");

				do
					if(cv==1)%�œ��Z�N�V������������
						disp("ply_ex");
						ply_data{ci}=ply_ex(ply_data{ci},1,90,C_T,mand_data(ci,1)./4.*1000);
						xplyc(ci,1)++;
					else%����ȊO�̃Z�N�V������������
						disp("ply_leng");
						cp=xplyc(ci,cv)+4;
						ply_data{ci}=ply_leng(ply_data{ci},cp,mand_data(ci,1)./4.*cv.*1000);
						xplyc(ci,cv)++;
					endif
					
					x_bend_disp

					smax=S_max_X(dp_num*(ci-1)+M_divp(1,cv)-1,1)
					% si=S_I(dp_num*(ci-1)+M_divp(1,cv),1)
					% xx=Xx(dp_num*(ci-1)+M_divp(1,cv),1)
					% rr=ply_data{ci}(rows(ply_data{ci}),1)

					% secend=input("PRESS ENTER---");

					
				until(S_max_X(dp_num*(ci-1)+M_divp(1,cv)-1,1)<XLS)
			
			endif
			
		endfor

	endfor


	figure(18);
	plot(S_mat(:,1),XLS,"-",S_mat(:,1),S_max_X(:,1),"-");
	xlabel('y[m]');
	ylabel('XS[Mpa]');
	grid on;
	% % print("output/defspar.png","-dpng","-r100")

	secend=input("PRESS ENTER---");



	% %�S���ɂ킽��0�x�w��brazier���ply
	% for cw=1:w_div_c
	
	% 	braz_ply=round(((mand_data(cw,4)/D_t_stand)/2-0.125*4)/0.111)
	% 	if braz_ply<=0
	% 		braz_ply=0
	% 	end
		
	% 	%brazier�̎�����o���l������ply_ex�Őϑw�𑝂�
	% 	if(braz_ply>0)
	% 		ply_data{cw}=ply_ex(ply_data{cw},braz_ply,90,C_T,mand_data(cw,1)*1000);
	% 	endif
	% 	ply_data{cw}
	% endfor
	
	%�ϑw�\����\��
	for i=1:w_div_c
		ply_disp(ply_data{i},i);
	endfor

	sign=input("1:ok 0:return ---");
	
until(sign==1)

%����݋Ȑ��ɍ��킹�ď㉺�����ϑw�J�n

%�����ϑw�J�E���^�[������
%zply_c(x,1)�͕����ϑw�������邽�߂̂���(x,2:4)�͉��w���₵�����������
for x=1:w_div_c
	zply_c(x,1:length(M_divp))=0;
endfor

for ci=1:w_div_c

	%�ꗃ���Ƃɐϑw������
	%�����ϑw��5�w�݂̂���ȏ�͈���ϑw�őΉ�
	%�����l�������č\������
	
	printf("generating ply [%d wing] \n",ci)
	
	bend_disp
	
	% figure(18);
	% 	plot(S_mat(:,1),S_mat_swang(:,1),"-",S_mat(:,1),ID_swang(:,1),"-");
	% 	xlabel('y[m]');
	% 	ylabel('z[m]');
	% 	grid on;
	% 	axis equal;
	% print("output/defspar.png","-dpng","-r100")
	
	for cv=1:length(M_divp)%�Z�N�V�����ԍ�1-4
		
		printf("generating ply [%d section] \n",cv);
		
		if(S_mat_swang(dp_num*(ci-1)+M_divp(1,cv)-1,1)>ID_swang(dp_num*(ci-1)+M_divp(1,cv)-1,1))%���z���݂��傫����΃v���C�𑝂�
			do

				psign=0;

				if(cv==1)%�œ��Z�N�V������������
					if(zply_c(ci,1)<Cply_lim)%�W���ϑw���E���ȉ���������
						disp("ply_ex");
						r_ang=Rply_ang(1,zply_c(ci,1)+1);%�֐��̓���q�G���[�h�~��
						ply_data{ci}=ply_ex(ply_data{ci},1,r_ang,C_T,mand_data(ci,1)./4.*1000);
						zply_c(ci,1)++;
					else%���E���𒴂��Ă�����
						disp("ply_ins");
						ply_ins_end=mand_data(ci,1)./4.*1000;
						insp=xplyc(ci,cv)+4;
						ply_data{ci}=ply_ins(ply_data{ci},insp,C_T,ply_ins_end);
					endif
				else%����ȊO�̃Z�N�V������������

					cp=zply_c(ci,cv)+xplyc(ci,cv)+4;%�L�΂��ׂ��w�ԍ�

					if(cp>rows(ply_data{ci})-1&&cpdn<lpdn)%�w�ԍ����ŊO�w�ɗ��Ă��܂����Ƃ��͍œ��Z�N�V�����̐ϑw��ύX����

						if(zply_c(ci,1)<Cply_lim)%�W���ϑw���E���ȉ���������
							disp("re_ply_ex");
							r_ang=Rply_ang(1,zply_c(ci,1)+1);%�֐��̓���q�G���[�h�~��
							ply_data{ci}=ply_ex(ply_data{ci},1,r_ang,C_T,mand_data(ci,1)./4.*cv.*1000);
							zply_c(ci,1)++;
						else%���E���𒴂��Ă�����
							disp("re_ply_ins");
							ply_ins_end=mand_data(ci,1)./4.*cv.*1000
							insp=xplyc(ci,cv)+4;
							ply_data{ci}=ply_ins(ply_data{ci},insp,C_T,ply_ins_end);
						
						endif
					elseif cp<=(rows(ply_data{ci})-1)
						
						disp("ply_leng");	
						ply_data{ci}=ply_leng(ply_data{ci},cp,mand_data(ci,1)/4*cv*1000);
						zply_c(ci,cv)++;

					else
						disp("ply_end");
						psign=1;
						
					endif

				endif
				
				bend_disp

				swang=S_mat_swang(dp_num*(ci-1)+M_divp(1,cv),1)
				id=ID_swang(dp_num*(ci-1)+M_divp(1,cv),1)

				if ci==1
					lpdn=200;%�K���ɑ傫����
				else
					lpdn=rows(ply_data{ci-1});
				endif
				
				cpdn=rows(ply_data{ci});

				if cv==1&&cpdn>=lpdn%���[�ɍs���ɂ�ϑw��������悤�ɂ��낢��
					psign=1;
				endif
	
			until(S_mat_swang(dp_num*(ci-1)+M_divp(1,cv)-1,1)<ID_swang(dp_num*(ci-1)+M_divp(1,cv)-1,1)||psign==1)
		
		endif
		
	endfor
	
endfor

figure(19);
plot(S_mat(:,1),S_mat_swang(:,1),"-;swang;",S_mat(:,1),ID_swang(:,1),"-;id;");
xlabel('y[m]');
ylabel('z[m]');
grid on;
axis equal;
% print("output/defspar.png","-dpng","-r100")

figure(20);
plot(S_mat(:,1),S_max_Z(:,1),"-");
xlabel('y[m]');
ylabel('s_max_z[Mpa]');
grid on;
% print("output/defspar.png","-dpng","-r100")

%�ϑw�\����\��
for i=1:w_div_c
	ply_disp(ply_data{i},i);
endfor

%�蓮�C��
disp("�蓮�C��");
ss=input("1:yes 0:no ---");
if(ss==1)
	do
	ss2=input("�P�F�C���@�Q�F�����ϑw�����@�R�F����ϑw�����@�S�F�ϑw���� ---");

	if ss2==1
		rw=input("�C�����ԍ�=");
		rs=input("�w�ԍ�=");
		rl=input("�Ē�`�͈� 1:1-1 2:1-2 3:1-3 4:1-4 ---");
		ply_data{rw}=ply_leng(ply_data{rw},rs,mand_data(rw,1)/4*rl*1000);

	elseif ss2==2
		rw=input("�������ԍ�=");
		ra=input("�ϑw�p�x=");
		rl=input("��`�͈� 1:1-1 2:1-2 3:1-3 4:1-4 ---");

		ply_data{rw}=ply_ex(ply_data{rw},1,ra,C_T,mand_data(rw,1)./4*rl*1000)

	elseif ss2==3
		rw=input("�������ԍ�=");
		rs=input("�w�ԍ�=");
		rl=input("��`�͈� 1:1-1 2:1-2 3:1-3 4:1-4 ---");

		ply_data{rw}=ply_ins(ply_data{rw},rs,C_T,mand_data(rw,1)./4*rl*1000);

	elseif ss2==4
		rw=input("�������ԍ�=");
		rs=input("�w�ԍ�=");
		ply_data{rw}(rs,:)=[];
		for j=2:rows(ply_data{rw})
			ply_data{rw}(j,1)=ply_data{rw}(j-1,1)+C_T(lookup(C_T(:,1),ply_data{rw}(j,2)),2)*2;
		end

	endif

	%�ϑw�\����\��
	for i=1:w_div_c
		ply_disp(ply_data{i},i);
	endfor

	bend_disp
	x_bend_disp

	figure(19);
	plot(S_mat(:,1),S_mat_swang(:,1),"-;swang;",S_mat(:,1),ID_swang(:,1),"-;id;");
	xlabel('y[m]');
	ylabel('z[m]');
	grid on;
	axis equal;
	% print("output/defspar.png","-dpng","-r100")

	figure(43)
	plot(S_mat(:,1),S_mat_Dswang(:,1),"-;swang;");
	xlabel('y[m]');
	ylabel('z[m]');
	grid on;

	figure(20);
	plot(S_mat(:,1),S_max_Z(:,1),"-");
	xlabel('y[m]');
	ylabel('s_max_z[Mpa]');
	grid on;
	% print("output/defspar.png","-dpng","-r100")

	figure(21)
	plot(S_mat(:,1),D_mat(:,1),"-");
	xlabel("y[m]");
	ylabel("line density");
	grid on;

	%�����ʂ邩�ŏI�m�F

	for i=1:w_div_c
		ps_gap(i,1)=(real_thn_xcp(i,1)*1000-ply_data{i}(rows(ply_data{i}),1))/2-p_thn;
		if(ps_gap(i,1)>=lim_gap)
			printf("[%d��]�@��-�v�����N�ԗ]�T : %f \n",i,ps_gap(i,1));
		else
			printf("[%d��]�@���E�]�T�ȉ��@����`�s�\ : %f \n",i,ps_gap(i,1));
		end
	end

	W_Spar=W_Spar
	W_wing=W_wing

	ssign=input("next:1 end:0 ---")
		
	until(ssign==0)

endif	



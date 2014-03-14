%桁設計
%座屈応力を入力させてbrazierをもとにply数を決定する

%-----カーボン物性(2列目24t,3列目40t)kgfmm
	E=[0,24,40;0,13000,22000;45,1900,1900;90,900,800];
%-----プリプレグ密度(g/mm^3)
	C_D=[24,0.001496;40,0.001559];
%-----形成厚(mm)
	C_T=[24,0.125;40,0.111];
%-----二次部材密度(kg/m2)
	Snd_rho=0.4;
%-----ポアソン比
	Pois=0.3;
%%skyscraperから持ってきた

%桁の解析時の分割点
M_divp=[25 50 75 100];

%部分積層角度
Rply_ang=[45 40 35 30 25];

%部分積層の限界数
Cply_lim=5;

%プランク厚
p_thn=2.5;

%限界桁プランク間余裕
lim_gap=2;

do
	clear ply_data;

	% LBS=input("座屈応力(Mpa)=");
	XLS=input("X方向限界曲げ応力(Mpa)=");
	E_xy=22000*9.8;
	Poison_r=0.3;
	%D_t_stand=1/LBS*2*sqrt(2)/9*E_xy/(1-Poison_r^2);
	%d/tを導出

	%全翼にわたり90/45/45/90を積層
	for cw=1:w_div_c
		
		ply_data{cw}(1:4,2)=24;
		%24t
		ply_data{cw}(1:4,3)=90;
		%積層角度
		ply_data{cw}(1:3:4,4)=90;
		ply_data{cw}(2:3,4)=45;
		%積層方向
		ply_data{cw}(1:4,5)=0;
		ply_data{cw}(1:4,6)=divw_span(cw,1)*1000;
		
		ply_data{cw}(1,1)=mand_data(cw,4)+C_T(lookup(C_T(:,1),ply_data{cw}(1,2)),2)*2;
		for j=2:4
				ply_data{cw}(j,1)=ply_data{cw}(j-1,1)+C_T(lookup(C_T(:,1),ply_data{cw}(j,2)),2)*2;
		endfor
		%外形計算
		
	endfor
	disp("generated start ply");

	%X方向最大応力基準で一周積層を決定

	for xpc=1:w_div_c %X基準プライのカウント
		xplyc(xpc,1:length(M_divp))=0;

	endfor

	for ci=1:w_div_c

		%一翼ごとに積層を決定
		%桁を四分割して構成決定
		
		printf("generating x ply [%d wing] \n",ci);		
		x_bend_disp
		
		% figure(18);
		% 	plot(S_mat(:,1),XLS,"-",S_mat(:,1),S_max_X(:,1),"-");
		% 	xlabel('y[m]');
		% 	ylabel('XS[Mpa]');
		% 	grid on;
		% print("output/defspar.png","-dpng","-r100")

		% secend=input("PRESS ENTER---");

		for cv=1:length(M_divp)%セクション番号1-4
			
			printf("generating x ply [%d section] \n",cv);
			
			if(S_max_X(dp_num*(ci-1)+M_divp(1,cv)-1,1)>XLS)%X方向最大応力より大きければプライを増す

				disp("plying");

				do
					if(cv==1)%最内セクションだった時
						disp("ply_ex");
						ply_data{ci}=ply_ex(ply_data{ci},1,90,C_T,mand_data(ci,1)./4.*1000);
						xplyc(ci,1)++;
					else%それ以外のセクションだった時
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



	% %全翼にわたり0度層をbrazier基準でply
	% for cw=1:w_div_c
	
	% 	braz_ply=round(((mand_data(cw,4)/D_t_stand)/2-0.125*4)/0.111)
	% 	if braz_ply<=0
	% 		braz_ply=0
	% 	end
		
	% 	%brazierの式から出た値を元にply_exで積層を増す
	% 	if(braz_ply>0)
	% 		ply_data{cw}=ply_ex(ply_data{cw},braz_ply,90,C_T,mand_data(cw,1)*1000);
	% 	endif
	% 	ply_data{cw}
	% endfor
	
	%積層構成を表示
	for i=1:w_div_c
		ply_disp(ply_data{i},i);
	endfor

	sign=input("1:ok 0:return ---");
	
until(sign==1)

%たわみ曲線に合わせて上下部分積層開始

%部分積層カウンター初期化
%zply_c(x,1)は部分積層数を入れるためのもの(x,2:4)は何層増やしたか見るもの
for x=1:w_div_c
	zply_c(x,1:length(M_divp))=0;
endfor

for ci=1:w_div_c

	%一翼ごとに積層を決定
	%部分積層は5層のみそれ以上は一周積層で対応
	%桁を四分割して構成決定
	
	printf("generating ply [%d wing] \n",ci)
	
	bend_disp
	
	% figure(18);
	% 	plot(S_mat(:,1),S_mat_swang(:,1),"-",S_mat(:,1),ID_swang(:,1),"-");
	% 	xlabel('y[m]');
	% 	ylabel('z[m]');
	% 	grid on;
	% 	axis equal;
	% print("output/defspar.png","-dpng","-r100")
	
	for cv=1:length(M_divp)%セクション番号1-4
		
		printf("generating ply [%d section] \n",cv);
		
		if(S_mat_swang(dp_num*(ci-1)+M_divp(1,cv)-1,1)>ID_swang(dp_num*(ci-1)+M_divp(1,cv)-1,1))%理想撓みより大きければプライを増す
			do

				psign=0;

				if(cv==1)%最内セクションだった時
					if(zply_c(ci,1)<Cply_lim)%集中積層限界数以下だった時
						disp("ply_ex");
						r_ang=Rply_ang(1,zply_c(ci,1)+1);%関数の入れ子エラー防止策
						ply_data{ci}=ply_ex(ply_data{ci},1,r_ang,C_T,mand_data(ci,1)./4.*1000);
						zply_c(ci,1)++;
					else%限界数を超えていた時
						disp("ply_ins");
						ply_ins_end=mand_data(ci,1)./4.*1000;
						insp=xplyc(ci,cv)+4;
						ply_data{ci}=ply_ins(ply_data{ci},insp,C_T,ply_ins_end);
					endif
				else%それ以外のセクションだった時

					cp=zply_c(ci,cv)+xplyc(ci,cv)+4;%伸ばすべき層番号

					if(cp>rows(ply_data{ci})-1&&cpdn<lpdn)%層番号が最外層に来てしまったときは最内セクションの積層を変更する

						if(zply_c(ci,1)<Cply_lim)%集中積層限界数以下だった時
							disp("re_ply_ex");
							r_ang=Rply_ang(1,zply_c(ci,1)+1);%関数の入れ子エラー防止策
							ply_data{ci}=ply_ex(ply_data{ci},1,r_ang,C_T,mand_data(ci,1)./4.*cv.*1000);
							zply_c(ci,1)++;
						else%限界数を超えていた時
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
					lpdn=200;%適当に大きい数
				else
					lpdn=rows(ply_data{ci-1});
				endif
				
				cpdn=rows(ply_data{ci});

				if cv==1&&cpdn>=lpdn%翼端に行くにつれ積層数が減るようにいろいろ
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

%積層構成を表示
for i=1:w_div_c
	ply_disp(ply_data{i},i);
endfor

%手動修正
disp("手動修正");
ss=input("1:yes 0:no ---");
if(ss==1)
	do
	ss2=input("１：修正　２：部分積層増加　３：一周積層増加　４：積層減少 ---");

	if ss2==1
		rw=input("修正翼番号=");
		rs=input("層番号=");
		rl=input("再定義範囲 1:1-1 2:1-2 3:1-3 4:1-4 ---");
		ply_data{rw}=ply_leng(ply_data{rw},rs,mand_data(rw,1)/4*rl*1000);

	elseif ss2==2
		rw=input("増加翼番号=");
		ra=input("積層角度=");
		rl=input("定義範囲 1:1-1 2:1-2 3:1-3 4:1-4 ---");

		ply_data{rw}=ply_ex(ply_data{rw},1,ra,C_T,mand_data(rw,1)./4*rl*1000)

	elseif ss2==3
		rw=input("増加翼番号=");
		rs=input("層番号=");
		rl=input("定義範囲 1:1-1 2:1-2 3:1-3 4:1-4 ---");

		ply_data{rw}=ply_ins(ply_data{rw},rs,C_T,mand_data(rw,1)./4*rl*1000);

	elseif ss2==4
		rw=input("減少翼番号=");
		rs=input("層番号=");
		ply_data{rw}(rs,:)=[];
		for j=2:rows(ply_data{rw})
			ply_data{rw}(j,1)=ply_data{rw}(j-1,1)+C_T(lookup(C_T(:,1),ply_data{rw}(j,2)),2)*2;
		end

	endif

	%積層構成を表示
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

	%桁が通るか最終確認

	for i=1:w_div_c
		ps_gap(i,1)=(real_thn_xcp(i,1)*1000-ply_data{i}(rows(ply_data{i}),1))/2-p_thn;
		if(ps_gap(i,1)>=lim_gap)
			printf("[%d翼]　桁-プランク間余裕 : %f \n",i,ps_gap(i,1));
		else
			printf("[%d翼]　限界余裕以下　翼定義不可能 : %f \n",i,ps_gap(i,1));
		end
	end

	W_Spar=W_Spar
	W_wing=W_wing

	ssign=input("next:1 end:0 ---")
		
	until(ssign==0)

endif	



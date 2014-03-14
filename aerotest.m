%まず楕円翼として曲げモーメントを求める
disp("基準スパンの0.8~1.4で探索")
span=input("基準スパン(m)=");
h_span=span/2;
d_s=h_span/100/2;
%ΔSだよ

%たわみ量を決定する
do
	disp("たわみ量は固定で探索")
	max_tawami=input("最大たわみ(m)=");
	%スパン2mあたり0.154mくらいでいいんじゃないかな
	alpha=max_tawami/(h_span*h_span);
	a=0;
	d_count=1;
	f_v=0;

	%線素の端点のy座標を求める（原点なし）
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
	
	%線素の端点のｚ座標を求める（原点なし）
	do
		line_e(d_count,2)=alpha*((line_e(d_count,1))^2);
		d_count=d_count+1;
	
	until(d_count>100)

	%disp(line_e);
	
	%線素のそれぞれの長さを求める（１００分割なので１００個）
	line_e_d(1,1)=sqrt((line_e(1,1))^2+(line_e(1,2))^2);

	d_count=2;

	do
		line_e_d(d_count,1)=sqrt((line_e(d_count,1)-line_e(d_count-1,1))^2+(line_e(d_count,2)-line_e(d_count-1,2))^2);
		d_count=d_count+1;
	
	until(d_count>100)

	%disp(line_e_d);

	%test2=rows(line_e_d)
	
	printf("正規最大たわみ(m)=%f\n",line_e(100,2));
	sign=input("ok:1 return:0 ---");
	
until(sign==1)

%コントロールポイント（上半角を考慮した時）の座標を求める
line_e_cp(1,1)=line_e(1,1)/2;
line_e_cp(1,2)=line_e(1,2)/2;
d_count=2;
do
	line_e_cp(d_count,1)=(line_e(d_count,1)-line_e(d_count-1,1))+line_e(d_count-1,1);
	line_e_cp(d_count,2)=(line_e(d_count,2)-line_e(d_count-1,2))+line_e(d_count-1,2);
	d_count=d_count+1;
	
until(d_count==101)

%各線素における局所上半角（ラジアン）を求める
line_e_ang(1,1)=atan(line_e(1,2)/line_e(1,1));
d_count=2;
do
	line_e_ang(d_count,1)=atan((line_e(d_count,2)-line_e(d_count-1,2))/(line_e(d_count,1)-line_e(d_count-1,1)));
	d_count=d_count+1;

until(d_count==101)

%上半角を考慮しなかった時のCP（のy座標のみ）を求める
d_count=1;
do
	nang_cp(d_count,1)=d_s+(d_count-1)*2*d_s;
	d_count=d_count+1;

until(d_count==101)

%翼根の基準曲げモーメントを楕円翼を用いて算出する
do
	weight=input("全備重量(kg)=");
	d_vel=input("設計機速(m/s)=");
	a_density=1.164;
	nu=1.604*10^(-5);

	%楕円循環と仮定して循環分布を求める
	root_gamma=4*weight*9.81/pi/a_density/d_vel/span;
	d_count=1;
	do
		d_gamma(d_count,1)=sqrt(root_gamma^2-(nang_cp(d_count,1)*root_gamma/span*2)^2);
		d_count=d_count+1;
		
	until(d_count==101)

	%循環分布を上半角を考慮した時のCPに適応して揚力とモーメントをもとめる
	%まず局所揚力を求める
	d_count=1;
	do
		d_l(d_count,1)=a_density*d_vel*d_gamma(d_count,1)*cos(line_e_ang(d_count,1))*line_e_d(d_count,1);
		d_count=d_count+1;

	until(d_count==101)
	%line_e_dはd_sでも問題ない？
	e_l=2*sum(d_l(:,1));
	%半翼分しか算出していないので二倍
	printf("上半角を考慮した時の揚力（N）=%f\n",e_l);
	printf("重量換算後(kg)=%f\n",e_l/9.81);
	
	sign=input("ok:1 return:0 ---");

until(sign==1)

%局所モーメントを求める（材力におけるモーメントではない）
d_count=1;
do
	d_b(d_count,1)=a_density*d_vel*d_gamma(d_count,1)*(line_e_cp(d_count,1)*cos(line_e_ang(d_count,1))+line_e_cp(d_count,2)*sin(line_e_ang(d_count,1)))*line_e_d(d_count,1);
	d_count=d_count+1;

until(d_count==101)

bend_m=sum(d_b(:,1));
bend_m=bend_m*3*pi/2;
%ここまで楕円翼と仮定した計算

%誘導抗力解析（ウイングレットは含まない）
div_count=input("解析分割数 ---");
if(div_count!=0)
	for i=1:div_count
		s_span(i,1)=span*(0.8+0.6/(div_count-1)*(i-1));
	endfor

	eq_i=(div_count-1)/3+1;


	%指定範囲で誘導抗力の値を算出
	for s_count=1:div_count
		printf("count %d/%d \n",s_count,div_count);
		%再度スパンとたわみを定義して本計算の準備を行う（今までの楕円翼の時のデータは曲げモーメントを除き失われる）
		span=s_span(s_count,1);
		h_span=span/2;
		d_s=h_span/100/2;
		%ΔSだよ
		
		%たわみ量は楕円翼から固定
		alpha=max_tawami/(h_span*h_span);
		a=0;
		d_count=1;
		f_v=0;

		%線素の端点のy座標を求める（原点なし）
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
		
		%線素の端点のｚ座標を求める（原点なし）
		do
			line_e(d_count,2)=alpha*((line_e(d_count,1))^2);
			d_count=d_count+1;
		
		until(d_count>100)

		%disp(line_e);
		
		%線素のそれぞれの長さを求める（１００分割なので１００個）
		line_e_d(1,1)=sqrt((line_e(1,1))^2+(line_e(1,2))^2);

		d_count=2;

		do
			line_e_d(d_count,1)=sqrt((line_e(d_count,1)-line_e(d_count-1,1))^2+(line_e(d_count,2)-line_e(d_count-1,2))^2);
			d_count=d_count+1;
		
		until(d_count>100)

		%disp(line_e_d);

		%test2=rows(line_e_d)

		%コントロールポイント（上半角を考慮した時）の座標を求める
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

		%各線素における局所上半角（ラジアン）を求める
		line_e_ang(1,1)=atan(line_e(1,2)/line_e(1,1));
		d_count=2;
		do
			line_e_ang(d_count,1)=atan((line_e(d_count,2)-line_e(d_count-1,2))/(line_e(d_count,1)-line_e(d_count-1,1)));
			d_count=d_count+1;

		until(d_count==101)
		%dlmwrite( "ang.txt",line_e_ang(:,1),"delimiter","***","newline","\n");

		%上半角を考慮しなかった時のCP（のy座標のみ）を求める
		d_count=1;
		do
			nang_cp(d_count,1)=d_s+(d_count-1)*2*d_s;
			d_count=d_count+1;

		until(d_count==101)

		%Ciを求める
		disp("calculating Ci.");
		d_count=1;
		do
			c_i(d_count,1)=2*cos(line_e_ang(d_count,1))*d_s;
			d_count=d_count+1;

		until(d_count==101)

		%Biを求める
		disp("calculating Bi.");
		d_count=1;
		do
			b_i(d_count,1)=3*pi/2*(line_e_cp(d_count,1)*cos(line_e_ang(d_count,1))+line_e_cp(d_count,2)*sin(line_e_ang(d_count,1)))*d_s;
			d_count=d_count+1;
			
		until(d_count==101)

		%yd_ijとかを求める
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

		%R2_ijとかを求める
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

		%Q_ijを求める
		disp("calculating Q_ij.");
		for i=1:100
			for j=1:100
				Q_ij(i,j)=-1/2/pi*(((yd_ij(i,j)-line_e_d(j,1)./2)./R2_pij(i,j)-1*(yd_ij(i,j)+line_e_d(j,1)./2)./R2_mij(i,j))*cos(line_e_ang(i,1)-line_e_ang(j,1))+(zd_ij(i,j)./R2_pij(i,j)-zd_ij(i,j)./R2_mij(i,j))*sin(line_e_ang(i,1)-line_e_ang(j,1))+((ydd_ij(i,j)-line_e_d(j,1)./2)./Rd2_mij(i,j)-1*(ydd_ij(i,j)+line_e_d(j,1)./2)./Rd2_pij(i,j))*cos(line_e_ang(i,1)+line_e_ang(j,1))+(zdd_ij(i,j)./Rd2_mij(i,j)-zdd_ij(i,j)./Rd2_pij(i,j))*sin(line_e_ang(i,1)+line_e_ang(j,1)));
			endfor
		endfor
		%dlmwrite( "qij.txt",Q_ij(:,1),"delimiter","***","newline","\n");

		%A_ijを求める
		disp("calculating A_ij.");
		for i=1:100
			for j=1:100
				A_ij(i,j)=pi*Q_ij(i,j)*line_e_d(i,1)/2;
			endfor
		endfor

		%dlmwrite( "aij.txt",A_ij(:,1),"delimiter","***","newline","\n");

		%最適化行列を作る
		%まずA_ijからなる一部分を作る
		for i=1:100
			for j=1:100
				Op_p1(i,j)=A_ij(i,j)+A_ij(j,i);
			endfor
		endfor

		%サイドの一部分を作る
		Op_p2=[-c_i -b_i];

		%下の一部分を作る
		Op_p3=[-(c_i)' 0 0;-(b_i)' 0 0];

		%全部合わせる
		Op_mat=[Op_p1 Op_p2;Op_p3];

		%0000-1betaの行列を作る
		for i=1:100
			Op_zero(i,1)=0;
		endfor
		%設計揚力
		del=weight*9.81;
		Op_r=[Op_zero;-del;-bend_m];
		a=[1.0 2.0;3.0 4.0];

		%確認用に行列を出力
		%dlmwrite( "Opr.txt",Op_mat(:,1:10),"delimiter","***","newline","\n");

		%本計算giを求める
		disp("calculating g_i.");
		g_i=Op_mat\Op_r;

		%giを元に循環分布を求める
		for i=1:100
			Op_gamma(i,1)=g_i(i,1)/2/a_density/d_vel;
		endfor

		%垂直方向誘導速度
		disp("calculating Vn_i.");
		for i=1:100
			for j=1:100
				Vn_i_d(j,i)=Q_ij(i,j)*Op_gamma(j,1)./2;
			endfor
			Vn_i(i,1)=sum(Vn_i_d(:,i));
			%printf("%d/100\n",i);
		endfor

		%誘導抗力
		disp("calculating di_d.");
		for i=1:100
			di_d(i,1)=Op_gamma(i,1)*a_density*Vn_i(i,1)*line_e_d(i,1);
		endfor
		%スパンごとに誘導抗力を格納
		di(s_count,1)=2*sum(di_d(:,1));
		
	endfor

	%指定範囲の誘導抗力の値を出力
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
		
		%-r100は解像度を100dpiにする．
	cd("design_data")
	cd(pj_name)
	print("Di_span.png",'-dpng','-r100')
	cd(now_work)
	
endif

%グラフからスパンを選び正式に計算
%再度スパンとたわみを定義して本計算の準備を行う（今までの楕円翼の時のデータは曲げモーメントを除き失われる）
%翼素分割数は可変
do

	dw_divn=50;%分割翼分割数
	disp("これより本計算です");
	span=input("本スパン(m)=");
	w_div_c=input("主翼分割数=");
	we_div=w_div_c*dw_divn;%主翼平面部分割数

	for i=1:w_div_c
		disp(i);
		divw_span(i,1)=input("分割翼長さ(m)=");
	endfor

	h_span=span/2;
	% d_s=h_span/we_div/2
	%使うのやめたよ
	d_s(:,1)=divw_span(:,1)./dw_divn;

	%たわみ量を決定する
	do
		max_tawami=input("最大たわみ(m)=");
		%スパン2mあたり0.154mくらいでいいんじゃないかな
		alpha=max_tawami/(h_span*h_span);
		a=0;
		d_count=1;
		f_v=0;

		%線素の端点のy座標を求める（原点なし）
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
		
		%線素の端点のｚ座標を求める（原点なし）
		do
			line_e(d_count,2)=alpha*((line_e(d_count,1))^2);
			d_count=d_count+1;
		
		until(d_count>we_div)

		%disp(line_e);
		
		%線素のそれぞれの長さを求める（１００分割なので１００個）
		line_e_d(1,1)=sqrt((line_e(1,1))^2+(line_e(1,2))^2);

		d_count=2;

		do
			line_e_d(d_count,1)=sqrt((line_e(d_count,1)-line_e(d_count-1,1))^2+(line_e(d_count,2)-line_e(d_count-1,2))^2);
			d_count=d_count+1;
		
		until(d_count>we_div)

		%disp(line_e_d);

		%test2=rows(line_e_d)
		
		printf("正規最大たわみ(m)=%f\n",line_e(we_div,2));
		sign=input("ok:1 return:0 ---");
		
	until(sign==1)

	%コントロールポイント（上半角を考慮した時）の座標を求める
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

	%各線素における局所上半角（ラジアン）を求める
	line_e_ang(1,1)=atan(line_e(1,2)/line_e(1,1));
	d_count=2;
	do
		line_e_ang(d_count,1)=atan((line_e(d_count,2)-line_e(d_count-1,2))/(line_e(d_count,1)-line_e(d_count-1,1)));
		d_count=d_count+1;

	until(d_count==we_div+1)
	%dlmwrite( "ang.txt",line_e_ang(:,1),"delimiter","***","newline","\n");
	
	%ここからウイングレット作成部
	%ウイングレットを作成する
	wl_ang=input("ウイングレットの取り付け角度（たわみなし時、上方プラス）=");
	wl_span=input("ウイングレットの長さ（m）=");
	
	if(wl_span==0)
	wld_c=0;
	endif
	
	if(wl_span>0)
		%ウイングレットの分割数を元に翼素の大きさを求める、ここでは２０分割
		wld_c=20;
		wl_d=wl_span/wld_c;
		%上半角を考慮したときの絶対取り付け角を算出する
		wl_ang_a=line_e_ang(we_div,1)+wl_ang/180*pi;
		%翼素１つあたりのｙ，ｚ座標変化を求める
		wl_d_y=wl_d*cos(wl_ang_a);
		wl_d_z=wl_d*sin(wl_ang_a);
		%１：２０までの翼素の端点座標を求める
		for i=1:wld_c
			wl_e(i,1)=line_e(we_div,1)+wl_d_y*i;
			wl_e(i,2)=line_e(we_div,2)+wl_d_z*i;
		endfor
		%1:20までのCP座標を求める
		wl_e_cp(1,1)=line_e(we_div,1)+wl_d_y/2;
		wl_e_cp(1,2)=line_e(we_div,2)+wl_d_z/2;
		
		for i=2:wld_c
			wl_e_cp(i,1)=wl_e_cp(i-1,1)+wl_d_y;
			wl_e_cp(i,2)=wl_e_cp(i-1,2)+wl_d_z;
		endfor
		
		%翼素のデータを格納している変数にウィングレットのデータを格納する
		line_e((we_div+1):(wld_c+we_div),1)=wl_e(1:wld_c,1);
		line_e((we_div+1):(wld_c+we_div),2)=wl_e(1:wld_c,2);
		line_e_cp((we_div+1):(wld_c+we_div),1)=wl_e_cp(1:wld_c,1);
		line_e_cp((we_div+1):(wld_c+we_div),2)=wl_e_cp(1:wld_c,2);
		line_e_ang((we_div+1):(wld_c+we_div),1)=wl_ang_a;
		line_e_d((we_div+1):(wld_c+we_div),1)=wl_d;
		
	endif
	
	%上半角を考慮しなかった時のCP（のy座標のみ）を求める
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

	%Ciを求める
	disp("calculating Ci.");
	d_count=1;
	do
		c_i(d_count,1)=2*cos(line_e_ang(d_count,1))*line_e_d(d_count,1)/2;
		d_count=d_count+1;

	until(d_count==wld_c+we_div+1)

	%Biを求める
	disp("calculating Bi.");
	d_count=1;
	do
		b_i(d_count,1)=3*pi/2*(line_e_cp(d_count,1)*cos(line_e_ang(d_count,1))+line_e_cp(d_count,2)*sin(line_e_ang(d_count,1)))*line_e_d(d_count,1)/2;
		d_count=d_count+1;
		
	until(d_count==wld_c+we_div+1)

	%yd_ijとかを求める
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

	%R2_ijとかを求める
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

	%Q_ijを求める
	disp("calculating Q_ij.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Q_ij(i,j)=-1/2/pi*(((yd_ij(i,j)-line_e_d(j,1)./2)./R2_pij(i,j)-1*(yd_ij(i,j)+line_e_d(j,1)./2)./R2_mij(i,j))*cos(line_e_ang(i,1)-line_e_ang(j,1))+(zd_ij(i,j)./R2_pij(i,j)-zd_ij(i,j)./R2_mij(i,j))*sin(line_e_ang(i,1)-line_e_ang(j,1))+((ydd_ij(i,j)-line_e_d(j,1)./2)./Rd2_mij(i,j)-1*(ydd_ij(i,j)+line_e_d(j,1)./2)./Rd2_pij(i,j))*cos(line_e_ang(i,1)+line_e_ang(j,1))+(zdd_ij(i,j)./Rd2_mij(i,j)-zdd_ij(i,j)./Rd2_pij(i,j))*sin(line_e_ang(i,1)+line_e_ang(j,1)));
		endfor
	endfor
	%dlmwrite( "qij.txt",Q_ij(:,1),"delimiter","***","newline","\n");

	%A_ijを求める
	disp("calculating A_ij.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			A_ij(i,j)=pi*Q_ij(i,j)*line_e_d(i,1)/2;
		endfor
	endfor

	%dlmwrite( "aij.txt",A_ij(:,1),"delimiter","***","newline","\n");

	%最適化行列を作る
	%まずA_ijからなる一部分を作る
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Op_p1(i,j)=A_ij(i,j)+A_ij(j,i);
		endfor
	endfor

	%サイドの一部分を作る
	Op_p2=[-c_i -b_i];

	%下の一部分を作る
	Op_p3=[-(c_i)' 0 0;-(b_i)' 0 0];

	%全部合わせる
	Op_mat=[Op_p1 Op_p2;Op_p3];

	%0000-1betaの行列を作る
	for i=1:(wld_c+we_div)
		Op_zero(i,1)=0;
	endfor
	%設計揚力
	del=weight*9.81;
	Op_r=[Op_zero;-del;-bend_m];
	a=[1.0 2.0;3.0 4.0];

	%確認用に行列を出力
	%dlmwrite( "Opr.txt",Op_mat(:,1:10),"delimiter","***","newline","\n");

	%本計算giを求める
	disp("calculating g_i.");
	g_i=Op_mat\Op_r;

	%giを元に循環分布を求める
	for i=1:(wld_c+we_div)
		Op_gamma(i,1)=g_i(i,1)/2/a_density/d_vel;
	endfor
	dlmwrite( "log/Opgamma.txt",Op_gamma(:,1),"delimiter","***","newline","\n");

	%垂直方向誘導速度
	disp("calculating Vn_i.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Vn_i_d(j,i)=Q_ij(i,j)*Op_gamma(j,1)./2;
		endfor
		Vn_i(i,1)=sum(Vn_i_d(:,i));
		%printf("%d/100\n",i);
	endfor

	%誘導角度
	disp("calculating i_ang.");
	for i=1:(wld_c+we_div)
		i_ang(i,1)=atan(Vn_i(i,1)/d_vel)*180/pi;
	endfor
	
	%合成速度計算
	disp("calculating comp_vel");
	for i=1:(wld_c+we_div)
		comp_vel(i,1)=sqrt(d_vel^2+Vn_i(i,1)^2);
	endfor

	%誘導抗力
	disp("calculating di_d.");
	for i=1:(wld_c+we_div)
		di_d(i,1)=Op_gamma(i,1)*Vn_i(i,1)*line_e_d(i,1);
	endfor
	di=2*a_density*sum(di_d(:,1));

	%局所揚力を求める(材力用,上半なし)
	disp("calculating d_l.");
	for i=1:(wld_c+we_div)
		d_l(i,1)=Op_gamma(i,1)*d_vel*a_density*line_e_d(i,1);
	endfor

	%せん断力を求める(ウイングレットなし)
	disp("calculating Q.");
	Q(we_div,1)=d_l(we_div,1);
	for i=(we_div-1):-1:1
		Q(i,1)=Q(i+1,1)+d_l(i,1);
	endfor
	%dlmwrite( "Q.txt",Q(:,1),"delimiter","***","newline","\n");

	if(wl_span>0)
		%ウイングレット部のせん断力を求める
		Q(wld_c+we_div,1)=d_l(wld_c+we_div,1);
		for i=(we_div-1+wld_c):-1:we_div+1
			Q(i,1)=Q(i+1,1)+d_l(i,1);
		endfor
	endif
	
	%局所モーメントを求める(mmに変換,ウイングレットなし)
	disp("calculating d_m.");
	d_m(we_div,1)=Q(we_div,1)*line_e_d(we_div,1)/2*1000;
	for i=(we_div-1):-1:1
		d_m(i,1)=(Q(i,1)+Q(i+1))*line_e_d(i,1)/2*1000;
	endfor
	%dlmwrite( "dM.txt",d_m(:,1),"delimiter","***","newline","\n");
	
	if(wl_span>0)
		%ウイングレット部の局所モーメントを求める
		d_m(wld_c+we_div,1)=Q(wld_c+we_div,1)*line_e_d(wld_c+we_div,1)/2*1000;
		for i=(we_div-1+wld_c):-1:(we_div+1)
			d_m(i,1)=(Q(i,1)+Q(i+1))*line_e_d(i,1)/2*1000;
		endfor
		
		%ウイングレットのモーメントを求めて翼端にマイナス方向モーメント荷重として計算する
		M(wld_c+we_div,1)=d_m(wld_c+we_div,1);
		for i=(we_div-1+wld_c):-1:we_div+1
			M(i,1)=M(i+1,1)+d_m(i,1);
		endfor
		%M(we_div+1,1)が翼端にかかると考えられるモーメント荷重
	endif

	if(wl_span==0)	
		M(we_div+1,1)=0;
	endif
		
	%モーメントを求める(N*mm)
	disp("calculating M.");
	M(we_div,1)=d_m(we_div,1)-M(we_div+1,1);
	for i=(we_div-1):-1:1
		M(i,1)=M(i+1,1)+d_m(i,1);
	endfor
	%M(1,1)が翼根部曲げモーメント
	%dlmwrite( "M.txt",M(:,1),"delimiter","***","newline","\n");

	%揚力計算
	disp("calculating L.");
	for i=1:(wld_c+we_div)
		dd_l(i,1)=Op_gamma(i,1)*cos(line_e_ang(i,1))*line_e_d(i,1)/2;
	endfor
	re_def_l=4*a_density*d_vel*sum(dd_l(:,1));

	%楕円循環と仮定して循環分布を求める(循環分布を確認するため)
	root_gamma=4*weight*9.81/pi/a_density/d_vel/(span+wl_span*2);
	d_count=1;
	do
		d_gamma(d_count,1)=sqrt(root_gamma^2-(nang_cp(d_count,1)*root_gamma/(span+wl_span*2)*2)^2);
		d_count=d_count+1;
		
	until(d_count==wld_c+we_div+1)


	%%%試験のため循環分布を多角形近似する
	for i=1:5
		gamma_d(i,1)=(Op_gamma(i*40+1,1)-Op_gamma((i-1)*40+1,1))/40;
		
		for j=1:40
			Hex_gamma((i-1)*40+j,1)=Op_gamma((i-1)*40+1,1)+gamma_d(i,1)*(j-1);
		endfor
	endfor

	Hex_gamma(201:220,1)=Op_gamma(201:220,1);


	%垂直方向誘導速度
	disp("calculating Hex_Vn_i.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Hex_Vn_i_d(j,i)=Q_ij(i,j)*Hex_gamma(j,1)./2;
		endfor
		Hex_Vn_i(i,1)=sum(Hex_Vn_i_d(:,i));
		%printf("%d/100\n",i);
	endfor

	%誘導角度
	disp("calculating Hex_i_ang.");
	for i=1:(wld_c+we_div)
		Hex_i_ang(i,1)=atan(Hex_Vn_i(i,1)/d_vel)*180/pi;
	endfor
	
	%合成速度計算
	disp("calculating Hex_comp_vel");
	for i=1:(wld_c+we_div)
		Hex_comp_vel(i,1)=sqrt(d_vel^2+Hex_Vn_i(i,1)^2);
	endfor

	%誘導抗力
	disp("calculating Hex_di_d.");
	for i=1:(wld_c+we_div)
		Hex_di_d(i,1)=Hex_gamma(i,1)*Hex_Vn_i(i,1)*line_e_d(i,1);
	endfor
	Hex_di=2*a_density*sum(Hex_di_d(:,1));



	%計算値を出力
	disp("計算結果");
	printf("機速(m/s)=%f\n",d_vel);
	printf("スパン(m)=%f\n",span);
	printf("誘導抗力(N)=%f\n",di);
	printf("多角形誘導抗力(N)=%f\n",Hex_di);

	%出力
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
			
		%-r100は解像度を100dpiにする．
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
%循環を元に平面形を作る

%翼型解析データの読み込み
%読込フォルダへの移動

now_work=pwd;

%read_d=input("翼型データのあるディレクトリ名---");
read_d="FX76MP149_e66";
cd (read_d)

%-----レイノルズリスト読み込み
Relist={50000,'0.050';100000,'0.100'; 150000,'0.150'; 200000,'0.200'; 250000,'0.250'; 300000,'0.300';350000,'0.350'; 400000,'0.400';450000,'0.450'; 500000,'0.500';550000,'0.550';600000,'0.600';};
foillist=["FX76MP149";'fab_h_b';'fab_f_d';'fab_d_f';'fab_b_h';"e66"];

for n=1:6
%-----nは順にfoila';'fab_h_b';'fab_f_d';'fab_d_f';'fab_b_h';'foilb

	for m=1:12
	%-----mは順に50000,100000,150000,200000,250000,300000,350000,400000,450000,500000,550000,600000

		foilname=foillist(n,:);
		datafile_a='_T1_Re';
		Re=Relist{12+m};
		datafile_b='_M0.00_N9.0.txt';
		dataname=strcat(foillist,datafile_a,Re,datafile_b);

		fpr=fopen(dataname(n,:),'r');
		do
			buffer=fscanf(fpr,'%c',[1,1]);
		until(buffer(1,1)=='-')
		buffer=fgetl(fpr);
		adata{n,m}=(fscanf(fpr,'%f',[10,100]))';
		fclose(fpr);
	endfor
clear buffer 
endfor

%もとの作業ディレクトリへ
cd (now_work)

%翼型選定の前に循環分布整える

%循環やその他グラフに於いて翼端で大きな値の振れが生じているので端点の値を無視してつなげるべきである
%なおウイングレット開始部において値の振れが生じているので平面翼とウイングレットの循環は全く別のものとして扱う(変化部付近6点は特異点として除外する)

do

	%循環編集用に循環をコピーする
	ed_gamma=Op_gamma;

	nc_span=input("翼型、翼弦長無変化部スパン（半翼分）（ｍ）=");
	for nc_p=1:we_div
		if(nang_cp(nc_p,1)>=nc_span)
		break;
		endif
	endfor

	%循環の頂点部分をカットする
	for i=1:(nc_p-1)
		ed_gamma(i,1)=ed_gamma(nc_p,1);
	endfor

	%途中矩形部の範囲指定を行う
	%とりあえず循環変化率のグラフと変化率の変化率のグラフを出す
	del_gamma=gradient(Op_gamma(1:we_div,1));
	del_2_gamma=gradient(del_gamma(:,1));

	figure(6);
		subplot(2,1,1);
			plot(nang_cp(1:we_div,1),del_gamma(:,1),"-");
			xlabel('y[m]');
			ylabel('d_gamma');
			grid on;
			
		subplot(2,1,2);
			plot(nang_cp(1:we_div,1),del_2_gamma(:,1),"-");
			xlabel('y[m]');
			ylabel('d2_gamma');
			grid on;
			
		%-r100は解像度を100dpiにする．
	cd("design_data")
	cd(pj_name)
	print("del_gamma.png","-dpng","-r100")
	cd(now_work)

	%線形近似する部分の範囲を指定する

	rec_start=input("途中矩形部開始位置(m)=");
	rec_end=input("途中矩形部終了(m)=");

	for rec_sp=1:we_div
		if(nang_cp(rec_sp,1)>=rec_start)
		break;
		endif
	endfor

	for rec_ep=1:we_div
		if(nang_cp(rec_ep,1)>=rec_end)
		break;
		endif
	endfor

	alpha_g=(ed_gamma(rec_ep,1)-ed_gamma(rec_sp,1))/(nang_cp(rec_ep,1)-nang_cp(rec_sp,1));
	for i=(rec_sp+1):(rec_ep-1)
		ed_gamma(i,1)=alpha_g*(nang_cp(i,1)-nang_cp(rec_sp,1))+ed_gamma(rec_sp,1);
	endfor
	%線形近似終了

	%接合部処理
	%接合部周り６点を除外、スプラインで保管する
	if wld_c>0
		for i=2:-1:0
			ed_gamma(we_div-i,1)=spline(nang_cp(1:(we_div-3),1),Op_gamma(1:(we_div-3),1),nang_cp(we_div-i,1));

		endfor

		for i=1:3
			ed_gamma(we_div+i,1)=spline(nang_cp((we_div+4):we_div+wld_c,1),Op_gamma((we_div+4):we_div+wld_c,1),nang_cp(we_div+i,1));

		endfor
	endif
	%循環編集終了
	%翼端処理は後で追加する

	disp("recalating wing");

	%翼性能再計算
	%垂直方向誘導速度
	disp("calculating Vn_i.");
	for i=1:(wld_c+we_div)
		for j=1:(wld_c+we_div)
			Vn_i_d_ed(j,i)=Q_ij(i,j)*ed_gamma(j,1)./2;
		endfor
		Vn_i_ed(i,1)=sum(Vn_i_d_ed(:,i));
		%printf("%d/100\n",i);
	endfor

	%誘導角度
	disp("calculating i_ang.");
	for i=1:(wld_c+we_div)
		i_ang_ed(i,1)=atan(Vn_i_ed(i,1)/d_vel)*180/pi;
	endfor

	%誘導抗力
	disp("calculating di_d.");
	for i=1:(wld_c+we_div)
		di_d_ed(i,1)=ed_gamma(i,1)*Vn_i_ed(i,1)*line_e_d(i,1);
	endfor
	di_ed=2*a_density*sum(di_d_ed(:,1));

	%局所揚力を求める(材力用,上半なし)
	disp("calculating d_l.");
	for i=1:(wld_c+we_div)
		d_l_ed(i,1)=ed_gamma(i,1)*d_vel*a_density*line_e_d(i,1);
	endfor

	%せん断力を求める(ウイングレットなし)
	disp("calculating Q.");
	Q_ed(we_div,1)=d_l_ed(we_div,1);
	for i=(we_div-1):-1:1
		Q_ed(i,1)=Q_ed(i+1,1)+d_l_ed(i,1);
	endfor
	%dlmwrite( "Q.txt",Q(:,1),"delimiter","***","newline","\n");

	if(wl_span>0)
		%ウイングレット部のせん断力を求める
		Q_ed(wld_c+we_div,1)=d_l_ed(wld_c+we_div,1);
		for i=(we_div-1+wld_c):-1:we_div+1
			Q_ed(i,1)=Q_ed(i+1,1)+d_l_ed(i,1);
		endfor
	endif

	%局所モーメントを求める(mmに変換,ウイングレットなし)
	disp("calculating d_m.");
	d_m_ed(we_div,1)=Q_ed(we_div,1)*line_e_d(we_div,1)/2*1000;
	for i=(we_div-1):-1:1
		d_m_ed(i,1)=(Q_ed(i,1)+Q_ed(i+1))*line_e_d(i,1)/2*1000;
	endfor
	%dlmwrite( "dM.txt",d_m(:,1),"delimiter","***","newline","\n");

	if(wl_span>0)
		%ウイングレット部の局所モーメントを求める
		d_m_ed(wld_c+we_div,1)=Q(wld_c+we_div,1)*line_e_d(wld_c+we_div,1)/2*1000;
		for i=(we_div-1+wld_c):-1:(we_div+1)
			d_m_ed(i,1)=(Q_ed(i,1)+Q_ed(i+1))*line_e_d(i,1)/2*1000;
		endfor
		
		%ウイングレットのモーメントを求めて翼端にマイナス方向モーメント荷重として計算する
		M_ed(wld_c+we_div,1)=d_m_ed(wld_c+we_div,1);
		for i=(we_div-1+wld_c):-1:we_div+1
			M_ed(i,1)=M_ed(i+1,1)+d_m_ed(i,1);
		endfor
		%M(we_div+1,1)が翼端にかかると考えられるモーメント荷重
	endif

	if(wl_span==0)	
		M_ed(we_div+1,1)=0;
	endif
		
	%モーメントを求める(N*mm)
	disp("calculating M.");
	M_ed(we_div,1)=d_m_ed(we_div,1)-M_ed(we_div+1,1);
	for i=(we_div-1):-1:1
		M_ed(i,1)=M_ed(i+1,1)+d_m_ed(i,1);
	endfor
	%M(1,1)が翼根部曲げモーメント
	%dlmwrite( "M.txt",M(:,1),"delimiter","***","newline","\n");

	%揚力計算
	disp("calculating L.");
	for i=1:(wld_c+we_div)
		dd_l_ed(i,1)=ed_gamma(i,1)*cos(line_e_ang(i,1))*line_e_d(i,1)/2;
	endfor
	re_def_l_ed=4*a_density*d_vel*sum(dd_l_ed(:,1));

	%計算値を出力
	disp("再計算結果");
	printf("上半考慮揚力(N)=%f\n",re_def_l_ed);
	printf("誘導抗力(N)=%f\n",di_ed);
	printf("翼根曲げモーメント(N*mm)=%f\n",M_ed(1,1));
		
	figure(7);
		plot(nang_cp(:,1),ed_gamma(:,1),"-;ed_gamma;",nang_cp(:,1),Op_gamma(:,1),"-;Op_gamma;");
		xlabel('y[m]');
		ylabel('gamma');
		xlim([0,span/2+wl_span+1]);
		axis("equal");
		grid on;
		
		%-r100は解像度を100dpiにする．
	cd("design_data")
	cd(pj_name)
	print("Op_ed_gamma.png","-dpng","-r100")
	cd(now_work)

	figure(8);
		subplot(2,1,1);
			plot(nang_cp(:,1),Vn_i_ed(:,1),"-");
			xlabel('y[m]');
			ylabel('Vn_i_ed[m/s]');
			xlim([0,span/2+wl_span+1]);
			grid on;
			
		subplot(2,1,2);
			plot(nang_cp(:,1),i_ang_ed(:,1),"-");
			xlabel('y[m]');
			ylabel('i_angle_ed');
			xlim([0,span/2+wl_span+1]);
			grid on;
			
		%-r100は解像度を100dpiにする．
	cd("design_data")
	cd(pj_name)
	print("Vn_iang_ed.png","-dpng","-r100")
	cd(now_work)
	
	sign=input("ok:1 return:0---");
	
until(sign==1)
%設計用循環準備完了

%簡略化のため垂直方向誘導速度を近似する
%中央翼では一定、それ以外では線形近似、接合部で分割

disp("approximation Vn_i&i_angle");

%線形近似部
%主翼
%近似位置定義
ap_p=ceil((rec_sp-nc_p)/2+nc_p);

alpha_v=(Vn_i_ed(we_div-9,1)-Vn_i_ed(ap_p,1))/(nang_cp(we_div-9,1)-nang_cp(ap_p,1));
for i=nc_p:we_div
	Vn_i_ed_ap(i,1)=alpha_v*(nang_cp(i,1)-nang_cp(ap_p,1))+Vn_i_ed(ap_p,1);
endfor

%中央翼
for i=1:nc_p
	% Vn_i_ed_ap(i,1)=sum(Vn_i_ed(1:nc_p,1))/nc_p;
	Vn_i_ed_ap(i,1)=Vn_i_ed_ap(nc_p,1);
		
endfor

%ウイングレット
alpha_vl=(Vn_i_ed(we_div+wld_c,1)-Vn_i_ed(we_div+8,1))/(nang_cp(we_div+wld_c,1)-nang_cp(we_div+8,1));
for i=(we_div+1):(we_div+wld_c)
	Vn_i_ed_ap(i,1)=alpha_vl*(nang_cp(i,1)-nang_cp(we_div+8,1))+Vn_i_ed(we_div+8,1);
endfor

%誘導角度再々計算
disp("calculating i_ang_ed_ap");
for i=1:(wld_c+we_div)
	i_ang_ed_ap(i,1)=atan(Vn_i_ed_ap(i,1)/d_vel)*180/pi;
endfor

%合成速度計算
disp("calculating comp_vel");
for i=1:(wld_c+we_div)
	comp_vel_ed_ap(i,1)=sqrt(d_vel^2+Vn_i_ed_ap(i,1)^2);
endfor

figure(9);
	subplot(2,1,1);
		plot(nang_cp(:,1),Vn_i_ed(:,1),"-;Vn_i_ed;",nang_cp(:,1),Vn_i_ed_ap(:,1),"-;Vn_i_ed_ap;");
		xlabel('y[m]');
		ylabel('Vn_i[m/s]');
		xlim([0,span/2+wl_span+1]);
		grid on;
		
	subplot(2,1,2);
		plot(nang_cp(:,1),i_ang_ed(:,1),"-;i_ang_ed;",nang_cp(:,1),i_ang_ed_ap(:,1),"-;i_ang_ed_ap;");
		xlabel('y[m]');
		ylabel('i_angle');
		xlim([0,span/2+wl_span+1]);
		grid on;
		
	%-r100は解像度を100dpiにする．
cd("design_data")
cd(pj_name)
print("Vn_iang_ed_ap.png","-dpng","-r100")
cd(now_work)


%平面形作成開始

%cl*cを定義
for i=1:(wld_c+we_div)
	c_cl(i,1)=2*ed_gamma(i,1)/comp_vel_ed_ap(i,1);
endfor

mt_ang=input("取り付け角=");

root_aoa=mt_ang-i_ang_ed_ap(1,1);

%Re数を450kと仮定し翼弦長を算出
root_chord=2*ed_gamma(1,1)/comp_vel_ed_ap(1,1)/(spline(adata{1,9}(:,1),adata{1,9}(:,2),root_aoa));

%Re数を算出
root_re=root_chord*comp_vel_ed_ap(1,1)/nu;
root_cl=ed_gamma(1,1)*2/comp_vel_ed_ap(1,1)/root_chord;

%取り付け角変更
%誘導迎え角分あらかじめ増やしておく
aoa=Re_spline(root_re,1,2,root_cl,1,adata);
mt_ang(1,1)=i_ang_ed_ap(1,1)+aoa;

disp("generating A wing");
%中央翼を定義
for m=1:10
	root_adata(1,m)=Re_spline(root_re,1,2,root_cl,m,adata);
endfor

for m=1:nc_p
	f_data(m,1:3)=[nang_cp(m,1),root_chord,1];
	f_data(m,4:13)=root_adata(1,:);
	re(m,1)=root_re;
	mt_ang(m,1)=mt_ang(1,1);
endfor

%第一テーパー部計算
disp("generating B wing");

%矩形部終了と同じ混合率で次の翼素を決定し、循環分布に合わせるための翼弦長を計算。翼弦長の式作成。代表レイノルズ数計算
%一つ目の翼弦長を計算
f_data(nc_p+1,2)=c_cl(nc_p+1,1)/f_data(nc_p,5);
alpha_t1=(f_data(nc_p+1,2)-f_data(nc_p,2))/(nang_cp(nc_p+1,1)-nang_cp(nc_p,1));
	%t_ratio1=c_cl(nc_p+1,1)/f_data(nc_p,5)/f_data(nc_p,2);
for t1=nc_p+1:ceil((rec_sp-nc_p)/2+nc_p)
	disp(t1);
	f_data(t1,2)=alpha_t1*(nang_cp(t1,1)-nang_cp(nc_p,1))+f_data(nc_p,2);
		%f_data(t1,2)=f_data(t1-1,2)*t_ratio1;
	cl(t1,1)=c_cl(t1,1)/f_data(t1,2);
	f_data(t1,1)=nang_cp(t1,1);
	re(t1,1)=comp_vel_ed_ap(t1,1)*f_data(t1,2)/nu;
	f_data(t1,3:8)=foil_mix(re(t1,1),cl(t1,1),(mt_ang(1,1)-i_ang_ed_ap(t1,1)),adata);
	mt_ang(t1,1)=mt_ang(1,1);
	if(f_data(t1,3)<0.8)
		break;
	endif
endfor

%第二テーパー部計算
disp("generating C wing");

%矩形部終了と同じ混合率で次の翼素を決定し、循環分布に合わせるための翼弦長を計算。翼弦長の式作成。代表レイノルズ数計算
	%t_ratio2=c_cl(t1+1,1)/f_data(t1,5)/f_data(t1,2);
f_data(t1+1,2)=c_cl(t1+1,1)/f_data(t1,5);
alpha_t2=(f_data(t1+1,2)-f_data(t1,2))/(nang_cp(t1+1,1)-nang_cp(t1,1));
for t2=t1+1:(rec_sp-1)
	disp(t2);
	f_data(t2,2)=alpha_t2*(nang_cp(t2,1)-nang_cp(t1,1))+f_data(t1,2);
		%f_data(t2,2)=f_data(t2-1,2)*t_ratio2;
	cl(t2,1)=c_cl(t2,1)/f_data(t2,2);
	f_data(t2,1)=nang_cp(t2,1);
	re(t2,1)=comp_vel_ed_ap(t2,1)*f_data(t2,2)/nu;
	f_data(t2,3:8)=foil_mix(re(t2,1),cl(t2,1),(mt_ang(1,1)-i_ang_ed_ap(t2,1)),adata);
	mt_ang(t2,1)=mt_ang(1,1);
	
endfor

%途中矩形部計算
disp("generating D wing");

%捻り下げ部
for m1=rec_sp:rec_ep
	disp(m1);
	f_data(m1,2)=f_data(t2,2);
	cl(m1,1)=c_cl(m1,1)/f_data(m1,2);
	f_data(m1,5)=cl(m1,1);
	re(m1,1)=comp_vel_ed_ap(m1,1)*f_data(m1,2)/nu;
	f_data(m1,1)=nang_cp(m1,1);
	%レイノルズ数、C、Cl、を固定しαを探索する
	f_data(m1,3:8)=aoa_def(re(m1,1),cl(m1,1),f_data(t2,3),adata);
	mt_ang(m1,1)=f_data(m1,4)+i_ang_ed_ap(m1,1);
	
endfor

%第三テーパー部計算
disp("generating E wing");
%矩形部終了と同じ混合率で次の翼素を決定し、循環分布に合わせるための翼弦長を計算。翼弦長の式作成。代表レイノルズ数計算
	%f_data(m1+1,2)=c_cl(m1+1,1)/f_data(m1,5);
	%alpha_t3=(f_data(m1+1,2)-f_data(m1,2))/(nang_cp(m1+1,1)-nang_cp(m1,1));	
t_ratio3=(c_cl(m1+1,1)/f_data(m1,5)/f_data(m1,2))^0.5;
for t3=m1+1:we_div
	disp(t3);
	f_data(t3,2)=f_data(t3-1,2)*t_ratio3;
		%f_data(t3,2)=alpha_t3*(nang_cp(t3,1)-nang_cp(m1,1))+f_data(m1,2);
	cl(t3,1)=c_cl(t3,1)/f_data(t3,2);
	f_data(t3,1)=nang_cp(t3,1);
	re(t3,1)=comp_vel_ed_ap(t3,1)*f_data(t3,2)/nu;
	f_data(t3,3:8)=foil_mix(re(t3,1),cl(t3,1),f_data(m1,4)+i_ang_ed_ap(m1,1)-i_ang_ed_ap(t3,1),adata);
	mt_ang(t3,1)=f_data(t3,4)+i_ang_ed_ap(t3,1);
	if(f_data(t3,3)<0.00)
		break;
	endif
	
endfor

%第四テーパー捻り下げ部計算
%混合率ゼロで接合部までテーパー
%等揚力を目指したい（願望）
disp("generating F wing");

%まずC*CLを台形化する
sum_ccl=sum(c_cl(t3:we_div));
e_ccl=sum_ccl*2/rows(c_cl(t3:we_div))-c_cl(t3);
%e_cclは翼端部のC*CL

%第四テーパ開始部の翼弦長を求める
% t4_root_c=zerofoil_def(c_cl(t3+1,1),mt_ang(t3,1)-i_ang_ed_ap(t3+1,1),comp_vel_ed_ap(t3+1,1),adata);
% f_data(t3+1,2)=t4_root_c;

%翼端部翼弦長を求める
e_c=zerofoil_def(e_ccl,mt_ang(t3,1)-i_ang_ed_ap(we_div,1),comp_vel_ed_ap(we_div,1),adata);

%翼端部翼翼弦長をもとに各翼素の弦長を定義
for i=t3+1:we_div
	t_ratio4=(e_c-f_data(t3,2))/(nang_cp(we_div,1)-nang_cp(t3,1));
	f_data(i,2)=t_ratio4*(nang_cp(i,1)-nang_cp(t3,1))+f_data(t3,2);
	
endfor

for t4=t3+1:we_div
	disp(t4);
	re(t4,1)=comp_vel_ed_ap(t4,1)*f_data(t4,2)/nu;
	f_data(t4,1)=nang_cp(t4,1);
	mt_ang(t4,1)=mt_ang(t3,1);
		%レイノルズ数、C、Cl、を固定しαを探索する
		%f_data(m1,3:8)=aoa_def(re(t4,1),cl(t4,1),f_data(t3,3),adata);
	f_data(t4,3)=0;
	f_data(t4,4)=mt_ang(t4,1)-i_ang_ed_ap(t4,1);
	f_data(t4,5)=Re_spline(re(t4),6,1,f_data(t4,4),2,adata);
	f_data(t4,6)=Re_spline(re(t4),6,1,f_data(t4,4),3,adata);
	f_data(t4,7)=Re_spline(re(t4),6,1,f_data(t4,4),5,adata);
	f_data(t4,8)=Re_spline(re(t4),6,1,f_data(t4,4),10,adata);
	mt_ang(t4,1)=f_data(t4,4)+i_ang_ed_ap(t4,1);
	
endfor

%ウイングレット部計算
%等揚力で翼端から３０ｃｍまで作る
disp("generating winglet");
	%t_ratio5=(c_cl(we_div+1,1)/f_data(we_div,5)/f_data(we_div,2))*0.;

	% f_data(t4+1,2)=c_cl(t4+1,1)/f_data(t4,5);
	% alpha_wl=(f_data(t4+1,2)-f_data(t4,2))/(nang_cp(t4+1,1)-nang_cp(t4,1))^0.01;

%30cmが何点目か
wlend_l=span/2+wl_span-0.3;
for wlend_p=we_div:we_div+wld_c
		if(nang_cp(wlend_p,1)>wlend_l)
			wlend_p=wlend_p-1;
			break;
		endif
endfor

%ウイングレット根弦長を定義
wl_root_c=zerofoil_def(c_cl(we_div+1,1),mt_ang(t4,1)-i_ang_ed_ap(we_div+1,1),comp_vel_ed_ap(we_div+1,1),adata);
f_data(we_div+1,2)=wl_root_c;

%まずC*CLを台形化する
sum_ccl_wl=sum(c_cl(we_div+1:wlend_p));
e_ccl_wl=sum_ccl_wl*2/rows(c_cl(we_div+1:wlend_p))-c_cl(we_div+1);
%e_cclは翼端部のC*CL

%翼端部翼弦長を求める
wl_e_c=zerofoil_def(e_ccl_wl,mt_ang(t4,1)-i_ang_ed_ap(wlend_p,1),comp_vel_ed_ap(wlend_p,1),adata);

%翼端部翼翼弦長をもとに各翼素の弦長を定義
for i=we_div+1:wlend_p
	t_ratio_wl=(wl_e_c-f_data(we_div+1,2))/(nang_cp(wlend_p,1)-nang_cp(we_div+1,1));
	f_data(i,2)=t_ratio_wl*(nang_cp(i,1)-nang_cp(we_div+1,1))+f_data(we_div+1,2);
	
endfor
	
for wl=we_div+1:wlend_p

	disp(wl);
		%f_data(wl,2)=f_data(wl-1,2)*t_ratio5;
		%f_data(wl,2)=alpha_wl*(nang_cp(wl,1)-nang_cp(t4,1))+f_data(t4,2);
	f_data(wl,1)=nang_cp(wl,1);
	re(wl,1)=comp_vel_ed_ap(wl,1)*f_data(wl,2)/nu;
	mt_ang(wl,1)=mt_ang(t3,1);
	f_data(wl,3)=0;
	f_data(wl,4)=mt_ang(wl,1)-i_ang_ed_ap(wl,1);
	f_data(wl,5)=Re_spline(re(wl),6,1,f_data(wl,4),2,adata);
	f_data(wl,6)=Re_spline(re(wl),6,1,f_data(wl,4),3,adata);
	f_data(wl,7)=Re_spline(re(wl),6,1,f_data(wl,4),5,adata);
	f_data(wl,8)=Re_spline(re(wl),6,1,f_data(wl,4),10,adata);
	
endfor

disp("generating winglet_end");
%ウイングレット翼端処理

disp("generating chord length");
for i=wl+1:we_div+wld_c
	f_data(i,2)=zerofoil_def(c_cl(i,1),mt_ang(wl,1)-i_ang_ed_ap(i,1),comp_vel_ed_ap(i,1),adata);
endfor

disp("define foil data");
for wle=wl+1:we_div+wld_c
	disp(wle);
	f_data(wle,1)=nang_cp(wle,1);
	% f_data(wle,2)=(t_ratio4*(wle-(t4+1))+f_data(t4+1,2))*sqrt(1-((wle-(t4+1))/((wld_c+we_div)-(t4+1)))^2);
	re(wle,1)=comp_vel_ed_ap(wle,1)*f_data(wle,2)/nu;
	mt_ang(wle,1)=mt_ang(t4,1);
	f_data(wle,3)=0;
	f_data(wle,4)=mt_ang(wle,1)-i_ang_ed_ap(wle,1);
	f_data(wle,5)=Re_spline(re(wle,1),6,1,f_data(wle,4),2,adata);
	cl(wle,1)=f_data(wle,5);
	f_data(wle,6)=Re_spline(re(wle,1),6,1,f_data(wle,4),3,adata);
	f_data(wle,7)=Re_spline(re(wle,1),6,1,f_data(wle,4),5,adata);
	f_data(wle,8)=Re_spline(re(wle,1),6,1,f_data(wle,4),10,adata);

endfor

%出力結果を元に最終的な循環を決定(誘導速度を解析していないため参考程度に)
re_gamma(:,1)=f_data(:,5).*f_data(:,2).*comp_vel_ed_ap(:,1)./2;

%翼面積を求める
for i=1:we_div+wld_c
	s_d(i,1)=f_data(i,2)*line_e_d(i,1);
endfor

mw_s=sum(s_d)*2;

figure(10)
	plot(f_data(1:we_div,1),f_data(1:we_div,2),f_data(we_div+1:we_div+wld_c,1),f_data(we_div+1:we_div+wld_c,2));
	ylabel("chord");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("re_chord.png",'-dpng','-r100')
	cd(now_work)

figure(11)
	plot(f_data(:,1),f_data(:,3));
	ylabel("mixture");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("re_mixture.png",'-dpng','-r100')
	cd(now_work)

figure(12)
	plot(f_data(:,1),mt_ang(:,1));
	ylabel("mt_ang");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("re_mt_ang.png",'-dpng','-r100')
	cd(now_work)

figure(13)
	plot(f_data(1:we_div,1),f_data(1:we_div,2),f_data(we_div+1:we_div+wld_c,1),f_data(we_div+1:we_div+wld_c,2));
	ylabel("chord");
	axis("equal");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("re_chord_equal.png",'-dpng','-r100')
	cd(now_work)

figure(14)
	plot(f_data(:,1),f_data(:,2),"-;chord;",f_data(:,1),f_data(:,3),"-;mixture;");
	xlabel('y[m]');
	grid on;
	cd("design_data")
	cd(pj_name)
	print("c_mix.png",'-dpng','-r100')
	cd(now_work)

figure(15)
	plot(f_data(:,1),ed_gamma(:,1),"-;ed_gamma;",f_data(:,1),re_gamma(:,1),"-;re_gamma;");
	xlabel('y[m]');
	grid on;
	cd("design_data")
	cd(pj_name)
	print("ed_re_gamma.png",'-dpng','-r100')
	cd(now_work)
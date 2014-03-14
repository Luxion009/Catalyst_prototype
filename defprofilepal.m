%翼型作成を補助する

%cl*cを定義
for i=1:(wld_c+we_div)
	c_cl(i,1)=2*ed_Hex_gamma(i,1)/ed_Hex_comp_vel(i,1);
endfor


do

	%主翼限界レイノルズ数定義
	lim_mre=input("主翼限界レイノルズ数=");

	%限界レイノルズ数から最外翼弦長を定義
	chord(we_div,1)=lim_mre*nu/d_vel;
	disp("最外翼端弦長");
	disp(chord(we_div,1));

	%中央翼取り付け角
	ang_cw=input("中央翼有効迎角＝");
	mt_ang=ang_cw+ed_Hex_i_ang(1,1)

	%中央翼よく弦長定義
	chord(1:dw_divn,1)=input("中央翼翼弦長(m)=");

	%ない翼たん
	chord(dw_divn*2,1)=input("内翼翼弦長(m)=");

	%外翼端
	chord(dw_divn*3,1)=input("外翼端翼弦長(m)=");

	%順最外翼端
	chord(dw_divn*4,1)=input("準最外翼端弦長(m)=");
	
	disp("翼弦長傾斜");
	for i=2:w_div_c
		i
		chord_alpha(i,1)=(chord(dw_divn*i,1)-chord(dw_divn*(i-1),1))/divw_span(i,1);
		disp(chord_alpha(i,1));

	endfor

	%翼弦長分布を作成
	disp("generating chord");
	for i=2:w_div_c
		chord_d(i,1)=(chord(dw_divn*(i),1)-chord(dw_divn*(i-1),1))/dw_divn;
		chord((i-1)*dw_divn+1,1)=chord((i-1)*dw_divn,1)+chord_d(i,1);
		for j=2:dw_divn
			chord((i-1)*dw_divn+j,1)=chord((i-1)*dw_divn+j-1,1)+chord_d(i,1);
		endfor
	endfor

	%レイノルズ数分布作成
	disp("generating re");
	chord_re(:,1)=chord(:,1).*d_vel/nu;

	%Cl分布作成
	disp("generating cl");
	chord_cl(1:we_div,1)=c_cl(1:we_div,1)./chord(1:we_div,1);


	figure(10);
	plot(nang_line(1:we_div),chord(1:we_div),"-");

	for i=1:w_div_c
		printf("%d翼\n",i)
		disp("Re数");
		disp(chord_re(dw_divn*i,1));
		disp("cl");
		disp(chord_cl(dw_divn*i,1));
		disp("誘導迎え角")
		disp(ed_Hex_i_ang(dw_divn*i,1));
		disp("有効迎角")
		disp(mt_ang-ed_Hex_i_ang(dw_divn*i,1))
		disp("---------------------")
	endfor

	chord_sig=input("1:ok 0:return ---");

until(chord_sig==1)

%各翼端必要データを算出
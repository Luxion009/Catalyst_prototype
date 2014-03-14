%つくったひと
%WASA14年設計主任　なぽー
%twitter:Luxion009

printf("\n");
disp("人力飛行機統合設計プログラム")
disp("---Catalyst prototype---");
printf("\n");
disp("");

data_load

do
	printf("[1]空力解析\n[2]翼型作成\n[3]平面型作成\n[4]空力解析\n[5]データ保管\n[6]桁設計\n[6.5]桁再定義\n[7]慣性モーメント定義\n[7.5]安定微係数計算\n[8]運動解析\n[9]データ出力\n[10]終了\n");
	mode=input("モードセレクト ---");

	if(mode==1)
		%循環を求める
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

		%翼型作成を補助する
		defprofilepal

		savesig=0;
		savesig=input("save data ok:1 cancel:0 ---");
		if savesig==1
			data_save
		endif

		% %作成した中央翼データから各部分よく端翼型を生成する
		% defprofile]

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=3;
		endif

	endif

	if(mode==3)

		%平面形を作る
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
		%空力解析
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
		% %座標点データ取得
		% foil_p_input
		% %翼厚定義
		% defthn
		%マンドレル設計
		defmand
		%データ補間
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
		
		%材料計算
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

		%部分的材料計算
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

		%慣性モーメント定義
		definertia

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=7.5;
		endif
	endif

	if(mode==7.5)
		%運動解析に必要な微係数の推算
		%上半角効果
		defclbeta
		%ロール角速度
		def_p
		%ヨー角速度
		def_r

		con=input("continue ok:1 end:0 ---");
		if(con==1)
			mode=8;
		endif
	endif


	if(mode==8)
		
		%テール定義
		deftail

		%水平尾翼定義
		deftwing

		%垂直尾翼定義
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
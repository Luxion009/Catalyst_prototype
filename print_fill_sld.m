%ソリッド用データ保管
clear ribposdata line_e_rib line_e_rib_d

%リブ間配列（中央ーリブ１が要素１）
%中央翼
ribposdata=[0.214666/2];
for i=2:8
	ribposdata(i,1)=0.214666;
endfor
for i=9:49
	ribposdata(i,1)=0.23;
endfor
for i=50:58
	ribposdata(i,1)=0.2;
endfor

hspan_rib=12.84;

max_tawami_rib=max(S_mat_swang(:,1));

%たわみ量は楕円翼から固定
alpha_rib=max_tawami_rib/(hspan_rib*hspan_rib);
a=0;
f_v=0;

ribposdata(:,1)=ribposdata(:,1)./2;

%線素の端点のy座標を求める（原点なし）
for i=1:rows(ribposdata(:,1))

	line_e_rib(i,1)=newton_met(a,alpha_rib,ribposdata(i,1),f_v);
	a=line_e_rib(i,1);
	f_v=a;

endfor

line_e_rib(:,2)=alpha_rib.*((line_e_rib(:,1)).^2);

rib_dif(:,1)=2.*alpha_rib.*line_e_rib(:,1);

rib_theta(:,1)=atan(rib_dif(:,1));

% line_e_rib_d(1,1)=sqrt((line_e_rib(1,1))^2+(line_e_rib(1,2))^2);

% for d_count=2:rows(ribposdata(:,1))
% 	line_e_rib_d(d_count,1)=sqrt((line_e_rib(d_count,1)-line_e_rib(d_count-1,1))^2+(line_e_rib(d_count,2)-line_e_rib(d_count-1,2))^2);
% endfor

% figure(111)
% 	plot(line_e_rib(:,1),line_e_rib(:,2),"-o",[0 12.84],[0 0]);

%翼型再読み込み
fp=fopen('use_foila.txt','r');
use_foil_ps{1}=fgetl(fp);
fclose(fp);

fp=fopen('use_foilb.txt','r');
use_foil_ps{2}=fgetl(fp);
fclose(fp);

fp=fopen('use_foilc.txt','r');
use_foil_ps{3}=fgetl(fp);
fclose(fp);

fp=fopen('use_foild.txt','r');
use_foil_ps{4}=fgetl(fp);
fclose(fp);

disp("input profiledata_ps");

for i=1:4
	%-----翼型を読込
	i;
	fp=fopen(strcat(use_foil_ps{i},'.dat'));
	buff=fgetl(fp);
	foil_ps{i}=(fscanf(fp,'%f',[2,350]))';
	foil_ps{i}=(foil_ps{i});
	fclose(fp);
	fs(i)=0;

	%-----翼の前縁を探索
	do	
		fs(i)=fs(i)+1;
	until(foil_ps{i}(fs(i),1)<=min(foil_ps{i}(:,1)))

end
	
%-----翼型データを整理
%-----x方向分割数
n_xdiv_l=80;
n_xdiv_t=60;

disp("generating xy_foil_ps");

x_foil_ps(1:n_xdiv_t+n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([n_xdiv_t+n_xdiv_l-1:-1:1]').^2;	
x_foil_ps(n_xdiv_t+n_xdiv_l,:)=0;
x_foil_ps(n_xdiv_t+n_xdiv_l+1:2*n_xdiv_t+2*n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([1:n_xdiv_t+n_xdiv_l-1]').^2;

for i=1:4
	y_foil_ps(1:n_xdiv_l+n_xdiv_t,i)=interp1(foil_ps{i}(1:fs(i),1),foil_ps{i}(1:fs(i),2),x_foil_ps(1:n_xdiv_l+n_xdiv_t,1),"extrap");
	y_foil_ps(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,i)=interp1(foil_ps{i}(fs(i):rows(foil_ps{i}),1),foil_ps{i}(fs(i):rows(foil_ps{i}),2),x_foil_ps(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,1),"extrap");
	y_foil_ps(n_xdiv_l+n_xdiv_t,i)=0;
end


end_pp=print_sld(0,0,dwing_chord(1,1),dwing_aoa(1,1),0,dwing_spos(1,1),[x_foil_ps(:,1) y_foil_ps(:,1)],'C',1,pj_name,now_work);

%A翼
for i=2:9
	end_pp=print_sld(line_e_rib(i-1,1),line_e_rib(i-1,2),dwing_chord(1,1),dwing_aoa(1,1),rib_theta(i-1,1),dwing_spos(1,1),[x_foil_ps(:,1) y_foil_ps(:,1)],'C',i,pj_name,now_work);
endfor

%B翼
disp("B wing")
for i=2:15
	i
	buf_yfoil=y_foil_ps(:,1).*print_mix(2,i).+y_foil_ps(:,2).*(1-print_mix(2,i));
	end_pp=print_sld(line_e_rib(i+7,1),line_e_rib(i+7,2),print_chord(2,i),print_aoa(2,i),rib_theta(i+7,1),print_spos(2,i),[x_foil_ps(:,1) buf_yfoil(:,1)],'N',i,pj_name,now_work);
	clear buf_yfoil
endfor

%C翼
disp("C wing")
for i=2:15
	i
	buf_yfoil=y_foil_ps(:,2).*print_mix(3,i).+y_foil_ps(:,3).*(1-print_mix(3,i));
	end_pp=print_sld(line_e_rib(i+21,1),line_e_rib(i+21,2),print_chord(3,i),print_aoa(3,i),rib_theta(i+21,1),print_spos(3,i),[x_foil_ps(:,1) buf_yfoil(:,1)],'G',i,pj_name,now_work);
	clear buf_yfoil
endfor

%D翼
disp("D wing")
for i=2:14
	i
	buf_yfoil=y_foil_ps(:,3).*print_mix(4,i).+y_foil_ps(:,4).*(1-print_mix(4,i));
	end_pp=print_sld(line_e_rib(i+35,1),line_e_rib(i+35,2),print_chord(4,i),print_aoa(4,i),rib_theta(i+35,1),print_spos(4,i),[x_foil_ps(:,1) buf_yfoil(:,1)],'FSG',i,pj_name,now_work);
	clear buf_yfoil
endfor

disp("E wing")
for i=2:10
	i
	buf_yfoil=y_foil_ps(:,4);
	end_pp=print_sld(line_e_rib(i+48,1),line_e_rib(i+48,2),print_chord(5,i),print_aoa(5,i),rib_theta(i+48,1),print_spos(5,i),[x_foil_ps(:,1) buf_yfoil(:,1)],'SSG',i,pj_name,now_work);
	clear buf_yfoil
endfor

%—ƒŒúŒvZŠÖ”

function res_thn=calc_thn(foila,foilb,at,bt,mix,pos)
	aa=spline(foila(1:at,1),foila(1:at,2),pos);
	ad=spline(foila(at+4:rows(foila),1),foila(at+4:rows(foila),2),pos);
	ba=spline(foilb(1:bt,1),foilb(1:bt,2),pos);
	bd=spline(foilb(bt+4:rows(foilb),1),foilb(bt+4:rows(foilb),2),pos);

	fstthn=abs(aa-ad);
	secthn=abs(ba-bd);

	res_thn=fstthn*mix+secthn*(1-mix);
endfunction
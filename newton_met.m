function n_res=newton_met(a,alpha,d_s,f_v)
	max=5000;
	eps=0.000001;
	count=0;
	do
		n_v=f_v-1*((alpha^2)*(f_v^4)+(1-2*(alpha^2)*(a^2))*(f_v^2)-2*a*f_v+a^2+(alpha^2)*(a^4)-(2*d_s)^2)/(4*(alpha^2)*(f_v^3)+2*(1-2*(alpha^2)*(a^2)*f_v));
		if(abs(n_v-f_v)<eps)
			break
		endif
		f_v=n_v;
		count=count+1;
	until(count>max)
	
	n_res=n_v;
endfunction
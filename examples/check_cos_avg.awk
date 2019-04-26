awk '
# Return gaussian random number
function grand(mu, sigma){
		
	do {
		u = 2.0 * rand() - 1.0
		v = 2.0 * rand() - 1.0
		r = u * u + v * v
			
	} while ( r > 1 )
		
	fc = sqrt( -2.0 * log(r) / r )
	z0 = u * fc

	return z0 * sigma + mu
}
# Minimum
function min(a){
	
	if ( rand() < a ) return 1
	else return 0

}
# Pick value with a given proportion 
function pick(a, b, c){
	
	s = min(1.0)
	d = min(0.5)
		
	g = int(d) * b + (1 - int(d)) * c 
	f = int(s) * g + (1 - int(s)) * a
	
	return f	

}
# Main function
function main(num1, num2){
	
	PI = atan2(0, -1)

	srand(123)
	
	for (i=0; i<=300000; i++) {
		
		mu1 = 0.0
		mu2 = PI/3.0
		mu3 = -PI/3.0
	
		a1 = grand(mu1, PI/12.0)
		a2 = grand(mu2, PI/12.0)
		a3 = grand(mu3, PI/12.0)

		X = pick(a1, a2, a3)
		print 180 * X / PI

		cpsi  =  cos(X)
		spsi  =  sin(X)

		avg2 +=  cpsi*cpsi
		acps +=  cpsi
		asps +=  spsi

		n++
		
	}

	psi = atan2(asps/n, acps/n)
	
	p2_1 = 0.5 * ( 3 * (avg2/n) - 1)
	p2_2 = 0.5 * ( 3 * cos(psi) - 1)

	print "# average cosine "p2_1
	print "# average angles "p2_2
}
# Script execution starts here
BEGIN{
	main()
}'

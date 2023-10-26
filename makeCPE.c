// Program to create a SPICE subcircuit using the Morrison-Scott-Seshadri method
// JBS January 2023, with Marcus, Vance & Logan

#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <math.h>
#include    <errno.h>
#include    <complex.h>


//#define MIN(A,B) (((A)<(B))?(A):(B))
#define PI 3.141592653589793
			
//3.141592654
#define TWOPI 6.283185307179586
//6.283185308

void err(char *error_text)
{
	fprintf(stderr,"%s\n",error_text);
	exit(1);
	}

#define MAXIMUM(A,B) (((A)>(B))?(A):(B))
static char *engstr(double x, int digits)	// return char* in engineering format
{
    int dp;
    int power=0;
    static short strptr=0;
    static char str[9][32];

    strptr = (strptr + 1) % 9;
    if(x==0) {
        sprintf(str[strptr], " 0.");
		for(dp=0;dp<digits-1;dp++) str[strptr][dp+3]='0';
		str[strptr][dp+3]='e';
		str[strptr][dp+4]='0';
		str[strptr][dp+5]='\0';
        return(str[strptr]);
        }
    while ( fabs(x)<1 ) {
        x = x *1000;
        power+=(-3);
        }
    while ( fabs(x)>1000 ) {
        x = x /1000.0;
        power+=3;
        }
    dp = digits -1 -(fabs(x)>=10) -(fabs(x)>=100) ;
    dp = MAXIMUM(0,dp); // extend if asked for <3 and x is big
    sprintf(str[strptr], "%.*lfe%+d", dp, x, power);
	// for const width: "%*.*lfe%+d", 3+digits, dp, x, power);
    return(str[strptr]);
    }
    
static char *sengstr(double x, int digits)	// similar to engstr, but letters for multiplier, a la SPICE
{
    int dp;
    int power=0;
	char c;
    static short strptr=0;
    static char str[9][32];

    strptr = (strptr + 1) % 9;
    if(x==0) {
        sprintf(str[strptr], " 0.");
		for(dp=0;dp<digits-1;dp++) str[strptr][dp+3]='0';
		str[strptr][dp+3]='\0';
        return(str[strptr]);
        }
    while ( fabs(x)<1 ) {
        x = x *1000;
        power+=(-3);
        }
    while ( fabs(x)>1000 ) {
        x = x /1000.0;
        power+=3;
        }
    dp = digits -1 -(fabs(x)>=10) -(fabs(x)>=100) ;
    dp = MAXIMUM(0,dp); // extend if asked for <3 and x is big
	switch(power) {
		case -18: c='a'; break;
		case -15: c='f'; break;
		case -12: c='p'; break;
		case -9: c='n'; break;
		case -6: c='u'; break;
		case -3: c='m'; break;
		case 0: c=' '; break;
		case 3: c='k'; break;
		case 6: c='M'; break;
		case 9: c='G'; break;
		case 12: c='T'; break;
		case 15: c='P'; break;
		default: c='?';
	}
    sprintf(str[strptr], "%.*lf%c", dp, x, c);
	// for const width: "%*.*lf%c", 3+digits, dp, x, c);
    return(str[strptr]);
    }
#undef MAXIMUM



int main(int argc, char* argv[])
{

	//double magY;			// the target impedance magnitude at f0
	double f0;//home frequency 
	double f_min, f_max;//frequency limits
	double k; //resistor spacing
	double k_f;//frequency spacing
    char *CPEName; //spice subcircuit name
    char *fileName; //output file, should end in .inc
	double m;
	double ytheta;			// from Morrison (49?)
	double thisC, thisR, thisf;
	double C0, R0;//home branch values
	double Ctotal, Rtotal;
	int i;
	int nh, nl;//number of branches
	double C_T, R_T;//terminations
	int narg=0, intnode;
	FILE *inc;
	double alpha, Z0;
	double km1; //morrison k^(m-1) capacitor spacing
	double C_f;						// the CPE mag constant (pseudo-capacitance)
	double k_f_perDecade;//shows if accuracy is bandwidth or branch limited, from testing, a value of 6.1 appears to give best accuracy
	double suggested_k_f;//aproximatly the maximum k_f that does not limit accuracy
    //complex double Ytheta,Ztheta;
    complex double ctmp;
	double minC, maxC, minR, maxR;

	// version 1.00: replacing/fixing SS python version
	// tested by comparison with Farrow version, separately coded
	// 1.02: Added comments 
	// 1.03: allow C_f parameter, if # params is 5; fix units of w0 to be radians
	// 1.04: include terminating components
	// 1.05: print largest and smallest R & C
	// 1.06: fix bug with radians and hertz specifying Cf
	// 2.00: all new approach with help from Marcus, checked against Vance & Logan
	// 2.01: precision improvements, now 17 digit to max out the accuracy of double precision
	// 2.02: minor debug information improvements
    // 2.02: added CPEName and fileName to support simulating multiple CPEs made by the same program in the same circuit
    float version = 2.03;    
    if (argc<5+1 || argc>6+3) { // ??
        fprintf(stderr,"makeCPE V%.2f jbs Jan 2023\n", version);
        fprintf(stderr,"Creates a SPICE subckt to simulate a CPE of given parameters.\n");
        fprintf(stderr,"Usage: makeCPE alpha |Z_f0| f0 f_min f_max k_f fileName CPEname\n");
        fprintf(stderr,"...or: makeCPE alpha C_f f_min f_max k_f fileName CPEname\n");
        fprintf(stderr,"where-  alpha is the order of the derivative in I = C_F d^(alpha)V/dt^(alpha);\n");
		fprintf(stderr,"        |Z_f0| is the magnitude of CPE impedance at f0 (in Hertz);\n");
		fprintf(stderr," or     C_f is the CPE constant from the fit;\n");
        fprintf(stderr,"        f_min and f_max are the frequency limits of the approximation;\n");
        fprintf(stderr,"        k>1 is the accuracy parameters, more accurate for k closer to 1.\n");
        fprintf(stderr,"        CPEname is the Spice subcircuit name"); 
		fprintf(stderr,"        fileName is the output file name"); 
        fprintf(stderr,"Writes file <fileName> containing the SPICE subcircuit code \n");
        fprintf(stderr,"Value of alpha must be 0<alpha<1, and f_min<f0<f_max.\n");
        fprintf(stderr,"Set alpha=0.5 for a Warburg element.\n");
        fprintf(stderr,"It is recommended to choose frequency limits >10x beyond use range.\n");
        fprintf(stderr,"fileName and CPEname default to cpe1.inc and CPE_1 if neither are provided.\n");
        fprintf(stderr,"\n");
        exit(1);
    }
    
    // process input arguments
    alpha = atof(argv[++narg]);
	if(alpha<0.0001) err("alpha is too small");//untested values, must be greater than 0
	if(alpha>0.9999) err("alpha is too large");//untested values, must be less than 1

	if(alpha<0.01) 	fprintf(stderr,"Warning, low alpha value, to avoid excessive rounding error, place the first node of the CPE as close to 0V as possible\n");
	if(alpha>0.99) 	fprintf(stderr,"Warning, high alpha value, to avoid excessive rounding error, place the second node of the CPE as close to 0V as possible\n");
	//noise floor = precision* (branches^alpha) * (frequency/min_branch_freq)^(1-alpha)
	//reverse oreintation noise floor = precision* (branches^(1-alpha)) * (frequency/max_branch_freq)^(alpha)

	if(((argc==8)||(argc==6))){	// user specified alpha & Cf
		C_f = atof(argv[++narg]);
		if(C_f<=0.00) err("C_f must be >0.");
	}else{			// user specified mag & freq_radians
	    Z0 = atof(argv[++narg]);
		if(Z0<=0.00) err("Z0 must be >0.");
	    f0 = atof(argv[++narg]);
		if(f0<=0.00) err("f0 must be positive.");
	}

    f_min = atof(argv[++narg]);
	if(((argc==9)||(argc==7)) && f_min>=f0) err("f_min must be lower than f0");

    f_max = atof(argv[++narg]);
	if(f_max<=f_min) err("f_min >= f_max");
	if(((argc==9)||(argc==7)) && f_max<=f0) err("f_max must exceed f0.");

    k_f = atof(argv[++narg]);//frequency based spacing of branches
	if(k_f<1.00000001) err("k must exceed 1.0000000.");
	k_f_perDecade=(k_f-1)* log10(f_max/f_min); //the highest accuracy k_f value usually gives a value of around 6.1
	suggested_k_f=1+(7/log10(f_max/f_min));
	if(k_f_perDecade>7) fprintf(stderr,"unnecessarily low k value, consider increasing k or increasing the bandwidth\nrecommended k value is %s\n",engstr(suggested_k_f,4));
	k = pow(k_f,alpha);	// a la Morrison



	fprintf(stderr,"-- makeCPE v%.2f \n",version);
   if((argc==9)||(argc==8)){
        fileName = argv[++narg];
    CPEName = argv[++narg];
   }
   else{
    fileName = "cpe1.inc";
    CPEName = "CPE_1";
   }
    
	// now open output file
	inc = fopen(fileName,"w+");						// inc file open 
	if(inc==NULL) err("Cannot open cpe.inc");

	// echo information into the include file as comments
	fprintf(inc,"* CPE alpha = %s\n", engstr(alpha,4) );
	if(((argc==8)||(argc==6))){	// user specified alpha & Cf
		fprintf(inc,"* given C_f = %s\n", engstr(C_f,4) );
		f0 = sqrt(f_min*f_max);	
	}else{			// user gave Z at freq
		fprintf(inc,"* given CPE |Z|   = %s Ohms\n", engstr(Z0,4) );
		fprintf(inc,"*        at f0    = %sHz or %s rad/s\n", sengstr(f0,4), engstr(f0*TWOPI,4) );
	}
	fprintf(inc,"* CPE f_min/f_max = %s/%s\n", sengstr(f_min,4), sengstr(f_max,4) );
	fprintf(inc,"* CPE k = %.6f\n", k );
	fprintf(inc,"* CPE k_f = %s, recommended k_f value = %s\n",engstr(k_f,6),engstr(suggested_k_f,6));

	// print the opposite parameter set from that given, as a check	
	if((argc==9)||(argc==7)){	// given Z@w so must find C_f
		// know that Z_cpe = 1/(C_f * (j w)^alpha), so C_f=abs( 1/(Z_cpe * (jw)^alpha ) )
		ctmp = CMPLX(0.00,TWOPI*f0);
		C_f = cabs( 1.0/(Z0 * cpow(ctmp,alpha) ) );
		fprintf(inc,"* check CPE C_f = %.6f\n", C_f );
	}else{			// given C_f so echo mag at w0
		ctmp = CMPLX(0.00,TWOPI*f0);
		Z0 = cabs( 1.0/(C_f * cpow(ctmp,alpha) ) );
		fprintf(inc,"* check Z0 = %.6f at f0 = %sHz or %s rad/s\n", Z0, sengstr(f0,4), engstr(f0*TWOPI,4) );
	}

	// fprintf(inc,"* noise floor if node cpe2 is ground = precision* (branches^alpha) * (frequency/min_branch_freq)^(1-alpha)\n");
	// fprintf(inc,"* approximately = precision* (branches^%s) * (frequency/%s)^(%s)\n",engstr(alpha,4),engstr(f_min,4),engstr((1-alpha),4));
	// fprintf(inc,"* noise floor if node cpe1 is ground = precision* (branches^(1-alpha)) * (frequency/max_branch_freq)^(alpha)\n");
	// fprintf(inc,"* approximately = precision* (branches^%s) * (frequency/%s)^(%s)\n",engstr((1-alpha),4),engstr(f_max,4),engstr(alpha,4));

	m = 1.00/alpha;
	km1 = pow(k,m-1);
	
	// calculate number of branches
	nh = floor(log10(f_max/f0)/log10(k_f));
	nl = floor(log10(f0/f_min)/log10(k_f));
	fprintf(stderr,"Subcircuit has 1+%d+%d+2 = %d branches including terminating branches\n",nl, nh, nl+nh+3);

	
	// now find home branch at f0
	ytheta = PI / ( cos((1.0-2.0/m)*PI/2.0) * (m * log(k)) );
	R0 = Z0 * ytheta;
	C0 = 1.0/(TWOPI*R0*f0);
	
	// now we can generate the SPICE subcircuit, except for terminations
	fprintf(inc,".subckt %s cpe1 cpe2\n",CPEName);
    
	// now print out Morrison branches, sum up some stats
	intnode=1;	// node number for the branch 
	// write home branch first
	fprintf(inc,"* home branch\n");
	fprintf(inc,"C%d cpe2 n%d %s\n", intnode, intnode, engstr(C0,17) );
	fprintf(inc,"R%d cpe1 n%d %s\n", intnode, intnode, engstr(R0,17) );
	
	Rtotal = R0;
	Ctotal = C0;
	// summation in high frequency direction
	for(i=1;i<=nh;i++){
		intnode+=1;
		thisR = R0 / pow(k,i);
		thisC = C0 / pow(km1,i);
		thisf = 1.0/(TWOPI*thisR*thisC);
		fprintf(inc,"* branch %d, centre frequency = %s\n", intnode, engstr(thisf,6) );
		fprintf(inc,"C%d cpe1 n%d %s\n", intnode, intnode, engstr(thisC,17) );
		fprintf(inc,"R%d cpe2 n%d %s\n", intnode, intnode, engstr(thisR,17) );
		Rtotal = 1.0/(1.0/Rtotal + 1.0/thisR);
		Ctotal += thisC;
	}
	minC = thisC;
	minR = thisR;
	// the HF termination
	C_T = thisC/(km1-1);
	Ctotal += C_T;
	fprintf(inc,"CtermHF cpe1 cpe2 %s\n", engstr(C_T,17) );
	// summation in low direction
	for(i=1;i<=nl;i++){
		intnode+=1;
		thisR = R0 * pow(k,i);
		thisC = C0 * pow(km1,i);
		thisf = 1.0/(TWOPI*thisR*thisC);
		fprintf(inc,"* branch %d, centre frequency = %s\n", intnode, engstr(thisf,6) );
		fprintf(inc,"C%d cpe1 n%d %s\n", intnode, intnode, engstr(thisC,17) );
		fprintf(inc,"R%d cpe2 n%d %s\n", intnode, intnode, engstr(thisR,17) );
		Ctotal += thisC;
		Rtotal = 1.0/(1.0/Rtotal + 1.0/thisR);
	}
	maxC = thisC;
	maxR = thisR;
	// the LF termination
	R_T = thisR * (k-1);
	Rtotal = 1.0/(1.0/Rtotal + 1.0/R_T);
	fprintf(inc,"RtermLF cpe1 cpe2 %s\n", engstr(R_T,17) );

	fprintf(inc,".ends\n");
	
	// values should make sense
	fprintf(stderr,"min C = %s, max C = %s, ratio = %s\n", engstr(minC,3), engstr(maxC,3), engstr(maxC/minC,3) );
	fprintf(stderr,"min R = %s, max R = %s, ratio = %s\n", engstr(minR,3), engstr(maxR,3), engstr(maxR/minR,3) );
	fprintf(stderr,"#branches=%d, total Cap = %s, total R// = %s\n", nh+nl+3, 
		sengstr(Ctotal,4), sengstr(Rtotal,4));

	fprintf(stderr,"C0= %s, R0=%s\n",engstr(C0,17),engstr(R0,17));
	fprintf(stderr,"CTerm= %s, RTerm=%s\n",engstr(C_T,17),engstr(R_T,17));
	
	// fprintf(stderr,"noise floor = precision* (branches^alpha) * (frequency/min_branch_freq)^(1-alpha)\n");
	// fprintf(stderr,"noise floor = %g\n", precision* (branches^alpha) * (frequency/min_branch_freq)^(1-alpha)\n);
	// fprintf(stderr,"reverse oreintation noise floor = precision* (branches^(1-alpha)) * (frequency/max_branch_freq)^(alpha)\n");
	fclose(inc);
}


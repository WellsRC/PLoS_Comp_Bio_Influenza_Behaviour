/*
 *  Vacc_SS.c
 *  
 *
 *  Created by Chad Wells on 1/10/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define PI 3.14159265

float urand()   /* returns random number between 0 and 1 */
{
	float answer;
	answer=(float)rand()/RAND_MAX;
	return(answer);
}
int Vacday(int day,float phi0,float phi1,float phi2,float phi3)
{
	float ftemp=urand();
	int Vday=-1;
	if(day<30)
	{
		if(ftemp<(phi0/(phi0+phi1+phi2+phi3))) //vaccination in Sept
		{
			Vday=rand()%(30-day) + day; //randomly chooses a day for the individual to vaccinate 
		}
		else if (ftemp<((phi1+phi0)/(phi0+phi1+phi2+phi3))) //Vaccination in Oct
		{
			Vday=rand()%31+30; //randomly chooses a day for the individual to vaccinate 
		}
		else if (ftemp<((phi1+phi0+phi2)/(phi0+phi1+phi2+phi3))) //Vaccination in Nov
		{
			Vday=rand()%30+61; //randomly chooses a day for the individual to vaccinate 
		}
		else  //Vaccination in Dec
		{
			Vday=rand()%31+91; //randomly chooses a day for the individual to vaccinate 
		}
	}
	else if(day<61)
	{
		if (ftemp<((phi1)/(phi1+phi2+phi3))) //Vaccination in Oct
		{
			Vday=rand()%(31-(day-30))+day; //randomly chooses a day for the individual to vaccinate 
		}
		else if (ftemp<((phi1+phi2)/(phi1+phi2+phi3))) //Vaccination in Nov
		{
			Vday=rand()%30+61; //randomly chooses a day for the individual to vaccinate 
		}
		else //Vaccination in Dec
		{
			Vday=rand()%31+91; //randomly chooses a day for the individual to vaccinate 
		}
	}
	else if(day<91)
	{
		if (ftemp<((phi2)/(phi2+phi3))) //Vaccination in Oct
		{
			Vday=rand()%(30-(day-61))+day; //randomly chooses a day for the individual to vaccinate 
		}
		else //Vaccination in Dec
		{
			Vday=rand()%31+91; //randomly chooses a day for the individual to vaccinate 
		}
	}
	else
	{
		Vday=rand()%(31-(day-91))+day; //randomly chooses a day for the individual to vaccinate
	}
	return(Vday);
}
float lognorm(float ev,float stdev, float dx)
{
	float answer;
	//parameters for lognormal distribtion with given expected value (ev) and variance (var)
	float mu=log(ev)-log(1+(pow(stdev,2)/(pow(ev,2))))/2; 
	float sigmasq=log(1+pow(stdev,2)/(pow(ev,2)));
	//to determine the number
	int jj=1;
	float temp=urand();
	float logn=exp(-pow(log(dx*jj)-mu,2)/(2*sigmasq))/(dx*jj*sqrt(2*sigmasq*PI))*dx;
	while(temp>logn)
	{
		jj++;
		logn=logn+exp(-pow(log(dx*jj)-mu,2)/(2*sigmasq))/(dx*jj*sqrt(2*sigmasq*PI))*dx;
	}
	if(dx*jj>1)
	{
		answer=1;
	}
	else
	{
		answer=dx*jj;
	}
	return(answer);	
}
long double factorial(long double num)
{
	long double answer;
	if(num<=1)
	{
		answer=1;
	}
	else
	{
		answer=num*factorial(num-1);
	}
	return(answer);
}
int SS(double Ro95, long double k, long double p, int delta)
{
	long double sum;
	sum=k*(1-pow((1-p),(double)delta));
	if (sum>Ro95)
	{
		return(2);
	}
	else
	{
		return(0);
	}
}
typedef struct
{
	int state; // what the current state of the individual is in
	/* The possible states
	 -1 -Vaccinated but suscpetible
	 0- Susceptible
	 1- Infected
	 2- Recovered
	 3- Vaccinated */
	int Q; //number of contacts
	int Qcn; // number of contacts (used in creating network)
	int super; // relation or not to super-spreader
	/* Possible relation
	 0- no relation to a super-spreader
	 1- contact of super-spreader
	 2- super-spreader 
	 3- super spreader and contact of super spreader*/
	int ctv; // if the individual has been choosen to vaccinate this season 
	/*
	 0- individual has been chosen for vaccination
	 1- chasen and vaccinated
	 2- not choosen for vaccination
	 3- has been infected and will not vaccinate
	 */
	float peps;// perceived vaccine efficacy
	float MI;// monetary incentive
	int TItemp; // used in the imiation 
	int TCtemp; //used in the imitation
	float pepstemp;//used in the imiation 
	int TI; //time since the individuals last infection
	int TC; //time sicne the individuals last vaccine complication
	float p; //the individual's infectability
	float avgp; //the individuals average infectiousness
	float sp; // the individuals variability in infectiousness
	float avgpeps; // average perceived vaccine efficacy
	double pn; // proabaility of becoming infected
	int delta; // the individual's infectious period
	int dc;// count the infectious period
	int Vday; // the day in whic hthe individual will vaccinate
	int VacLast; // determine whether or not vaccinated last season or not
	int* contact; // the nodes contacts} node_descript;
}node_descript;
int main()
{
	int inc=0; // determine whether incentives are used or not
	// inc =0 - no incetives
	//inc=1 - incentives used
	int strat=0; // the vaccination strategy being used 
	/*
	 0-passive
	 1- random 
	 2- nearest neighbour
	 3- page rank
	 4 - improved nearest neighbour
	 */		
	int htran=0; //heterogenity in transmission 0- No 1- Yes
	int hip=0; //heterogenity in infectious period 0- No 1- Yes	
	FILE *f2;
	f2=fopen("Parameters-Poisson.txt","w");
	FILE *f3I;
	f3I=fopen("Infection-daily-Poisson.txt","w");
	FILE *f3V;
	f3V=fopen("Vaccination-daily-Poisson.txt","w");
	FILE *f4I;
	f4I=fopen("TotalInfection-Season-Poisson.txt","w");
	FILE *f4V;
	f4V=fopen("TotalVaccination-Season-Poisson.txt","w");
	FILE *f6I;
	f6I=fopen("TotalSSInfection-Season-Poisson.txt","w");
	FILE *f6V;
	f6V=fopen("TotalSSVaccination-Season-Poisson.txt","w");
	FILE *f8V;
	f8V=fopen("SS-SelfVac-Poisson.txt","w");
	FILE *f9;
	f9=fopen("Degree-Infection-Poisson.txt","w");
	FILE *f10;
	f10=fopen("Degree-Vaccination-Poisson.txt","w");
	FILE *f11;
	f11=fopen("NumSS-NumSSCont-Poisson.txt","w");
	FILE *f12;
	f12=fopen("NumberDegree-Poisson.txt","w");
	FILE *f13;
	f13=fopen("SS-SelectedRandom-Poisson.txt","w");
	FILE *f14;
	f14=fopen("SS-NNSelected-Poisson.txt","w");
	FILE *f15;
	f15=fopen("Degree-avgpeps-Poisson.txt","w");
	FILE *f16;
	f16=fopen("SSChar-Poisson.txt","w");
	FILE *f17;
	f17=fopen("Incentives-Poisson.txt","w");
	srand(time(NULL));
	int N=10000; //population size
	node_descript ind[N]; //creation of the individuals
	int TransYear=125; // years to run out transient
	fprintf(f2,"TransYear %d \n",TransYear);
	int Nyears=25; // number of years to run for
	fprintf(f2,"NYears %d \n",Nyears);
	int NSims=400; //number of simulations to run
	fprintf(f2,"NSims %d \n",NSims);
	float q=0.44; // for Ro and the calc of avgp
	fprintf(f2,"q %f \n",q);
	float avgp=q/39; //average probability of infection along one edge per day corresponding to the average degree of network random network	
	float sp=0.03; //seasonal change in the transmission rate corresponding to the uniform random network	
	fprintf(f2,"sp %f \n",sp);
	int avgdelta=5; //the average duration of infection
	fprintf(f2,"avgdelta %d \n",avgdelta);
	float rho=0.25; //the average recovery time from the disease until susceptible
	fprintf(f2,"rho %f \n",rho);
	float omega=0.5; //the average waning period for vaccination
	fprintf(f2,"omega %f \n",omega);
	float eps=0.7;// the vaccine efficacy
	fprintf(f2,"eps %f \n",eps);
	float emin=0.65; // minimum perceived vaccine efficacy
	fprintf(f2,"emin %f \n",emin);
	float emax=0.90; // maximum perceived vaccine efficacy
	fprintf(f2,"emax %f \n",emax);
	float lambda=0.5; // the weight assignned to personal experiences 
	fprintf(f2,"lambda %f \n",lambda);
	float gamma=0.01; // probability of vaccine comoplications
	fprintf(f2,"gamma %f \n",gamma);
	float QALY=50000;
	fprintf(f2,"QALY %f \n",QALY);
	float maxcinf=0.15*0.023*QALY; // maximum cost of infection
	fprintf(f2,"maxcinf %f \n",maxcinf);
	float maxcvac=QALY*0.0023; //additional cost of vaccination
	fprintf(f2,"maxcvac %f \n",maxcvac);
	float mincvac=45; //minimum cost of vaccination
	fprintf(f2,"mincvac %f \n",mincvac);
	float m=0.3; //the memeory decay rate
	fprintf(f2,"m %f \n",m);
	float xi=15; //the memory decy factor when not vaccinate
	fprintf(f2,"m %f \n",xi);
	float sigma=0.5; //the imitation probability
	fprintf(f2,"sigma %f \n",sigma);
	float nu=0.5;// the preferecne of imitating contact
	fprintf(f2,"nu %f \n",nu);
	float alpha_m=0.00035; //mean probability of experiencing influenza like illness
	fprintf(f2,"alpha_s;	%.15f; \n", alpha_m);
	float alpha_s=0.000035; //stdev of probability of experiencing influenza like illness
	float alpha; // the probability of experiencing ili
	float beta=0.5; //probability individual thinks they experienced influenza when actually it was ili
	fprintf(f2,"beta;	%.15f; \n", beta);	
	float psi=0.7; // the probabiloity that symptoms appear for influenza (i.e. Ti decreases to 0 and don't vaccinate this year)
	fprintf(f2,"psi;	%.15f; \n", psi);	
	int Svac=0; //start of vaccination period
	fprintf(f2,"Svac %d \n",Svac);
	int Evac=121; //end of vaccination period
	fprintf(f2,"Evac %d \n",Evac);
	int dnVac=20; // the number of individuals to be contacted each day for random and nearest neighbout vaccination 
	fprintf(f2,"dnVac %d \n",dnVac);
	float VMI=0;// the value of the monetary incentives 
	if(inc==1)
	{
		VMI=20; 
	}
	else if(inc==2)
	{
		VMI=50;
		inc=1;
	}
	fprintf(f2,"VMI %f \n",VMI);
	float NSS; //number of super-spreaders
	float NSSC; // number of super-spreader contacts
	int SSNN; // number of super-spreaders not already contacted for vaccination as individuals nearest neighbour
	int SSRV; // number of super-spreaders selected as random
	int year; // year counter
	int day; //day counter
	int Sim; // simulation counter
	int ii; //counter
	int jj; //counter
	int kk; //counter	
	int zz; // counter
	float ppv=0; // probability of seeking vaccination this year (changes value in the code)
	float ftemp;// temporary float
	float ftemp3; //temporary float
	float ftemp2; // temporary float
	long double dtemp;// temporary double
	long double dtemp2; 
	long double dtemp3;
	int itemp; //temporary integer
	int itemp2; //temporary integer
	int itemp3; //temporary integer
	int itemp4; //temporar yinteger
	float PVPN; //Payoff to vaccinate minus the payoff not to vaccinate
	float SSV; // number of super-spreaders vaccinating each season
	float SSI; // number of super spreaders infected each Season
	float SSSV; // number of super-spreader contacts vaccinating each season
	float I; //total number of infected in the population after each season
	float V; //total number of vaccinating in the population after each season
	float NumI; // Number of infected in population
	float NumV; // Number of people that vaccinated in the population
	float LastI;// last years total infection in the population
	float k; //paramter used in creation of the contact network
	int temp;
	int countloop;
	int arac;	
	int pr[dnVac];
	int Ro95;
	float phi[4];
	float dj=0.5;// probability of recommened a friend for the page rank process
	fprintf(f2,"dj %f \n",dj);	
	Ro95=0;	
	
	ftemp2=(pow(q*avgdelta,Ro95)*exp(-q*avgdelta))/factorial((double)Ro95);
	while(ftemp2<0.95)
	{
		Ro95++;
		ftemp2=ftemp2+((pow(q*avgdelta,Ro95)*exp(-q*avgdelta))/factorial((double)Ro95));
	}
	for(Sim=0;Sim<NSims;Sim++)
	{
		//remove /* if would like to apply multivariate sensitivity analysis
		/*emin=lognorm(.35,sqrt(0.001),0.0001); 
		emax=1-lognorm(0.1,sqrt(0.003),0.0001);
		lambda=lognorm(0.5, sqrt(0.001),0.0001);
		sigma=lognorm(0.5, sqrt(0.001),0.0001);
		nu=lognorm(0.5, sqrt(0.001),0.0001);
		omega=lognorm(0.5, sqrt(0.001),0.0001);
		maxcinf=0.15*lognorm(0.023, 0.0005, 0.00001)*QALY;
		maxcvac=lognorm(0.002, 0.0001, 0.00001)*QALY;
		mincvac=45*100*lognorm(0.01, 0.001, 0.00001);
		eps=1-lognorm(.3,sqrt(0.001),0.0001);
		gamma=lognorm(.01,sqrt(0.000001),0.00001);
		m=lognorm(.25,sqrt(0.00095),0.0001);
		psi=1-lognorm(.3,sqrt(0.001),0.0001);*/
		k=39;
		for(ii=0;ii<N;ii++)
		{
			if(Sim>0)
			{
				free(ind[ii].contact);
			}
			ftemp=urand();
			itemp=0;
			ftemp2=(pow(k,itemp)*exp(-k))/factorial((double)itemp);
			while(ftemp2<ftemp)
			{
				itemp++;
				ftemp2=ftemp2+((pow(k,itemp)*exp(-k))/factorial((double)itemp));
			}
			ind[ii].Q=itemp;
			ind[ii].contact=(int*) malloc(itemp*sizeof(int));
			ind[ii].Qcn=0;			
		}	
		for (ii=0;ii<N;ii++) // build social network
		{
			countloop=0;
			while ((ind[ii].Qcn<ind[ii].Q) && (countloop <=50))
			{
				countloop=0;
				temp=rand()%(N-(ii))+ii;
				arac=0;	 // used to see if someone is already a contact
				// to see if someone is already a contact
				for (jj=0;jj<ind[ii].Qcn;jj++)
				{
					if (temp==ind[ii].contact[jj])
						arac=1;			// the person choosen is already a contact
				}
				// keeps choosing a new contact until it is not itself, already a contact and if the person choosen already has the max number of contacts
				while (((temp == ii) || (arac==1) || (ind[temp].Qcn==ind[temp].Q)) && (countloop<=50))
				{	
					countloop++;
					arac=0;
					temp=rand()%(N-ii)+ii;
					if (temp !=ii)
					{
						for (kk=0;kk<ind[ii].Qcn;kk++)
						{
							if (temp==ind[ii].contact[kk])
								arac=1;
						}
					}					
				}
				if (countloop<=50)
				{
					ind[ii].contact[ind[ii].Qcn]=temp;
					ind[temp].contact[ind[temp].Qcn]=ii;
					// reducing the number of contacts the two have
					ind[ii].Qcn++;
					ind[temp].Qcn++;
				}
			}
		}			
		// Connect remaing nodes that don't have required degree
		for(ii=0;ii<N;ii++) //build the social contact network
		{
			for(jj=0;jj<N;jj++)
			{
				if((ii !=jj) && (ind[ii].Q != ind[ii].Qcn) && (ind[jj].Q !=ind[jj].Qcn))
				{
					itemp=0;
					itemp2=ind[ii].Qcn;
					for(kk=0;kk<itemp2;kk++)
					{
						if(ind[ii].contact[kk]==jj)
						{
							itemp=1;
						}
					}
					if(itemp==0)
					{
						ind[ii].contact[itemp2]=jj;
						itemp2++;
						ind[ii].Qcn=itemp2;
						itemp2=ind[jj].Qcn;
						ind[jj].contact[itemp2]=ii;
						itemp2++;
						ind[jj].Qcn=itemp2;
					}
				}
			}
		}
		//Set intital conidtions for simulation
		LastI=0; // number of infections last season
		NSS=0; //set number of suoper-spreaders to zero
		NSSC=0; // set number of super-spreader contacts to zero
		NumI=0;	 //number of current infectious people
		//Double check the number of contacts for each individual in case we where not able to satisfy all the degree numbers
		for(ii=0;ii<N;ii++)
		{
			ind[ii].Q=ind[ii].Qcn;
			fprintf(f12,"%d \n",ind[ii].Q);
		}
		// Determine who are super-spreaders				
    	for (ii=0;ii<N;ii++)
		{
			//set the parameters for each individual
			ind[ii].peps=eps;
			if(hip==1)
			{
				itemp=0;
				while((itemp==0)||(itemp>=15))
				{
					itemp=0;
					ftemp=urand();
					ftemp2=(pow(avgdelta,itemp)*exp(-avgdelta))/factorial((double)itemp);
					while(ftemp2<ftemp)
					{
						itemp++;
						ftemp2=ftemp2+((pow(avgdelta,itemp)*exp(-avgdelta))/factorial((double)itemp));
					}		
				}
				ind[ii].delta=itemp;
			}
			else
			{
				ind[ii].delta=avgdelta;
			}
			if(htran==1)
			{
				ind[ii].avgp=lognorm(avgp,avgp*0.4,0.00001);
			}
			else
			{
				ind[ii].avgp=avgp;
			}
			ind[ii].sp=sp;		
			ind[ii].TI=365*100000;
			ind[ii].TC=365*100000;
			ind[ii].state=0;
			ind[ii].ctv=0;
			ind[ii].avgpeps=0;
			ind[ii].VacLast=0; // set last season vaccination to zero
			ind[ii].super=SS((double)Ro95, (double)ind[ii].Q,(double)ind[ii].avgp,ind[ii].delta); //input (Ro, # contacts, transmission for individual, delta period for individual) avg infection determine whether or not the individual is a super-spreader
			if(ind[ii].super==2)
			{
				fprintf(f16,"%d;%f;%d \n", ind[ii].Q, ind[ii].avgp, ind[ii].delta);
			}
		}
		for(ii=0;ii<N;ii++)
		{
			if(ind[ii].super>=2)
			{ 
				NSS++;
				for(jj=0;jj<ind[ii].Q;jj++) // Determine who are contacts of super-spreaders
				{
					if(ind[ind[ii].contact[jj]].super==0) //Change if there is no current realtion to a super-spreader
					{
						NSSC++;
						ind[ind[ii].contact[jj]].super=1; //now a contact of a super spreader
					}
					if(ind[ind[ii].contact[jj]].super==2) //Change if there is no current realtion to a super-spreader
					{
						NSSC++;
						ind[ind[ii].contact[jj]].super=3; //now a contact of a super spreader, as well as being a super-spreader
					}
				}
			}
		}
		fprintf(f11,"%f %f \n",NSS,NSSC);
		for(year=0;year<(Nyears+TransYear);year++)
		{
			NumV=0; //number of people that vaccinated (set to zero each season because vaccination starts over again each year
			SSV=0; // return vaccination of super-spreaders back to zero for this season
			SSI=0; //return number of super-spreaders infected this season back to zero
			SSSV=0; // return numer self vaccination supper spreaders to zero 
			V=0; // return number vaccinating this seaosn back to zero			
			I=0; // return number infected this seaosn back to zero	
			SSNN=0; // number of super-spreaders not alreadty contacted for vaccination as individuals nearest neighbour
			SSRV=0; //number of super-spreaders selected at random for vaccination that have already not vaccinated
			for(ii=0;ii<N;ii++)
			{
				if( (ind[ii].state==2)&& (urand()<rho)) //check to see if the individual loses natural immunity
				{
					ind[ii].state=0;	
				}
				else if((ind[ii].state==3)&&(urand()<omega)) //check to see if individual loses vaccine immunity
				{
					ind[ii].state=0;
				}
				if(ind[ii].VacLast==1)
				{
					ind[ii].peps=ind[ii].peps*exp(-m)+(1-exp(-m))*emax;
					ind[ii].pepstemp=ind[ii].peps;
				}
				if(ind[ii].VacLast==0)
				{
					ind[ii].peps=ind[ii].peps*exp(-m/xi)+(1-exp(-m/xi))*emax;
					ind[ii].pepstemp=ind[ii].peps;
				}
				ind[ii].TI++;
				ind[ii].TC++;
				ind[ii].TItemp=ind[ii].TI;
				ind[ii].TCtemp=ind[ii].TC;
				ind[ii].pepstemp=ind[ii].peps;
				ind[ii].VacLast=0;		
				ind[ii].MI=0;
				ind[ii].Vday=-1;//
				ind[ii].ctv=2; // set to two to determine that have not been choosen to vaccinate
			}	
			//imitation (learning process) and determine whether or not they vaccinate
			if(year>0)
			{
				phi[0]=lognorm(0.0593,sqrt(0.00015),0.0001);
				phi[1]=lognorm(0.1842,sqrt(0.00072),0.0001);
				phi[2]=lognorm(0.1012,sqrt(0.00038),0.0001);
				phi[3]=lognorm(0.0283,sqrt(0.00004),0.0001);
				for(ii=0;ii<N;ii++)
				{
					PVPN=ind[ii].peps*(lambda*maxcinf*exp(-m*ind[ii].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[ii].TC/365);
					ftemp2=PVPN;
					if(urand()<sigma) 
					{	
						// for imitation only of the random netwrok individual				
						itemp=rand()%N;
						while(itemp==ii)
						{
							itemp=rand()%N;
						}
						//randomly chosen contact imitation
						itemp2=ind[ii].contact[rand()%ind[ii].Q];
						ind[ii].TI=(1-sigma)*ind[ii].TItemp+sigma*((1-nu)*ind[itemp].TItemp+nu*ind[itemp2].TItemp);
						ind[ii].TC=(1-sigma)*ind[ii].TCtemp+sigma*((1-nu)*ind[itemp].TCtemp+nu*ind[itemp2].TCtemp);
						ind[ii].peps=(1-sigma)*ind[ii].pepstemp+sigma*((1-nu)*ind[itemp].pepstemp+nu*ind[itemp2].pepstemp);
						PVPN=ind[ii].peps*(lambda*maxcinf*exp(-m*ind[ii].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[ii].TC/365);						
					}	
					ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
					if(urand()<ppv)
					{
						ftemp=urand();
						if(ftemp<(phi[0]/(phi[0]+phi[1]+phi[2]+phi[3]))) //vaccination in Sept
						{
							ind[ii].ctv=0; // the individual will vaccinate this season
							ind[ii].Vday=rand()%30; //randomly chooses a day for the individual to vaccinate 
						}
						else if (ftemp<((phi[0]+phi[1])/(phi[0]+phi[1]+phi[2]+phi[3]))) //Vaccination in Oct
						{
							ind[ii].ctv=0; // the individual will vaccinate this season
							ind[ii].Vday=rand()%31+30; //randomly chooses a day for the individual to vaccinate 
						}
						else if (ftemp<((phi[0]+phi[1]+phi[2])/(phi[0]+phi[1]+phi[2]+phi[3]))) //Vaccination in Nov
						{
							ind[ii].ctv=0; // the individual will vaccinate this season
							ind[ii].Vday=rand()%30+61; //randomly chooses a day for the individual to vaccinate 
						}
						else //Vaccination in Dec
						{
							ind[ii].ctv=0; // the individual chooses to vaccinate this season
							ind[ii].Vday=rand()%31+91; //randomly chooses a day for the individual to vaccinate 
						}
						if(ind[ii].super>1) //number of super-spreaders willing to vaccinate
						{
							SSSV++;
						}
					}
				}		
				for (ii=0;ii<N;ii++) //update the individual's new information into the temproary terms
				{
					ind[ii].TItemp=ind[ii].TI;
					ind[ii].TCtemp=ind[ii].TC;
					ind[ii].pepstemp=ind[ii].peps;
				}				
			}
			countloop=0; // to count the weeks to inculate people with infection
			for(day=0;day<365;day++)
			{ 
				if((day-68==0)||(day-89==0)||(day-96==0)||(day-99==0)||(day-117==0)||(day-119==0)||(day-121==0))//randomly infecting people starting first part of november
				{
					itemp=rand()%N;
					while(ind[itemp].state!=0)
					{
						itemp=rand()%N;
					}
					if(urand()<psi)
					{
						ind[itemp].TI=-1; //set to -1 so that when it comes time to add TI=0								
						ind[itemp].ctv=3; // indicates that the individual has been infected this year and wont vaccinate this year									
						if(ind[itemp].VacLast==1)
						{
							ind[itemp].peps=ind[itemp].peps*(1-emin);
							ind[itemp].pepstemp=ind[itemp].peps;
							ind[itemp].VacLast=-1; // indicating that I vaccinated and it didn't work this year
						}
					}
					I++;
					if(ind[itemp].super>1)
					{
						SSI++;
					}
					ind[itemp].dc=1;
					ind[itemp].state=1;
					NumI++;					
					countloop++;
				}
				for(ii=0;ii<N;ii++) //set seasonal transmission for each individual and vaccination
				{
					
					ind[ii].TI++;
					ind[ii].TC++;
					ind[ii].p=ind[ii].avgp*(1+ind[ii].sp*sin(2*PI*(day-76)/365));
					ind[ii].TItemp=ind[ii].TI;
					ind[ii].TCtemp=ind[ii].TC;
					ind[ii].pepstemp=ind[ii].peps;
				}
				if((day>=Svac)&&(day<=Evac)&&(year>0))
				{
					if(strat==1) //random vaccination
					{
						itemp=0;
						while(itemp<dnVac)
						{
							itemp2=rand()%N;
							itemp++;
							if(ind[itemp2].ctv==2)
							{
								ind[itemp2].MI=ind[itemp2].MI+VMI;
								PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
								ftemp2=PVPN;
								if(urand()<sigma) 
								{	
									// for imitation only of the random netwrok individual				
									itemp3=rand()%N;
									while(itemp3==itemp2)
									{
										itemp3=rand()%N;
									}
									//randomly chosen contact imitation
									itemp4=ind[itemp2].contact[rand()%ind[itemp2].Q];
									ind[itemp2].TI=(1-sigma)*ind[itemp2].TItemp+sigma*((1-nu)*ind[itemp3].TItemp+nu*ind[itemp4].TItemp);
									ind[itemp2].TC=(1-sigma)*ind[itemp2].TCtemp+sigma*((1-nu)*ind[itemp3].TCtemp+nu*ind[itemp4].TCtemp);
									ind[itemp2].peps=(1-sigma)*ind[itemp2].pepstemp+sigma*((1-nu)*ind[itemp3].pepstemp+nu*ind[itemp4].pepstemp);
									PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
								}	
								PVPN=PVPN+ind[itemp2].MI;
								ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
								if(urand()<ppv)
								{
									ind[itemp2].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
									if(ind[itemp2].Vday>=0)
									{
										ind[itemp2].ctv=0; // the individual chooses to vaccinate this season
									}
									if(ind[itemp2].super>=2)
									{
										SSRV++;
									}
								}								
							}
							else if ((ind[itemp2].ctv==0)&&(inc==1))
							{
								ind[itemp2].MI=ind[itemp2].MI+VMI;
							}
						}
					}
					else if(strat==2) //nearest neighbour vaccination
					{
						itemp=0;
						while(itemp<dnVac)
						{
							itemp2=rand()%N;
							if(ind[itemp2].ctv==2)
							{
								PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
								ftemp2=PVPN;
								if(urand()<sigma) 
								{	
									// for imitation only of the random netwrok individual				
									itemp3=rand()%N;
									while(itemp3==itemp2)
									{
										itemp3=rand()%N;
									}
									//randomly chosen contact imitation
									itemp4=ind[itemp2].contact[rand()%ind[itemp2].Q];
									ind[itemp2].TI=(1-sigma)*ind[itemp2].TItemp+sigma*((1-nu)*ind[itemp3].TItemp+nu*ind[itemp4].TItemp);
									ind[itemp2].TC=(1-sigma)*ind[itemp2].TCtemp+sigma*((1-nu)*ind[itemp3].TCtemp+nu*ind[itemp4].TCtemp);
									ind[itemp2].peps=(1-sigma)*ind[itemp2].pepstemp+sigma*((1-nu)*ind[itemp3].pepstemp+nu*ind[itemp4].pepstemp);
									PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);									
								}	
								PVPN=PVPN+ind[itemp2].MI; //add on MI because they could have had incentives given previously to them in the season
								ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
								if(urand()<ppv)
								{
									ind[itemp2].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
									if(ind[itemp2].Vday>=0)
									{
										ind[itemp2].ctv=0; // the individual chooses to vaccinate this season
									}
									if(ind[itemp2].super>=2)
									{
										SSRV++;
									}
								}								
							}
							itemp=itemp+2;
							itemp3=rand()%ind[itemp2].Q; //contact of random individual	
							itemp3=ind[itemp2].contact[itemp3];
							if(ind[itemp3].ctv==2)
							{
								ind[itemp3].MI=ind[itemp3].MI+VMI;
								PVPN=ind[itemp3].peps*(lambda*maxcinf*exp(-m*ind[itemp3].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp3].TC/365);
								ftemp2=PVPN;
								if(urand()<sigma) 
								{	
									// for imitation only of the random netwrok individual				
									itemp2=rand()%N;
									while(itemp2==itemp3)
									{
										itemp2=rand()%N;
									}
									//randomly chosen contact imitation
									itemp4=ind[itemp3].contact[rand()%ind[itemp3].Q];
									ind[itemp3].TI=(1-sigma)*ind[itemp3].TItemp+sigma*((1-nu)*ind[itemp2].TItemp+nu*ind[itemp4].TItemp);
									ind[itemp3].TC=(1-sigma)*ind[itemp3].TCtemp+sigma*((1-nu)*ind[itemp2].TCtemp+nu*ind[itemp4].TCtemp);
									ind[itemp3].peps=(1-sigma)*ind[itemp3].pepstemp+sigma*((1-nu)*ind[itemp2].pepstemp+nu*ind[itemp4].pepstemp);
									PVPN=ind[itemp3].peps*(lambda*maxcinf*exp(-m*ind[itemp3].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp3].TC/365);
								}	
								PVPN=PVPN+ind[itemp3].MI;
								ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
								if(urand()<ppv)
								{
									ind[itemp3].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
									if(ind[itemp3].Vday>=0)
									{
										ind[itemp3].ctv=0; // the individual chooses to vaccinate this season
									}
									if(ind[itemp3].super>=2)
									{
										SSNN++;
									}
								}
							}
							else if ((ind[itemp3].ctv==0)&&(inc==1))
							{
								ind[itemp3].MI=ind[itemp3].MI+VMI;
							}
						}
					}
					else if(strat==3) // page rank vaccination
					{
						itemp=0;
						while(itemp<dnVac)
						{
							if(year==1)
							{
								itemp2=rand()%N;
								if(ind[itemp2].ctv==2)
								{
									ind[itemp2].MI=ind[itemp2].MI+VMI;
									PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
									ftemp2=PVPN;
									if(urand()<sigma) 
									{	
										// for imitation only of the random netwrok individual				
										itemp3=rand()%N;
										while(itemp3==itemp2)
										{
											itemp3=rand()%N;
										}
										//randomly chosen contact imitation
										itemp4=ind[itemp2].contact[rand()%ind[itemp2].Q];
										ind[itemp2].TI=(1-sigma)*ind[itemp2].TItemp+sigma*((1-nu)*ind[itemp3].TItemp+nu*ind[itemp4].TItemp);
										ind[itemp2].TC=(1-sigma)*ind[itemp2].TCtemp+sigma*((1-nu)*ind[itemp3].TCtemp+nu*ind[itemp4].TCtemp);
										ind[itemp2].peps=(1-sigma)*ind[itemp2].pepstemp+sigma*((1-nu)*ind[itemp3].pepstemp+nu*ind[itemp4].pepstemp);
										PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
									}	
									PVPN=PVPN+ind[itemp2].MI;
									ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
									if(urand()<ppv)
									{
										ind[itemp2].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
										if(ind[itemp2].Vday>=0)
										{
											ind[itemp2].ctv=0; // the individual chooses to vaccinate this season
										}
										if(ind[itemp2].super>=2)
										{
											SSRV++;
										}
									}								
								}
								else if ((ind[itemp2].ctv==0)&&(inc==1))
								{
									ind[itemp2].MI=ind[itemp2].MI+VMI;
								}
								pr[itemp]=itemp2;
							}
							else
							{
								if(urand()<dj)
								{
									itemp2=rand()%ind[pr[itemp]].Q;
									itemp2=ind[pr[itemp]].contact[itemp2];
									if(ind[itemp2].ctv==2)
									{
										ind[itemp2].MI=ind[itemp2].MI+VMI;
										PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
										ftemp2=PVPN;
										if(urand()<sigma) 
										{	
											// for imitation only of the random netwrok individual				
											itemp3=rand()%N;
											while(itemp3==itemp2)
											{
												itemp3=rand()%N;
											}
											//randomly chosen contact imitation
											itemp4=ind[itemp2].contact[rand()%ind[itemp2].Q];
											ind[itemp2].TI=(1-sigma)*ind[itemp2].TItemp+sigma*((1-nu)*ind[itemp3].TItemp+nu*ind[itemp4].TItemp);
											ind[itemp2].TC=(1-sigma)*ind[itemp2].TCtemp+sigma*((1-nu)*ind[itemp3].TCtemp+nu*ind[itemp4].TCtemp);
											ind[itemp2].peps=(1-sigma)*ind[itemp2].pepstemp+sigma*((1-nu)*ind[itemp3].pepstemp+nu*ind[itemp4].pepstemp);
											PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
										}	
										PVPN=PVPN+ind[itemp2].MI;
										ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
										if(urand()<ppv)
										{
											ind[itemp2].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
											if(ind[itemp2].Vday>=0)
											{
												ind[itemp2].ctv=0; // the individual chooses to vaccinate this season
											}
											if(ind[itemp2].super>=2)
											{
												SSNN++;
											}
										}								
									}
									else if ((ind[itemp2].ctv==0)&&(inc==1))
									{
										ind[itemp2].MI=ind[itemp2].MI+VMI;
									}									
									pr[itemp]=itemp2;
								}
								else 
								{
									itemp2=rand()%N;
									if(ind[itemp2].ctv==2)
									{
										ind[itemp2].MI=ind[itemp2].MI+VMI;
										PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
										ftemp2=PVPN;
										if(urand()<sigma) 
										{	
											// for imitation only of the random netwrok individual				
											itemp3=rand()%N;
											while(itemp3==itemp2)
											{
												itemp3=rand()%N;
											}
											//randomly chosen contact imitation
											itemp4=ind[itemp2].contact[rand()%ind[itemp2].Q];
											ind[itemp2].TI=(1-sigma)*ind[itemp2].TItemp+sigma*((1-nu)*ind[itemp3].TItemp+nu*ind[itemp4].TItemp);
											ind[itemp2].TC=(1-sigma)*ind[itemp2].TCtemp+sigma*((1-nu)*ind[itemp3].TCtemp+nu*ind[itemp4].TCtemp);
											ind[itemp2].peps=(1-sigma)*ind[itemp2].pepstemp+sigma*((1-nu)*ind[itemp3].pepstemp+nu*ind[itemp4].pepstemp);
											PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);											
										}	
										PVPN=PVPN+ind[itemp2].MI;
										ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
										if(urand()<ppv)
										{
											ind[itemp2].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
											if(ind[itemp2].Vday>=0)
											{
												ind[itemp2].ctv=0; // the individual chooses to vaccinate this season
											}
											if(ind[itemp2].super>=2)
											{
												SSRV++;
											}
										}								
									}
									else if ((ind[itemp2].ctv==0)&&(inc==1))
									{
										ind[itemp2].MI=ind[itemp2].MI+VMI;
									}
									pr[itemp]=itemp2;
								}
							}
							itemp++;							
						}
					}
					else if(strat==4) // nearest neighbour improved vaccination
					{
						itemp=0;
						while(itemp<dnVac)
						{
							itemp=itemp+2;
							itemp2=rand()%N;
							if(ind[itemp2].ctv==2)
							{
								ind[itemp2].MI=ind[itemp2].MI+VMI;
								PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
								ftemp2=PVPN;
								if(urand()<sigma) 
								{	
									// for imitation only of the random netwrok individual				
									itemp3=rand()%N;
									while(itemp3==itemp2)
									{
										itemp3=rand()%N;
									}
									//randomly chosen contact imitation
									itemp4=ind[itemp2].contact[rand()%ind[itemp2].Q];
									ind[itemp2].TI=(1-sigma)*ind[itemp2].TItemp+sigma*((1-nu)*ind[itemp3].TItemp+nu*ind[itemp4].TItemp);
									ind[itemp2].TC=(1-sigma)*ind[itemp2].TCtemp+sigma*((1-nu)*ind[itemp3].TCtemp+nu*ind[itemp4].TCtemp);
									ind[itemp2].peps=(1-sigma)*ind[itemp2].pepstemp+sigma*((1-nu)*ind[itemp3].pepstemp+nu*ind[itemp4].pepstemp);
									PVPN=ind[itemp2].peps*(lambda*maxcinf*exp(-m*ind[itemp2].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp2].TC/365);
								}	
								PVPN=PVPN+ind[itemp2].MI;
								ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
								if(urand()<ppv)
								{
									ind[itemp2].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
									if(ind[itemp2].Vday>=0)
									{
										ind[itemp2].ctv=0; // the individual chooses to vaccinate this season
									}
									if(ind[itemp2].super>=2)
									{
										SSRV++;
									}
								}
							}
							else if ((ind[itemp2].ctv==0)&&(inc==1))
							{
								ind[itemp2].MI=ind[itemp2].MI+VMI;
							}
							
							dtemp=0;
							for(kk=0;kk<ind[itemp2].Q;kk++)
							{
								dtemp=dtemp+ind[ind[itemp2].contact[kk]].Q;
							}
							dtemp2=urand();
							dtemp3=0;
							itemp3=-1;
							while(dtemp3<dtemp2)
							{
								itemp3++;
								dtemp3=dtemp3+(long double)ind[ind[itemp2].contact[itemp3]].Q/dtemp;
							}
							itemp3=ind[itemp2].contact[itemp3];
							if(ind[itemp3].ctv==2)
							{
								ind[itemp3].MI=ind[itemp3].MI+VMI;
								PVPN=ind[itemp3].peps*(lambda*maxcinf*exp(-m*ind[itemp3].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp3].TC/365);
								ftemp2=PVPN;
								if(urand()<sigma) 
								{	
									// for imitation only of the random netwrok individual				
									itemp2=rand()%N;
									while(itemp2==itemp3)
									{
										itemp2=rand()%N;
									}
									//randomly chosen contact imitation
									itemp4=ind[itemp3].contact[rand()%ind[itemp3].Q];
									ind[itemp3].TI=(1-sigma)*ind[itemp3].TItemp+sigma*((1-nu)*ind[itemp2].TItemp+nu*ind[itemp4].TItemp);
									ind[itemp3].TC=(1-sigma)*ind[itemp3].TCtemp+sigma*((1-nu)*ind[itemp2].TCtemp+nu*ind[itemp4].TCtemp);
									ind[itemp3].peps=(1-sigma)*ind[itemp3].pepstemp+sigma*((1-nu)*ind[itemp2].pepstemp+nu*ind[itemp4].pepstemp);
									PVPN=ind[itemp3].peps*(lambda*maxcinf*exp(-m*ind[itemp3].TI/365)+(1-lambda)*maxcinf*LastI)-mincvac-maxcvac*exp(-m*ind[itemp3].TC/365);
								}	
								PVPN=PVPN+ind[itemp3].MI;
								ppv=(atan(1360*PVPN/QALY)+(PI/2))/PI;
								if(urand()<ppv)
								{
									ind[itemp3].Vday=Vacday(day,phi[0],phi[1],phi[2],phi[3]);
									if(ind[itemp3].Vday>=0)
									{
										ind[itemp3].ctv=0; // the individual chooses to vaccinate this season
									}
									
									if(ind[itemp3].super>=2)
									{
										SSNN++;
									}
								}
							}
							else if ((ind[itemp3].ctv==0)&&(inc==1))
							{
								ind[itemp3].MI=ind[itemp3].MI+VMI;
							}
						}
					}					
					
					for(ii=0;ii<N;ii++) //vaccination occurring
					{
						if((ind[ii].ctv==0)&&(ind[ii].VacLast==0)) // if in vaccination period after first year and have been chosen to vaccinate but not vaccinated yet this year then can go through process
						{
							if(ind[ii].Vday==day) // see if the vaccinating or not
							{
								ind[ii].ctv=1;
								NumV++;
								V++;
								if(year>=TransYear)
								{
									fprintf(f10,"%d \n",ind[ii].Q); // record the degree of the vaccinated individual
									if(ind[ii].MI>0)
									{
										fprintf(f17,"%f; %d \n",ind[ii].MI,ind[ii].Q);
									}
								}
								ind[ii].MI=0;
								ind[ii].VacLast=1; //vaccinated last season
								if(ind[ii].super>1)
								{
									SSV++;
								}
								if((urand()<=eps)&&(ind[ii].state==0)) //see if vaccine was succesful
								{
									ind[ii].state=3; // indidivual suscesfully vaccinated
								}
								if(urand()<=gamma) //check to see if there is a complication from the vaccine
								{
									ind[ii].TC=-1; //individual suffered a complication, set to -1 so that when add next round TC=0
								}
							}							
						}
					}
				}				
				if(NumI>0)
				{
					for(ii=0;ii<N;ii++) // calculate probabilities of becoming infected
					{
						ind[ii].pn=1;
						for(jj=0;jj<ind[ii].Q;jj++)
						{
							itemp=ind[ii].contact[jj];
							if(ind[itemp].state==1)
							{
								ind[ii].pn=ind[ii].pn*(1-ind[itemp].p);
							}
						}
						ind[ii].pn=1-ind[ii].pn; //probability of indvidual becoming infected
					}
					for(ii=0;ii<N;ii++) //infect people
					{					
						if(ind[ii].state==0)
						{
							if(urand()<ind[ii].pn)
							{
								ind[ii].state=1; //changed to infected state since infected in previous loop
								if(urand()<psi)
								{
									ind[ii].TI=-1; //set to -1 so that when it comes time to add TI=0								
									ind[ii].ctv=3; // indicates that the individual has been infected this year and wont vaccinate this year									
									if(ind[ii].VacLast==1)
									{
										ind[ii].peps=ind[ii].peps*(1-emin);
										ind[ii].pepstemp=ind[ii].peps;
										ind[ii].VacLast=-1; // indicating that I vaccinated and it didn't work this year
									}
								}
								ind[ii].dc=0;
								if(year>=TransYear)
								{
									fprintf(f9,"%d \n",ind[ii].Q); //record the degree of the infected individual
								}
								NumI++;
								I++;
								if(ind[ii].super>1)
								{
									SSI++;
								}
							}
						}					
					}
					for(ii=0;ii<N;ii++) //remove individual to the recovered 
					{
						if((ind[ii].dc==ind[ii].delta)&&(ind[ii].state==1)) 
						{
							NumI--;
							ind[ii].state=2;
						}	
						ind[ii].dc++;
					}
				}
				if(year>=TransYear)
				{
					fprintf(f3I,"%f;",(float)(NumI/N));
					fprintf(f3V,"%f;",(float)(NumV/N));
				}
				alpha=lognorm(alpha_m,alpha_s,0.000001);
				for(ii=0;ii<N;ii++)
				{
					if(urand()<alpha)
					{
						if(urand()<beta) //determine whether or not they thought they had the flu
						{
							ind[ii].TI=-1;
							ind[ii].ctv=3; //thought it was flu so wont vaccinate this year
							if(ind[ii].VacLast==1)
							{
								ind[ii].peps=ind[ii].peps*(1-emin);
								ind[ii].pepstemp=ind[ii].peps;
								ind[ii].VacLast=-1; // indicating that I vaccinated and it didn't work this year
							}
						}
					}				
				}
			}
			LastI=exp(-m)*LastI+(float)I/N; // decaying the previous infections
			if(year>=TransYear)
			{
				for(ii=0;ii<N;ii++)
				{
					ind[ii].avgpeps=ind[ii].avgpeps+ind[ii].peps;
				}
				fprintf(f4I,"%f;",(float)(I/N));
				fprintf(f6I,"%f;",(float)(SSI/NSS));
				fprintf(f4V,"%f;",(float)(V/N));
				fprintf(f8V,"%f;",(float)(SSSV/NSS));
				fprintf(f6V,"%f;",(float)(SSV/NSS));
				fprintf(f13,"%f;",(float)(SSRV/NSS));
				fprintf(f14,"%f;",(float)(SSNN/NSS));
			}				
		} 
		
		for(ii=0;ii<N;ii++)
		{
			ind[ii].avgpeps=ind[ii].avgpeps/Nyears;
			fprintf(f15,"%f; %d \n",ind[ii].avgpeps,ind[ii].Q);
		}
		fprintf(f3I,"\n");
		fprintf(f4I,"\n");
		fprintf(f6I,"\n");
		fprintf(f3V,"\n");
		fprintf(f4V,"\n");
		fprintf(f6V,"\n");
		fprintf(f8V,"\n");
		fprintf(f13,"\n");
		fprintf(f14,"\n");
		
	}	
	fclose(f2);
	fclose(f3I);
	fclose(f3V);
	fclose(f4I);
	fclose(f4V);
	fclose(f6I);
	fclose(f6V);
	fclose(f8V);
	fclose(f9);
	fclose(f10);
	fclose(f11);
	fclose(f12);
	fclose(f13);
	fclose(f14);
	fclose(f15);
	fclose(f16);
	fclose(f17);
}
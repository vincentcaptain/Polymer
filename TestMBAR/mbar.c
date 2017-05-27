/*
 *
 *  author: David Limmer
 *  Written to preform MBAR anaylsis
 *  and compute psi and phi the large deviation functions
 *  Copyright 2017
 *
 */

#include<stdlib.h> 
#include<stdio.h> 
#include<math.h> 
#include<time.h> 
#define LINESIZE 500
#define SIDEMAX 100

/* Per Umbrella Variables */
double umb_beta[SIDEMAX];
double umb_spring[SIDEMAX];
double umb_loc[SIDEMAX];
double umb_num[SIDEMAX];
double f_old[SIDEMAX];
double f_new[SIDEMAX];
double exp_inner[SIDEMAX];
double exp_outer[SIDEMAX][200000];

/* Per Sample Variables */
double config_energy[SIDEMAX][200000];
double config_order[SIDEMAX][200000];
double config_order2[SIDEMAX][200000];
double config_denom[SIDEMAX][200000];
double config_weight[SIDEMAX][200000];

int hist[200];
double log_hist[200][200000];
double mino,maxo;
int bins,biastype;

/* State Point */
double temp_w,press,beta_w;
double tol=1e-7;

/* Functions */
int read_metafile(char *);
void read_data(char *,int ,int );
void printweights(char *,int);
void initial_f(int);
double compute_ui(int, int, int);
double mbar_sum(int);
double subtract_max_inner(int);
double subtract_max_outer(int);
double subtract_max_hist(int,int);
void printfi(int);
void printshifts(int);
void histogram(int,double,double,int);


/* MAIN PROGRAM */
int main(int argc, char** argv) {
	
	int num_windows;
	char metafile[LINESIZE];
	char weightsfile[LINESIZE];
	
    if(argc==7){
        sscanf(argv[1], "%s", &metafile);
        sscanf(argv[2], "%s", &weightsfile);
        sscanf(argv[3], "%lf", &mino);
        sscanf(argv[4], "%lf", &maxo);
        sscanf(argv[5], "%d", &bins);
        sscanf(argv[6], "%d", &biastype);
    
    
    
	num_windows = read_metafile(metafile);
	
	initial_f(num_windows);
	
	int iteration=0;
	double error=1;
	
	while(error>tol){
		
		iteration++;
		
		error=mbar_sum(num_windows);
		
		if(iteration%10==0) printf("Interation %d: Error %lf\n",iteration,error);
	}
	
	printweights(weightsfile,num_windows);
	if(biastype==1) printshifts(num_windows);
	printfi(num_windows);	
	
	histogram(num_windows,mino,maxo,bins);			
    
    }
    else{
        printf("Syntax: metafile.txt weightsfile.txt minQ maxQ bins bias_type (1=linear 2=harmonic) \n");
    }
    return 0;
}

/* Read In metadata file*/

int read_metafile(char *metafile){

	FILE * file;
	file=fopen(metafile,"r");
	
	char *line;
	char filename[LINESIZE];
	int num_windows,vals;
	double loc,spring,num_points,temp;
	
	num_windows = 0;
	line = malloc(sizeof(char)*LINESIZE);

	// make sure we're at the beginning of the file
	rewind(file);
	line = fgets(line,LINESIZE,file);

	while (line != NULL)
    {
        if(biastype==1){
            vals = sscanf(line, "%s %lf %lf", filename, &spring, &num_points);
            
            umb_beta[num_windows]=1.;
            umb_spring[num_windows]=spring;
            umb_num[num_windows]=num_points;
            
            read_data(filename,num_points,num_windows);
            num_windows++;
            line = fgets(line,LINESIZE,file);

        }
        if(biastype==2){
            vals = sscanf(line, "%s %lf %lf %lf", filename, &loc, &spring, &num_points);
		
            umb_beta[num_windows]=1.;
            umb_spring[num_windows]=spring;
            umb_loc[num_windows]=loc;
            umb_num[num_windows]=num_points;
            
            read_data(filename,num_points,num_windows);
            num_windows++;	    
            line = fgets(line,LINESIZE,file);
        }
    }
	free(line);
	
	return num_windows;
}

/* Read In sample file*/

void read_data(char *filename,int num_points,int num_windows){

	int i;
	char *line;
	double time,energy,volume,order,order2;
	double space,space1;
	
	FILE * file;
	file=fopen(filename,"r");
	line = malloc(sizeof(char)*LINESIZE);

	// make sure we're at the beginning of the file
	rewind(file);
	line = fgets(line,LINESIZE,file);

	for(i=0; i<num_points; i++){
	
		sscanf(line, "%lf %lf", &time, &order);
		
		config_energy[num_windows][i]=0.;
		config_order[num_windows][i]=order;
		config_order2[num_windows][i]=order;
		line = fgets(line,LINESIZE,file);
	}
		
	free(line);
}

/* Print weights file*/

void printshifts(int num_windows){

	int i;
	FILE * file;
	file=fopen("psi.txt","w");
	
	fprintf(file,"#BiasParameter Psi\n");
	
    double shift=0;
    for(i=0; i<num_windows; i++){
        if(umb_spring[i]==0.0) {
         shift=f_new[i];
            break;
        }
    }
	for(i=0; i<num_windows; i++){
		
		fprintf(file,"%lf %lf\n",umb_spring[i],-(f_new[i]-shift));
				
	}

}

/* Initial f_i's file*/

void initial_f(int num_windows){
	int i;
	for(i=0; i<num_windows; i++) f_old[i]=0.;
}

double compute_ui(int window_n, int window, int config){
	
	double bias,utot;
	
    if(biastype==1){
        bias=-umb_spring[window_n]*(config_order[window][config]);

    }
    if(biastype==2){
        bias=umb_spring[window_n]*(config_order[window][config]-umb_loc[window_n])*(config_order[window][config]-umb_loc[window_n]);
    }
    
	utot=umb_beta[window_n]*(bias);
	
	return utot;

}

/* Workhorse function */

double mbar_sum(int num_window){

	int i,j,k,n,np;
	double ui,uk;
	double logsum,denom;
	double error;
	
	for(j=0; j<num_window; j++){
		
		np=umb_num[j];
		for(n=0; n<np; n++){

			for(k=0; k<num_window; k++){
				uk=compute_ui(k,j,n);
				exp_inner[k]=log(umb_num[k])+f_old[k]-uk;
			}
			
			config_denom[j][n]=subtract_max_inner(num_window);
				
		}
	}
	
	for(i=0; i<num_window; i++){
	
		logsum=0;
		
		for(j=0; j<num_window; j++){
		
			np=umb_num[j];
			for(n=0; n<np; n++){
				
				ui=compute_ui(i,j,n);

				exp_outer[j][n]= -ui-config_denom[j][n];
				
			}
		}
		
		f_new[i]=-subtract_max_outer(num_window);
	
	}
	
	error=0;

	for(i=0; i<num_window; i++){
		
		f_new[i]=f_new[i]-f_new[0];
		
	}

	
	for(i=0; i<num_window; i++){
		
		error+=(f_new[i]-f_old[i])*(f_new[i]-f_old[i]);
		f_old[i]=f_new[i];
		
	}
	
	return error;

}

/* Print delta f's file*/

void printfi(int num_windows){

	int i;
	for(i=0; i<num_windows; i++) printf("%d %lf\n",i,f_old[i]);
	
}

/* Print weights file*/

void printweights(char *weightsfile,int num_window){

	int i,j,Np,k,n;
	double energy,order,weight,uk,order2;
	FILE * file;
	file=fopen(weightsfile,"w");
	
	for(i=0; i<num_window; i++){
		Np=umb_num[i];

		for(j=0; j<Np; j++){
			
			energy=config_energy[i][j];
			order=config_order[i][j];
			
			for(k=0; k<num_window; k++){
				uk=compute_ui(k,i,j);
				exp_inner[k]=log(umb_num[k])+f_new[k]-uk;
			}
			
			weight=-1.*subtract_max_inner(num_window);
			
			config_weight[i][j]=weight-umb_beta[i]*(config_energy[i][j]);
			
			fprintf(file,"%lf %lf %lf %lf\n",config_order[i][j],config_order2[i][j],umb_beta[i],config_weight[i][j]);
		}		
	}
		
}

/* Compute Projection of the free energy (bin weights) */

void histogram(int num_window, double mino, double maxo,int bins){
	
	int k,n,np,i,id,kk,ii,num_bin;
	double delo,q,bias,order;
    double logp[bins];
    
	FILE * file;
	file=fopen("phi.txt","w");
	
	fprintf(file,"# OrderP Phi\n");
	
	delo=(maxo-mino)/bins;

	for(i=0; i<bins; i++){ 
		hist[i]=0.;
		for(ii=0; ii<10000; ii++){
			log_hist[i][ii]=0;
		}			
	}	
	
	for(k=0; k<num_window; k++){	
		np=umb_num[k];
		
		for(n=0; n<np; n++){
	
			order=config_order2[k][n];
			if(order<maxo && order>mino){
				//id=round((order-mino-delo/2.)/delo);
				id=round((order-mino)/delo);
				
				log_hist[id][hist[id]]+=config_weight[k][n];
				hist[id]++;
					
			}
		}

	}
    
    double maxlp=-1000;
    
	for(i=0; i<bins; i++){
			
		//q=mino+i*delo+delo/2.;
		q=mino+i*delo;

		num_bin=hist[i]-1;
		logp[i]=subtract_max_hist(num_bin,i);
        if(logp[i]>maxlp) maxlp=logp[i];
		
	}
    
    for(i=0; i<bins; i++) {
      q=mino+i*delo;
      fprintf(file,"%lf %lf\n",q,logp[i]-maxlp);
    }
	
}

/* Functions for taming overflow */

double subtract_max_inner(int num_window){
	
	int k;
	double logsum=0.;
	double max_arg=0.;
	
	for(k=0; k<num_window; k++){
			if(exp_inner[k]>max_arg) max_arg=exp_inner[k];
	}
	
	for(k=0; k<num_window; k++){
			logsum+=exp(exp_inner[k]-max_arg);
	}

	
	logsum=log(logsum)+max_arg;
	return logsum;
}

double subtract_max_outer(int num_window){
	
	int k,n,np;
	double logsum=0.;
	double max_arg=0.;
	
	for(k=0; k<num_window; k++){
		np=umb_num[k];
		
		for(n=0; n<np; n++){	
			if(exp_outer[k][n]>max_arg) max_arg=exp_outer[k][n];
		}
		
	}
	
	for(k=0; k<num_window; k++){
		np=umb_num[k];
		
		for(n=0; n<np; n++){
			logsum+=exp(exp_outer[k][n]-max_arg);
		}
	}	
	
	logsum=log(logsum)+max_arg;
	
	return logsum;
}

double subtract_max_hist(int num_bin, int i){
	
	int k;
	double logsum=0.;
	double max_arg=0.;
	
	for(k=0; k<num_bin; k++){
			if(log_hist[i][k]>max_arg) max_arg=log_hist[i][k];
	}
	
	for(k=0; k<num_bin; k++){
			logsum+=exp(log_hist[i][k]-max_arg);
	}

	logsum=log(logsum)+max_arg;

	return logsum;
}

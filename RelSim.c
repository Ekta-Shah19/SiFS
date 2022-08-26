#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<sys/timeb.h>

struct feature{
			int rank;
			char gene[100];
			double relevance;
			double similarity;
			double objective;
			double *frequency;
			int *discrete;
};

double **realarray_allocation(int,int);
void make_normalize(double **,int,int);
void discretize_features(double *,int *,int );
double compute_jaccard_similarity(double **,int ,int ,int );
void calculate_marginal_frequency(int *,double *,int,int);
void calculate_joint_frequency(int *,int *,double **,int,int,int);
void getObjective(double **,int *,int ,int ,int ,int ,struct feature *);
double calculate_mutual_information(double *,double *,double **,int,int);
void realarray_deallocation(double **,int);

int main(int argc,char *argv[])
{
	int i,j,k,l;
	int maxindex;
	int counter,temp;
	int nnodes,nedges;
	int no_of_nodes,no_of_samples;
	int no_of_genes_selected;
	int max_gene_count,no_eqv_classes;
	int cluster_node_count;
	int *index,*classMat;
	int mergecount,classCount;
	int node1,node2,wt;
	double max,min;
	double sum;
	double alpha;
	double maximum;
	double similarity;
	double avg_rel;
	double **mat,**dataMat;
	char **gene_subset;
	char file[1000],gene[200];
	char comand[1000];
	time_t t;
    struct timeb ti,tf;
	struct feature *Feature;
	FILE *fp,*fp1;
	
	(void)ftime(&ti);
	fp=fopen(argv[1],"r");
	if(fp==NULL)
	{
		printf("Error opening file\n");
		exit(0);
	}
	else
	{
		fscanf(fp,"%d\t%d\t%d\n",&no_of_samples,&no_of_nodes,&classCount);
		
		dataMat=(double**)malloc(no_of_nodes*sizeof(double*));
		for(i=0;i<no_of_nodes;i++)
			dataMat[i]=(double *)calloc(no_of_samples,sizeof(double));
		classMat=(int*)malloc(no_of_samples*sizeof(int));
		
		no_eqv_classes=3;
		Feature=(struct feature *)malloc(sizeof(struct feature)*no_of_nodes);
		for(i=0;i<no_of_nodes;i++)
		{
			Feature[i].frequency=(double *)malloc(sizeof(double)*no_eqv_classes);
			Feature[i].discrete=(int *)malloc(sizeof(int)*no_of_samples);
			Feature[i].rank=-1;
		}
		
		for(i=0;i<no_of_nodes;i++)
		{
			fscanf(fp,"%s",Feature[i].gene);
			for(j=0;j<no_of_samples;j++)
				fscanf(fp,"%lf",&dataMat[i][j]);
		}
		for(i=0;i<no_of_samples;i++)
			fscanf(fp,"%d",&classMat[i]);
	}
	fclose(fp);
	
	make_normalize(dataMat,no_of_samples,no_of_nodes);
	getObjective(dataMat,classMat,no_of_samples,no_of_nodes,classCount,no_eqv_classes,Feature);
	printf("Relevance Values Computed\n");
	
	for(i=0;i<no_of_nodes;i++)
		free(dataMat[i]);
	free(dataMat);
	free(classMat);
	
	
	fp=fopen(argv[2],"r");
	fscanf(fp,"%d\t%d\n",&nnodes,&nedges);
	if(nnodes!=no_of_nodes)
	{
		fclose(fp);
		printf("Error Encountered\n");
		exit(0);
	}
	else
	{
		mat=(double **)malloc(no_of_nodes*sizeof(double *));
		for(i=0;i<no_of_nodes;i++)
			mat[i]=(double *)calloc(no_of_nodes,sizeof(double));
		
		for(i=0;i<nedges;i++)
		{
			fscanf(fp,"%d\t%d\t%d\n",&node1,&node2,&wt);
			mat[node1][node2]=wt*0.001;
			mat[node2][node1]=mat[node1][node2];
		}
		fclose(fp);
	}
	printf("Reading Similarity Matrix %s\n",argv[2]);
	
	
	sprintf(comand,"wc -l %s > temp.txt",argv[3]);
	system(comand);
	
	fp=fopen("temp.txt","r");
	fscanf(fp,"%d",&mergecount);
	fclose(fp);
	
	gene_subset=(char **)malloc(mergecount*sizeof(char *));
	for(i=0;i<mergecount;i++)
		gene_subset[i]=(char *)malloc(100*sizeof(char));
	
	fp=fopen(argv[3],"r");
	for(i=0;i<mergecount;i++)
		fscanf(fp,"%s\n",gene_subset[i]);
	fclose(fp);
	
	for(i=0;i<no_of_nodes;i++)
		for(j=0;j<mergecount;j++)
			if(strcmp(Feature[i].gene,gene_subset[j])==0)
			{
				Feature[i].rank=0;
				break;
			}
	
	avg_rel=0.0;
	for(i=0;i<no_of_nodes;i++)
		if(Feature[i].rank==0)
			avg_rel+=Feature[i].relevance;
	avg_rel=avg_rel/(double)mergecount;
	printf("Avg. Relevance: %lf\n",avg_rel);
	
	counter=0;
	for(i=0;i<no_of_nodes;i++)
	{
		if((Feature[i].rank==0)&&(Feature[i].relevance<avg_rel))
			Feature[i].rank=-1;
		if((Feature[i].rank==0)&&(Feature[i].relevance>=avg_rel))
			counter++;
	}
	printf("Reduced Set Size: %d\n",counter);
	
	max_gene_count=atoi(argv[4]);
	alpha=atof(argv[5]);
	
	printf("Maximum Number of Genes To Be Selected: %d\n",max_gene_count);
	printf("Alpha: %lf\n",alpha);
	
	index=(int *)malloc(max_gene_count*sizeof(int));
	
	i=0;
	while(Feature[i].rank!=0)
		i++;
	printf("First Feature With Rank 0: %d\n",i);
	
	no_of_genes_selected=0;
	maximum=Feature[i].relevance;
	maxindex=i;
	for(i=0;i<no_of_nodes;i++)
	{
		if((Feature[i].rank==0)&&(Feature[i].relevance>maximum))
		{
			maximum=Feature[i].relevance;
			maxindex=i;
		}
	}
	index[no_of_genes_selected]=maxindex;
	Feature[maxindex].rank=1;
	Feature[maxindex].objective=Feature[maxindex].relevance;
	no_of_genes_selected=1;
	printf("Gene: %s\tRelevance: %lf\n",Feature[maxindex].gene,maximum);
	
	while((no_of_genes_selected<max_gene_count)&&(counter>0))
	{
		for(j=0;j<no_of_nodes;j++)
		{
			if(Feature[j].rank==0)
			{
				for(i=0;i<no_of_nodes;i++)
				{
					if(Feature[i].rank==no_of_genes_selected)
					{
						similarity=compute_jaccard_similarity(mat,i,j,no_of_nodes);
						if(similarity>0.0)
						{
							Feature[j].similarity=(((no_of_genes_selected-1)*Feature[j].similarity)+similarity)/(double)no_of_genes_selected;
							Feature[j].objective=(alpha*Feature[j].relevance)+((1-alpha)*Feature[j].similarity);
							maximum=Feature[j].objective;
						}
						else
							Feature[j].rank=-1;
						break;
					}
				}
			}
		}
		
		counter=0;maximum=0.0;
		for(j=0;j<no_of_nodes;j++)
		{
			if(Feature[j].rank==0)
			{
				counter++;
				if(Feature[j].objective>maximum)
				{
					maximum=Feature[j].objective;
					maxindex=j;
				}
			}
		}
		
		printf("Iteration: %d\tCounter:%d\tGene:%s\tRelevance:%lf\tSimilarity:%lf\tObjective:%lf\n",no_of_genes_selected,counter,Feature[maxindex].gene,Feature[maxindex].relevance,Feature[maxindex].similarity,Feature[maxindex].objective);
		
		index[no_of_genes_selected]=maxindex;
		no_of_genes_selected++;
		Feature[maxindex].rank=no_of_genes_selected;
	}
	(void)ftime(&tf);
	printf("\nTOTAL TIME REQUIRED=%d millisec\n",(int)(1000.0*(tf.time-ti.time)+(tf.millitm-ti.millitm)));
	
	printf("Output File: %s\n",argv[6]);
	fp=fopen(argv[6],"w");
	for(i=0;i<no_of_genes_selected;i++)
	{
		j=index[i];
		fprintf(fp,"%s\n",Feature[j].gene);
	}
	fclose(fp);
	
	for(i=0;i<no_of_nodes;i++)
		free(mat[i]);
	free(mat);
	free(Feature);
	
	return 0;
}

double compute_jaccard_similarity(double **mat,int indx1,int indx2,int no_of_nodes)
{
	int i,j,k;
	double num,denom;
	double sim;
	
	num=0.0;denom=0.0;
	for(i=0;i<no_of_nodes;i++)
	{
		if(mat[indx1][i]<mat[indx2][i])
		{
			num+=mat[indx1][i];
			denom+=mat[indx2][i];
		}
		else
		{
			num+=mat[indx2][i];
			denom+=mat[indx1][i];
		}
	}
	
	if(num==0.0)
		sim=0.0;
	else
	{
		if(denom==0)
		{
			printf("Invalid Case for (%d,%d)\n",i,j);
			printf("Num: %lf\n",num);
			fflush(stdout);
			exit(0);
		}
		else
			sim=num/(double)denom;
	}
	
	return sim;
}

void realarray_deallocation(double **data,int row)
{
        int i;
        for(i=0;i<row;i++)
                free(data[i]);
        free(data);
}

double **realarray_allocation(int row,int column)
{
	int i;
	double **data;
	
	data=(double **)malloc(sizeof(double *)*row);
	for(i=0;i<row;i++)
		data[i]=(double *)malloc(sizeof(double)*column);
	return(data);
}


void getObjective(double **features_table,int *class_labels,int number_of_tuples,int number_of_features,int number_of_classes,int no_eqv_classes,struct feature *Feature)
{
	int i,j,k,l;
	double *object,*marginal_frequency;
	double **joint_frequency;
	double **sorted,temp,index;
	char fname[100];
	FILE *fp;

	marginal_frequency=(double *)malloc(sizeof(double)*number_of_classes);
	object=(double *)malloc(sizeof(double)*number_of_tuples);

	for(i=0;i<number_of_features;i++)
	{
		for(j=0;j<number_of_tuples;j++)
			object[j]=features_table[i][j];
		discretize_features(object,Feature[i].discrete,number_of_tuples);
		calculate_marginal_frequency(Feature[i].discrete,Feature[i].frequency,no_eqv_classes,number_of_tuples);
	}
	calculate_marginal_frequency(class_labels,marginal_frequency,number_of_classes,number_of_tuples);
	
	joint_frequency=realarray_allocation(no_eqv_classes,number_of_classes);
	for(i=0;i<number_of_features;i++)
	{
		calculate_joint_frequency(Feature[i].discrete,class_labels,joint_frequency,no_eqv_classes,number_of_classes,number_of_tuples);
		Feature[i].relevance=calculate_mutual_information(Feature[i].frequency, marginal_frequency, joint_frequency, no_eqv_classes, number_of_classes);
	}
	realarray_deallocation(joint_frequency,no_eqv_classes);
	
	free(object);
	free(marginal_frequency);
	
	return;
}

void discretize_features(double *object,int *discretized,int number_of_samples)
{
	int i;
	double mean,stddev;

	for(i=0;i<number_of_samples;i++)
	{
		if(object[i]>1)
			discretized[i]=0;
		else if(object[i]<-1)
			discretized[i]=1;
		else
			discretized[i]=2;
	}
}

void calculate_marginal_frequency(int *discrete,double *frequency,int number_of_classes,int number_of_tuples)
{
	int i;

	for(i=0;i<number_of_classes;i++)
		frequency[i]=0;
	for(i=0;i<number_of_tuples;i++)
		frequency[discrete[i]]++;
	for(i=0;i<number_of_classes;i++)
		frequency[i]/=number_of_tuples;
}

void calculate_joint_frequency(int *discrete,int *Discrete,double **joint_frequency,int no_eqv_classes,int number_of_classes,int number_of_tuples)
{
	int i,j;

	for(i=0;i<no_eqv_classes;i++)
		for(j=0;j<number_of_classes;j++)
			joint_frequency[i][j]=0;

	for(i=0;i<number_of_tuples;i++)
			joint_frequency[discrete[i]][Discrete[i]]++;
	for(i=0;i<no_eqv_classes;i++)
		for(j=0;j<number_of_classes;j++)
			joint_frequency[i][j]/=number_of_tuples;
}

double calculate_mutual_information(double *AFrequency,double *BFrequency,double **joint_frequency,int NoofAclass_of_samples,int NoofBclass_of_samples)
{
	int i,j,k;
	double AEntropy,BEntropy;
	double joint_entropy;
	double mutual_information;
	double normalized_mutual_information;
    double entropy_correlation_coefficient;

	AEntropy=0.0;
	for(i=0;i<NoofAclass_of_samples;i++)
		if(AFrequency[i]!=0.0)
			AEntropy=AEntropy-(AFrequency[i]*((double)log10((double)AFrequency[i])/(double)log10(2.0)));
	BEntropy=0.0;
	for(i=0;i<NoofBclass_of_samples;i++)
		if(BFrequency[i]!=0.0)
			BEntropy=BEntropy-(BFrequency[i]*((double)log10((double)BFrequency[i])/(double)log10(2.0)));
	
	joint_entropy=0.0;
	for(i=0;i<NoofAclass_of_samples;i++)
		for(j=0;j<NoofBclass_of_samples;j++)
			if(joint_frequency[i][j]!=0.0)
				joint_entropy=joint_entropy-(joint_frequency[i][j]*((double)log10((double)joint_frequency[i][j])/(double)log10(2.0)));
	mutual_information=AEntropy+BEntropy-joint_entropy;

	return(mutual_information);
}

void make_normalize(double **input_data,int number_of_tuples,int number_of_features)
{
   	int i,j;
	int choice;
   	double maximum,minimum;
	double mean,stddev;

	for(i=0;i<number_of_features;i++)
	{
		mean=0.0;
		for(j=0;j<number_of_tuples;j++)
			mean+=input_data[i][j];
		mean/=number_of_tuples;
		stddev=0.0;
		for(j=0;j<number_of_tuples;j++)
			stddev+=(input_data[i][j]-mean)*(input_data[i][j]-mean);
		stddev/=(number_of_tuples-1);
		stddev=(double)sqrt((double)stddev);
		for(j=0;j<number_of_tuples;j++)
			input_data[i][j]=(input_data[i][j]-mean)/stddev;
	}
	
	return;
}

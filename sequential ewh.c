#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define N 4     //No.of of hash tables
#define B 512    //bit size
#define NI 16    //(B/32) no.of integers
#define K 5     //no.of elements in key
#define M (int)pow(2,K)    
#define S (int)(1+K+(K*(K-1)/2))  //error vector hash index size

int ls_hash(int *n,int t,int val[N][K]);
int ham_dist(int *q,int *n);
void error_hash(int *q,int t,int *erhash,int val[N][K]);

int main()
{ 
 FILE *fp;
 int i,j,k,l,temp,temp2,n[NI],e=1,index,count,count2=0,co=0,s0=20,num=50,hamdist; 
 int q[NI]={1521,0,0,0,0,1213,1256,123,23,56,67,97,1223,256,4567,1};
 int **table,**ptr,**offset,*score,*data,*erhash,*alpha,*result,*res_dist;
 double start,end; 
 struct timeval time;
 int val[N][K]; 
 //printf("\nValue matrix\n");
  for(i=0;i<N;i++) {
    for(j=0;j<K;j++)
    {val[i][j]=rand()%B;
    //printf("%d ",val[i][j]); 
  }
  //printf("\n");
  }
 
 ptr=(int**)malloc(N*sizeof(int*));
 for(i=0;i<N;i++) 
   ptr[i]=(int*)malloc(M*sizeof(int));
 for(i=0;i<N;i++)
   for(j=0;j<M;j++)  ptr[i][j]=0;

 //fp=fopen("/home/mas/12/secbadar/project/text1.txt","r");
 //fscanf(fp,"%d",&count);
 //count=4194304;
 //count=1048576;
 count=67108864;
 //count=8192;
 data=(int*)malloc(count*sizeof(int));
 count=count/NI;//total 7 strings
 //while(fscanf(fp,"%d",&n)!=EOF)
 
 for(j=0;j<count;j++)
 {
  for(i=0;i<NI;i++) {
  //fscanf(fp,"%d",&n[i]);
  n[i]=rand()%10000;
  data[count2]=n[i]; count2++;
  }
  for(i=0;i<N;i++)
  {
   index=ls_hash(n,i,val);
   
  
   ptr[i][index]++;
  }
  
  //printf("\n");
 }//fclose(fp);
 

 
 for(i=0;i<N;i++) { 
  temp=0;
  for(j=0;j<M;j++) {
    temp2=ptr[i][j];
    ptr[i][j]=temp;
    temp+=temp2;
  }   
 }
 
 
 score=(int*)malloc(count*sizeof(int));    
 table=(int**)malloc(N*sizeof(int*));
 for(i=0;i<N;i++) 
   table[i]=(int*)malloc(count*sizeof(int));
 offset=(int**)malloc(N*sizeof(int*));
 for(i=0;i<N;i++) 
   offset[i]=(int*)malloc(M*sizeof(int)); 
 for(i=0;i<N;i++)
   for(j=0;j<M;j++)  offset[i][j]=0;   
 gettimeofday(&time,NULL);
 start=(time.tv_usec*1e-6)+time.tv_sec;    
 for(j=0;j<count;j++)
 {
  for(i=0;i<NI;i++)  n[i]=data[j*NI+i];
  for(i=0;i<N;i++)
  {
   index=ls_hash(n,i,val);
   table[i][ptr[i][index]+offset[i][index]]=j;
   offset[i][index]++;
  }
 }

 gettimeofday(&time,NULL);
 end=(time.tv_usec*1e-6)+time.tv_sec-start;
 printf("cpu time hashing taken is %lf\n",end);
 
 
 /*
 printf("\nData\n");
 for(i=0;i<count*NI;i++)
 printf("%d ",data[i]);
 */
 
 
 /*
 printf("\nTable\n");
 for(i=0;i<N;i++)
 {
 	for(j=0;j<count;j++)
 	{
 		printf("%d ",table[i][j]);	
 	}
 	printf("\n");
 
 }
 printf("\n");
 
 printf("\nOffset\n");
 for(i=0;i<N;i++)
 {
 	for(j=0;j<M;j++)
 	{
 		printf("%d ",offset[i][j]);
 	}
 	printf("\n");
 
 }
 
 printf("\n");
 printf("\nPtr\n");
 for(i=0;i<N;i++)
 {
 	for(j=0;j<M+1;j++)
 	{
 		printf("%d ",ptr[i][j]);
 	}
 	printf("\n");
 
 }
 */
 
 free(offset);
 
 /************Retrival***************/
 for(i=0;i<count;i++) score[i]=0;
 erhash=(int*)malloc(S*sizeof(int));
 alpha=(int*)malloc(S*sizeof(int));
 result=(int*)malloc(num*sizeof(int));
 res_dist=(int*)malloc(num*sizeof(int));
 alpha[0]=10;
 for(i=1;i<=K;i++) alpha[i]=5;
 for(i=K+1;i<=S;i++) alpha[i]=2;
 
 for(i=0;i<N;i++) 
 {
   error_hash(q,i,erhash,val);
    
    /*
    for(j=0;j<S;j++)
 	printf("%d ",erhash[j]);
 	printf("\n");
   */
   
   for(j=0;j<S;j++) {
     //index=ls_hash(q,i,val);
     for(k=ptr[i][erhash[j]];k<ptr[i][erhash[j]+1];k++)
       score[table[i][k]]+=alpha[j];
   }
 }
 
 /*
 for(i=0;i<count;i++) 
 printf("%d ",score[i]);

 printf("\n");
 */
 
 //for(i=0;i<count;i++) printf("%d\t",score[i]); printf("\n");
 gettimeofday(&time,NULL);
 start=(time.tv_usec*1e-6)+time.tv_sec;

 for(i=0;i<num;i++) res_dist[i]=B+2;
 for(i=0;i<count;i++) {
   if(score[i]>s0) { 
 
    hamdist=ham_dist(q,&data[i*NI]);
   // printf("%d ",data[i*NI]);
    
    k=0;
    while(hamdist>res_dist[k]) k++;
    for(j=num-2;j>=k;j--) { 
     result[j+1]=result[j]; 
     res_dist[j+1]=res_dist[j];
    }
    if(k<num) {
    result[k]=i;
    res_dist[k]=hamdist;
     }
  }
 }  
 
 
  //printf("\n\n");
 
 
 gettimeofday(&time,NULL);
 end=(time.tv_usec*1e-6)+time.tv_sec-start;
 printf("cpu time retrival is %lf\n",end);
 
 /*
 for (i=0;i<num;i++)
 printf("%d ",result[i]);
 
 printf("\n");
 */
 /*
 for (i=0;i<num;i++)
 printf("%d ",res_dist[i]);
 */
 

 //printf("Total count is %d small %d large %d\n",count,res_dist[0],res_dist[num-1]);
}

int ham_dist(int *q,int *n)
{
 int b[B],c[B],i,j,d=0,co=0;
 /*for(i=0;i<B;i++)
 {
  b[i]=q>>i & 1;
  c[i]=n>>i & 1;
 }*/
 for(i=0;i<NI;i++)
  for(j=0;j<32;j++) {
    b[co]=q[i]>>j & 1;
    c[co]=n[i]>>j & 1;
    co++;
  }
  /*
 for(i=0;i<B-1;i++)
 {
  b[i]=b[i]^b[i+1];
  c[i]=c[i]^c[i+1];
 }
 */
 for(i=0;i<B;i++)
 {
  if(b[i]!=c[i])
    d++;
 } 
  //printf("d is %d\n",d);
 return(d);
}

int ls_hash(int *n,int t,int val[N][K])
{ 
 int i,j,b[B],c=0;
   
 for(i=0;i<NI;i++)
  for(j=0;j<32;j++) {
    b[c]=n[i]>>j & 1;
    c++;
  }
  /*
 for(i=0;i<B-1;i++)
    b[i]=b[i]^b[i+1];
    */
    
 c=0;
 for(i=0;i<K;i++)
 { 
  if(b[val[t][i]]) 
    c=c+pow(2,i);
 }
 return(c);
}

void error_hash(int *q,int t,int *erhash,int val[N][K])
{
 int i,j,k,b[B],co=0,c[K],count=1+K;
 for(i=0;i<NI;i++)
  for(j=0;j<32;j++) {
    b[co]=q[i]>>j & 1;
    co++;
  }
  /*
 for(i=0;i<B-1;i++)
    b[i]=b[i]^b[i+1];
    */
 for(i=0;i<S;i++) erhash[i]=0;   
 for(i=0;i<K;i++)
 { 
  if(b[val[t][i]]) { c[i]=1; erhash[0]+=pow(2,i); }
  else   c[i]=0;
 } 
 for(i=0;i<K;i++) {
  for(j=0;j<K;j++)
  {
   if(i==j)  erhash[i+1]+=(!c[j])*pow(2,j);////Check this!!!!
   else erhash[i+1]+=c[j]*pow(2,j);
  } 
 }
 for(i=0;i<K;i++)
   for(j=0;j<i;j++) {
     for(k=0;k<K;k++)
     {
      if(k==i || k==j) erhash[count]+=(!c[k])*pow(2,k);
      else erhash[count]+=c[k]*pow(2,k);
     }
     count++;
   }  
 //for(i=0;i<S;i++) printf("%d\t",erhash[i]); printf("\n");
  
}

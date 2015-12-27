#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define N 4     //No.of of hash tables
#define B 512   //bit size
#define NI 16   //(B/32) no.of integers
#define K 5     //no.of elements in key
#define M 32    //(int)pow(2,K) No.of buckets in each hash table    
#define S 16    //(int)(1+K+(K*(K-1)/2))  //error vector hash index size



//#define SPT 64
 __constant__ int dB=B,dN=N,dNI=NI,dK=K,dM=M,dS=S;//,dSPT=SPT;

// Simple utility function to check for CUDA runtime errors

void checkCUDAError(const char* msg);

int ls_hash(int *n,int t,int val[N]);
void error_hash(int *q,int t,int *erhash,int *val,int *b);

__device__ int lsHash(int *n,int t,int *deVal);

__device__ volatile int sem = 0;

__device__ void acquire_semaphore(volatile int *lock){
  while (atomicCAS((int *)lock, 0, 1) != 0);
  }

__device__ void release_semaphore(volatile int *lock){
  *lock = 0;
  __threadfence();
  }

__global__ void Hash_table(int *devData,int *table,int *deVal,unsigned int *devPtr,unsigned int *devOff,int count,int dSPT)
{
	
  int i,j,n[16],index,gid=blockIdx.x*blockDim.x+threadIdx.x;
  
  unsigned int dummy;
  //extern  __shared__  int temp[];
  //int ptr[16],off[16];
  //for(i=0;i<dM;i++)  { ptr[i]=0; off[i]=0;} 
  dSPT=1;
  
  for(i=0;i<dSPT;i++)
  {
   for(j=0;j<dNI;j++)
   	//if(gid*dNI+i+j<2048)
      n[j]=devData[gid*dNI+i+j];
      
   for(j=0;j<dN;j++)
   {
	   index=lsHash(n,j,deVal);
	   dummy=atomicInc(&devPtr[index+j*(dM+1)],100000);  //devPtr[index]++;
  }
 }
  
}

__global__ void Hash_table2(int *devData,int *table,int *deVal,unsigned int *devPtr,unsigned int *devOff,int count,int dSPT)
{
	
int i,j,n[16],index,gid=blockIdx.x*blockDim.x+threadIdx.x;
  
  unsigned int dummy;
   
  dSPT=1;
  

  for(i=0;i<dSPT;i++)
  {
   for(j=0;j<dNI;j++)
   		//if(gid*dNI+i+j<4096)
      n[j]=devData[gid*dNI+i+j];
  
   
   
   for(j=0;j<dN;j++)
   {
	   index=lsHash(n,j,deVal);
	   
  
   //if((index+j*dM)<(dN*dM))
   dummy=atomicInc(&devOff[index+j*dM],100000);  //off[index]++;
   
   //if(((j*(dM+1)+index)<dN*(dM+1)) && (j*count+devPtr[j*(dM+1)+index]+dummy)< (dN*count))
   //table[t*count+devPtr[t*(dM+1)+index]+dummy]=gid*dSPT+i;    //table[i][ptr[i][index]+offset[i][index]]
   table[j*count+devPtr[j*(dM+1)+index]+dummy]=threadIdx.x+ blockIdx.x * blockDim.x;
  }
 }
  
}

__global__ void Score_update(int *table,unsigned int *devPtr,int *devScore,int *devErhash,int *devAlpha,int count)
{
 int i=0,b=(blockIdx.x)%dS,dummy, t=blockIdx.x/dS, c=blockDim.x;
 //t from 0 to N
 //b from 0 to s
 
 //if(t*dS+b <dS*dN)
 
 for(i=devPtr[(t*(dM+1)+devErhash[t*dS+b])%(N*(M+1))];i+threadIdx.x<devPtr[(t*(dM+1)+devErhash[t*dS+b]+1)%(N*(M+1))];i=i+c)          //ptr[i][erhash[j]]
 {
 //if(t*dS+b < dS)
 //if((t*count+i<dN*count))//&&(t*dS+b < dS))// && (table[t*count+i]<count))
 
  dummy=atomicAdd(&devScore[table[t*count+i+threadIdx.x]],devAlpha[t*dS+b]);    //score[table[i][k]]+=alpha[j];
  //dummy=atomicAdd(&devScore[table[t*count+i]],devAlpha[t*dS+b]);    //score[table[i][k]]+=alpha[j];
 }//*/
}




__global__ void Ham_dist(int *devData,int *devScore,int count,int s0,int *qHam,int *devRind,int *devHdist,int kn)
{
 int wid,gid=threadIdx.x+blockIdx.x*blockDim.x,off=gridDim.x*blockDim.x,tps=16,sn,str,b[32],sum=0,temp3,newvar;
 unsigned int i,j,dumy,last,counter=0;
 extern  __shared__  unsigned int res[];
  __shared__  unsigned int temp2[1500];
 
 
 __shared__ unsigned int s3,s[256];//check this size 
 if(threadIdx.x==0) {s3=0;s[0]=0;}
 __syncthreads();
 for(i=gid;i<count;i+=off) {                               
  if(devScore[i]>s0)  {                                    //If score is above s0, store the index of the string in shared memory
	 
   dumy=atomicInc(&s[0],count);
   
   temp2[dumy]=i;  
   
  }
 }
 
 __syncthreads();
 
 last=s[0];
 
 
 __syncthreads();
 off=blockDim.x/tps; wid=threadIdx.x%tps, newvar=threadIdx.x/tps;
 for(i=newvar;i<last;i+=off) {	
  
  sum=0;
  if(wid<dNI) 	
  {
   str=devData[temp2[i]*dNI+wid];	
   for(j=0;j<32;j++) b[j]=str>>j & 1;
 //  for(j=0;j<32-1;j++) b[j]=b[j]^b[j+1];
   for(j=0;j<32;j++) if(b[j]!=qHam[wid*32+j]) sum++;
  }
  

	
   s[threadIdx.x]=sum;                           //Find sum
   
   
  for (j=1;j<tps;j<<=1) {
     __syncthreads();
     temp3=(wid>=j) ? s[threadIdx.x-j]:0;
     __syncthreads();
     s[threadIdx.x]+=temp3;
   }
  
  
  
  __syncthreads();
  
  
   //sn=(threadIdx.x/tps)+(i*off);
 
 
  if(wid==tps-1)	///sum in end and to res[sn]
  {
 	sn=i; 
  res[sn]=s[threadIdx.x];
 // printf("%d ",res[sn]);
 }
 
 
 }
 
 
if (res[threadIdx.x]==0)
res[threadIdx.x]=dB+2;
 
 
 __syncthreads();
 for(i=0;i<last;i+=blockDim.x)
 {
  off=threadIdx.x+(i*blockDim.x);sum=0;	
  //off=threadIdx.x;sum=0;
  for(j=0;j<last;j++)
   if(res[off]>res[j]) sum++;
  if(sum<kn) 
  {
  devHdist[sum+kn*blockIdx.x+counter]=res[off];
  devRind[sum+kn*blockIdx.x+counter]=temp2[off];
  counter++;
 // printf("%d ",devHdist[sum+kn*blockIdx.x]); 
  }
 }
 

}


__global__ void Final_sort(int *devResult,int *devHdist,int *devRind,int kn)
{
 int i=0,j,sum,temp,gid=threadIdx.x+blockDim.x*blockIdx.x,strtemp;
 unsigned int counter=0;
 
 for(i=gid;i<kn*gridDim.x;i+=(blockDim.x*gridDim.x))
 {
  temp=devHdist[i];  sum=0;
  strtemp=devRind[i];
  for(j=i;j<kn*gridDim.x+i;j++)
     if(temp>devHdist[j%(kn*gridDim.x)]) 
     {	   
     sum++;
     }
     
  if(sum<kn) 
  {
  
  /*__syncthreads();
		if (threadIdx.x == 0)
		  acquire_semaphore(&sem);
		__syncthreads();
		  //begin critical section
		  // ... your critical section code goes here
  */while(devResult[sum]!=0) sum++; //check this line??
 devResult[sum+counter]=strtemp;
 counter++; 
  //printf("%d ",devResult[sum]);
  
  //end critical section
	/*	__syncthreads();
		if (threadIdx.x == 0)
		  release_semaphore(&sem);
		__syncthreads();
  */
  }
  
 }
}

int main( int argc, char** argv) 
{
cudaSetDevice(1);
 FILE *fp;
 int i,j,qHam[B],n[NI],count,count2=0,s0=20,num=50;   //index,temp,temp2,co=0,hamdist,k,l,e=1; 
 int q[NI]={1521,0,0,0,0,1213,1256,123,23,56,67,97,1223,256,4567,1};
 int *data,*erhash,*alpha,*result,*val,*score;//*res_dist;
 int *devData,*devScore,*table,*deVal,*devErhash,*devAlpha,*devqHam,*devHdist,*devRind,*devResult,*hosttable,*hostHdist;
 unsigned int *ptr,*devPtr,*off,*devOff,*hostoff,*hostptr;
 int n_threads,n_blocks,shared_mem_size,SPT,indexsid;
 double start;
 struct timeval time;
 val=(int*)malloc(sizeof(int)*N*K); 
 
 //#pragma omp parallel for default(none) private(j) shared(val) //schedule(auto) cannot have hre because of rand function!!
  for(i=0;i<N;i++) {
    for(j=0;j<K;j++)
    val[i*K+j]=rand()%B; //check openmp
  }
  
  
  
 ptr=(unsigned int*)malloc(N*(M+1)*sizeof(unsigned int*));
 off=(unsigned int*)malloc(N*(M+1)*sizeof(unsigned int*));
 
 hostoff=(unsigned int*)malloc(N*M*sizeof(unsigned int*));
 hostptr=(unsigned int*)malloc(N*(M+1)*sizeof(unsigned int*));
 hosttable=(int*)malloc(N*count*sizeof(int*));
 
 //fp=fopen("/home/mas/12/secbadar/project/text1.txt","r"); 
 //fscanf(fp,"%d",&count); 
 count=67108864;
 //count=1024;
 //count=8192;
 data=(int*)malloc(count*sizeof(int));
 count=count/NI;
 //while(fscanf(fp,"%d",&n)!=EOF)
 
 //check openmp here!
 for(j=0;j<count;j++)
 {
  for(i=0;i<NI;i++) {
    n[i]=rand()%10000;
  //fscanf(fp,"%d",&n[i]);
  data[count2]=n[i]; count2++;
  }
 }//fclose(fp);
   
 //n_threads=256;  n_blocks=N*(count/16);
 n_threads=1024;  n_blocks=4096;
 //n_threads=16;  n_blocks=4; //for 1024
 //n_threads=32;  n_blocks=16; //for 8192
 //For count=4194304 it is 256,1024  numthreads*numblocks=count/NI (new count)
 //For count=4096 it is 32,8
 //For count=1048576 it is 256,256
 SPT=(count*N)/(n_blocks*n_threads);
 dim3 grid_dim(n_blocks);
 dim3 block_dim(n_threads);

 cudaMalloc((void**)&devData,sizeof(int)*count*NI);        //allocating memory in device
 cudaMalloc((void**)&table,sizeof(int)*N*count);
 cudaMalloc((void**)&devPtr,sizeof(unsigned int)*N*(M+1));
 cudaMalloc((void**)&devOff,sizeof(unsigned int)*N*M);
 cudaMalloc((void**)&deVal,sizeof(int)*N*K);
 
// #pragma omp parallel for default(shared) schedule(auto)
 //for(i=0;i<N*(M+1);i++) ptr[i]=0; 
 memset(ptr,0,sizeof(ptr));
 
// #pragma omp parallel for default(shared) schedule(auto)
// for(i=0;i<N*M;i++) off[i]=0;   
 memset(off,0,sizeof(off));
 
 gettimeofday(&time,NULL);
 start=(time.tv_usec*1e-6)+time.tv_sec;

 cudaMemcpy(devPtr,ptr,sizeof(unsigned int)*N*(M+1),cudaMemcpyHostToDevice); 
 cudaMemcpy(devOff,off,sizeof(unsigned int)*N*M,cudaMemcpyHostToDevice);       
 cudaMemcpy(devData,data,sizeof(int)*count*NI,cudaMemcpyHostToDevice);   // copy host memory to device input array
 cudaMemcpy(deVal,val,sizeof(int)*N*K,cudaMemcpyHostToDevice);
 checkCUDAError("memcpy 1"); 
 //cudaThreadSynchronize();
 //shared_mem_size=sizeof(int)*(n_threads+1);
 //printf("\nGrid Dim or no. of blocks: %d\nBlock Dim of no. of threads: %d\nShared Memory Size: %d\n SPT: %d : %d : %d\n", n_blocks,n_threads,shared_mem_size,SPT, grid_dim.x, block_dim.x );
 
 
 
 int temp,temp2;
 /*
 for(j=0;j<count;j++)
 {
 for(i=0;i<NI;i++) {
  //fscanf(fp,"%d",&n[i]);
  n[i]= data[j*NI+i];
  }
 // #pragma omp parallel for default(shared) private(i,index,n) schedule(auto) 
  for(i=0;i<N;i++)
  {
   indexsid=ls_hash(n,i,val);
  // printf("%d ",indexsid);
   ptr[(i*(M+1))+ indexsid]++;
  }
 }
 
 for(i=0;i<N;i++) { 
  temp=0;
  for(j=0;j<M;j++) {
    temp2=ptr[i*(M+1)+j];
    ptr[i*(M+1)+j]=temp;
    temp+=temp2;
  }   
 }
 
 for(i=0;i<N;i++)
  for(j=0;j<M+1;j++)
  	printf("%d ",ptr[i*(M+1)+j]);
 */
 
 Hash_table<<<grid_dim,block_dim>>>(devData,table,deVal,devPtr,devOff,count,SPT);  //invoke gpu kernel
 
 checkCUDAError("kernel1-a invocation");
 
 /*
 #pragma omp parallel for default(shared) private(i,j,indexsid,n,temp,temp2) schedule(static) 


 for(j=448;j<count;j++)
 {
 for(i=0;i<NI;i++) {
  //fscanf(fp,"%d",&n[i]);
  n[i]= data[j*NI+i];
  }
 // #pragma omp parallel for default(shared) private(i,index,n) schedule(auto) 
  for(i=0;i<N;i++)
  {
   indexsid=ls_hash(n,i,val);
  // printf("%d ",indexsid);
  #pragma omp atomic	//awesome!!
   ptr[(i*(M+1))+ indexsid]++;
  }
 }
 
 */
 
 cudaThreadSynchronize();
 
 
 cudaMemcpy(hostptr,devPtr,sizeof(unsigned int)*N*(M+1),cudaMemcpyDeviceToHost);              
 cudaThreadSynchronize();
 
 
 for(i=0;i<N;i++)
  {
  	for(j=0;j<M;j++)
  	{
  		hostptr[(i*(M+1))+ j]+=ptr[(i*(M+1))+ j];
  		
  	}
 }
 
 
 
 


//unc#pragma omp parallel for default(shared) private(temp,temp2,j) schedule(auto)
 for(i=0;i<N;i++) { 
  temp=0;
  for(j=0;j<M;j++) {
    temp2=hostptr[i*(M+1)+j];
    hostptr[i*(M+1)+j]=temp;
    temp+=temp2;
  }   
 }
 
 cudaThreadSynchronize();
 cudaMemcpy(devPtr,hostptr,sizeof(unsigned int)*N*(M+1),cudaMemcpyHostToDevice); 
 //cudaThreadSynchronize();
 






 Hash_table2<<<grid_dim,block_dim>>>(devData,table,deVal,devPtr,devOff,count,SPT);  //invoke gpu kernel
 
 checkCUDAError("kernel1-b invocation");
 
 
 /*
 #pragma omp parallel for default(shared) private(i,j,indexsid,n) schedule(static) 
 for(j=448;j<count;j++)
 {
  for(i=0;i<NI;i++)  n[i]=data[j*NI+i];
  for(i=0;i<N;i++)
  {
   indexsid=ls_hash(n,i,val);
   #pragma omp barrier
   hosttableomp[(i*count)+hostptr[(i*(M+1))+indexsid]+off[(i*M) + indexsid]]=j;	///this line cannot be parallelised in openmp
   #pragma omp atomic
   off[i*M + indexsid]++;
  }
 }
*/
 
 
 
 
 
 cudaThreadSynchronize();
 
	
 gettimeofday(&time,NULL);
 start=(time.tv_usec*1e-6)+time.tv_sec-start;
 printf("gpu hashing time taken is %lf\n",start);
 
 //cudaMemcpy(hosttable,table,sizeof(int)*N*count,cudaMemcpyDeviceToHost);
 cudaMemcpy(hosttable,table,sizeof(int)*N*count,cudaMemcpyDeviceToHost);
 
 checkCUDAError("memcpy 2"); 
 cudaMemcpy(hostoff,devOff,sizeof(unsigned int)*N*M,cudaMemcpyDeviceToHost);
 checkCUDAError("memcpy 3"); 
 cudaMemcpy(hostptr,devPtr,sizeof(unsigned int)*N*(M+1),cudaMemcpyDeviceToHost);              
 checkCUDAError("memcpy 4"); 
 //cudaThreadSynchronize();
 
 
 /*
 for(i=0;i<N;i++)
  {
  	for(j=0;j<M;j++)
  	{
  		hostoff[(i*(M))+ j]+=off[(i*(M))+ j];
  		
  	}
 }
 
 
 for(i=0;i<N;i++)
  {
  	for(j=0;j<count;j++)
  	{
  		if(!hosttableomp[i*count+j])
  		hosttable[(i*(count))+ j]=hosttableomp[(i*(count))+ j];
  		
  	}
 }
 */
 
 
 
 /* Checking correctness of hashing
 */
 /*printf("\nTable\n");
 for(i=0;i<N;i++)
 {
 	for(j=0;j<count;j++)
 	{
 		printf("%d ",hosttable[i*count+j]);
 	}
 	printf("\n");
 
 }
 
 printf("\nOffset\n");
 for(i=0;i<N;i++)
 {
 	for(j=0;j<M;j++)
 	{
 		printf("%d ",hostoff[i*M+j]);
 	}
 	printf("\n");
 
 }
 
 
 printf("\nPtr\n");
 for(i=0;i<N;i++)
 {
 	for(j=0;j<M+1;j++)
 	{
 		printf("%d ",hostptr[i*(M+1)+j]);
 	}
 	printf("\n");
 
 }
  
  
  */
 

//*******************************************************************************************
/***************Retrival*******************/
 score=(int*)malloc(count*sizeof(int));
 
 //#pragma omp parallel for default(shared) schedule(auto)
 //for(i=0;i<count;i++)  score[i]=0;
 
 memset(score,0,sizeof(score));
 
 alpha=(int*)malloc(S*N*sizeof(int));
 erhash=(int*)malloc(S*N*sizeof(int));
 
 int m;
 
//unc #pragma omp parallel for default(none) private(i,m) shared(alpha) 
 for(m=0;m<N;m++)
 {
 
 alpha[0+m*S]=10;
//#pragma omp parallel for default(none) schedule(auto) shared(alpha,m) private(i)
 for(i=1;i<=K;i++) alpha[i+m*S]=5;
// #pragma omp parallel for default(none) schedule(auto) shared(alpha,m) private(i)
 for(i=K+1;i<=S;i++) alpha[i+m*S]=2;
 }
 
//unc #pragma omp parallel for default(none) shared(q,erhash, qHam, val) schedule(auto) private(i)
 for(i=0;i<N;i++)  error_hash(q,i,&erhash[i*S],val,qHam);
 
 /*To print error hash
 for(i=0;i<N;i++)
 {for(j=0;j<S;j++)
 {printf("%d ",erhash[i*S+j]);}
 printf("\n");
 }
 */
 
 
 

 cudaMalloc((void**)&devAlpha,sizeof(int)*S*N);
 cudaMalloc((void**)&devErhash,sizeof(int)*S*N);
 cudaMalloc((void**)&devScore,sizeof(int)*count);

  gettimeofday(&time,NULL);
 start=(time.tv_usec*1e-6)+time.tv_sec;
 cudaMemcpy(devAlpha,alpha,sizeof(int)*S*N,cudaMemcpyHostToDevice);
 cudaMemcpy(devErhash,erhash,sizeof(int)*S*N,cudaMemcpyHostToDevice);
 cudaMemcpy(devScore,score,sizeof(int)*count,cudaMemcpyHostToDevice); 
 //cudaThreadSynchronize();
 
 //n_threads=256;  n_blocks=N*S*(count/8);
 //n_threads=32;  n_blocks=64;//blocks N*S
 n_threads=1024;  n_blocks=64;//blocks N*S
 dim3 grid_dim2(n_blocks);
 dim3 block_dim2(n_threads);

 
 
Score_update<<<grid_dim2,block_dim2>>>(table,devPtr,devScore,devErhash,devAlpha,count);
 checkCUDAError("kernel2 invocation");
 
 
 
 ///openmp not required here pretty less work been done!!
 /*
 for(i=0;i<N;i++) 
 {
   error_hash(q,i,erhash,val);
    
    /*
    for(j=0;j<S;j++)
 	printf("%d ",erhash[j]);
 	printf("\n");
     
   for(j=0;j<S;j++) {
     //index=ls_hash(q,i,val);
     for(k=ptr[i][erhash[j]];k<ptr[i][erhash[j]+1];k++)
       score[table[i][k]]+=alpha[j];
   }
 }
 */
 
 
 
 cudaThreadSynchronize();
 
 
 /* To print score
 cudaMemcpy(score,devScore,sizeof(int)*count,cudaMemcpyDeviceToHost); 
 for(i=0;i<count;i++) 
 	printf("%d ",score[i]);
printf("\n");
 */
 
 
  
//*******************************************************************************************
 
 
 
 //n_blocks=(count/8);
 n_threads=1024;  n_blocks=4096;
 
 dim3 grid_dim3(n_blocks);
 dim3 block_dim3(n_threads);
 
 hostHdist=(int*)malloc(n_blocks*num*sizeof(int));
 
 cudaMalloc((void**)&devqHam,sizeof(int)*B);
 cudaMalloc((void**)&devRind,sizeof(int)*n_blocks*num);
 cudaMalloc((void**)&devHdist,sizeof(int)*n_blocks*num);
 cudaMemcpy(devqHam,qHam,sizeof(int)*B,cudaMemcpyHostToDevice);
 cudaMemcpy(devHdist,hostHdist,sizeof(int)*n_blocks*num,cudaMemcpyHostToDevice); 
// cudaThreadSynchronize();
 
 shared_mem_size=sizeof(int)*4000;
 
 Ham_dist<<<grid_dim3,block_dim3,shared_mem_size>>>(devData,devScore,count-100,s0,devqHam,devRind,devHdist,num); //have to change count
 checkCUDAError("kernel3 invocation");
 cudaThreadSynchronize();
 
 
 
 /*
 for(i=count-100;i<count;i++) {
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
 
 */
 
  
 cudaMemcpy(hostHdist,devHdist,sizeof(int)*n_blocks*num,cudaMemcpyDeviceToHost);
 //cudaThreadSynchronize();
 checkCUDAError("Hamm error");
 
//unc #pragma omp parallel for default(none) schedule(auto) shared(hostHdist, n_blocks, num) private(i)
 for(i=0;i<n_blocks*num;i++)
 if(!hostHdist[i])
 hostHdist[i]=B+2;
 //cudaThreadSynchronize();
 
 /*
 printf("\nHost Hamming Distance\n");
 
 for(i=0;i<n_blocks*num;i++)
 if(hostHdist[i]!=514)printf("%d ",hostHdist[i]);
 printf("\n");
 
 */
 
 
 /*printf("aftr hamdist updata\n"); 
 gettimeofday(&time,NULL);
 end=(time.tv_usec*1e-6)+time.tv_sec-start;
 printf("gpu scoreupdate time taken is %lf\n",end);


 gettimeofday(&time,NULL);
 start=(time.tv_usec*1e-6)+time.tv_sec;*/


 cudaMalloc((void**)&devResult,sizeof(int)*num);
  cudaMemcpy(devHdist,hostHdist,sizeof(int)*n_blocks*num,cudaMemcpyHostToDevice); 
 // cudaThreadSynchronize();



//*******************************************************************************************


Final_sort<<<grid_dim3,block_dim3>>>(devResult,devHdist,devRind,num);
checkCUDAError("kernel4 invocation");
cudaThreadSynchronize();


 result=(int*)malloc(sizeof(int)*num);
 cudaMemcpy(result,devResult,sizeof(int)*num,cudaMemcpyDeviceToHost);
// cudaThreadSynchronize();
 checkCUDAError("memcpy final sort"); 
 
 
 gettimeofday(&time,NULL);
 start=(time.tv_usec*1e-6)+time.tv_sec-start;
 printf("gpu sort taken is %lf\n",start);

/*
//printf("Final Sorted Result\n"); 
 for (i=0;i<num;i++)
 printf("%d ",result[i]);
 printf("\n");
 */
 
	cudaFree(devData);
	cudaFree(devScore);
	cudaFree(deVal);
	cudaFree(devErhash);
	cudaFree(devAlpha);
	cudaFree(devqHam);
	cudaFree(devHdist);
	cudaFree(devRind);
	cudaFree(devResult);
	cudaFree(devPtr);
	cudaFree(devOff);


}



//*******************************************************************************************
//*******************************************************************************************


__device__ int lsHash(int *n,int t,int *deVal)
{
 int i,j,k,b[512],c=0,s;
 //for(i=0;i<dB;i++)   b[i]=n>>i & 1;
 for(i=0;i<dNI;i++)
  for(j=0;j<32;j++) {
    b[c]=n[i]>>j & 1;
    c++;
  }
  
 //for(i=0;i<dB-1;i++) b[i]=b[i]^b[i+1];
 c=0;
 for(i=0;i<dK;i++)
 { 
  if(b[deVal[t*dK+i]]) {
    s=1;
    for(k=1;k<=i;k++) s*=2;
    c+=s;
  }
 }
 return(c);
}

void error_hash(int *q,int t,int *erhash,int *val,int *b)
{
 int i,j,k,co=0,c[K],count=1+K;
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
  if(b[val[t*K+i]]) { c[i]=1; erhash[0]+=pow(2,i); }
  else   c[i]=0;
 } 
 for(i=0;i<K;i++) {
  for(j=0;j<K;j++)
  {
   if(i==j)  erhash[i+1]+=(!c[j])*pow(2,j);
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

void checkCUDAError(const char *msg)
{
    cudaError_t err=cudaGetLastError();
    if(cudaSuccess!=err)
    {
        fprintf(stderr,"Cuda error: %s: %s.\n", msg, cudaGetErrorString( err));
        exit(EXIT_FAILURE);
    }
}


int ls_hash(int *n,int t,int val[N])
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
  if(b[val[t*K+i]]) 
    c=c+pow(2,i);
 }
 return(c);
}



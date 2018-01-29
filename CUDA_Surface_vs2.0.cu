#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<cuda.h>
#include<time.h>


//CUDA error wrapping

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/*
**************************************************************************************************
Author: Maxwel Gama Monteiro Junior
Contact: maxweljr@gmail.com

Description: Obtains the surface of a given crystal using its coordination number as order
parameter.
**************************************************************************************************
                       ;\
                      _' \_
                    ,' '  '`.
                   ;,)       \
                  /          :
                  (_         :
                   `--.       \			~ CUDA Crystal Surface Cutter vs2.0
                      /        `.

Version History:

1.0 - Does stuff

2.0 - Does stuff WITH MORE GLITTER! oh and also carries on counting on GPU only - see
deprecated function at bottom

**************************************************************************************************
**************************************************************************************************
**************************************************************************************************
*/

//Double approximated reciprocal square root function (drsqrt)

    __device__ __forceinline__ double drsqrt (double a)
    {
      double y, h, l, e;
      unsigned int ilo, ihi, g, f;
      int d;

      ihi = __double2hiint(a);
      ilo = __double2loint(a);
      if (((unsigned int)ihi) - 0x00100000U < 0x7fe00000U){
        f = ihi | 0x3fe00000;
        g = f & 0x3fffffff;
        d = g - ihi;
        a = __hiloint2double(g, ilo); 
        y = rsqrt (a);
        h = __dmul_rn (y, y);
        l = __fma_rn (y, y, -h);
        e = __fma_rn (l, -a, __fma_rn (h, -a, 1.0));
        // Round as shown in Peter Markstein, "IA-64 and Elementary Functions"
        y = __fma_rn (__fma_rn (0.375, e, 0.5), e * y, y);
        d = d >> 1;
        a = __hiloint2double(__double2hiint(y) + d, __double2loint(y));
      } else if (a == 0.0) {
        a = __hiloint2double ((ihi & 0x80000000) | 0x7ff00000, 0x00000000);
      } else if (a < 0.0) {
        a = __hiloint2double (0xfff80000, 0x00000000);
      } else if (isinf (a)) {
        a = __hiloint2double (ihi & 0x80000000, 0x00000000);
      } else if (isnan (a)) {
        a = a + a;
      } else {
        a = a * __hiloint2double (0x7fd00000, 0);
        y = rsqrt (a);
        h = __dmul_rn (y, y);
        l = __fma_rn (y, y, -h);
        e = __fma_rn (l, -a, __fma_rn (h, -a, 1.0));
        // Round as shown in Peter Markstein, "IA-64 and Elementary Functions"
        y = __fma_rn (__fma_rn (0.375, e, 0.5), e * y, y);
        a = __hiloint2double(__double2hiint(y) + 0x1ff00000,__double2loint(y));
      }
      return a;
    }


__global__ void parameter_counter(double *x_, double *y_, double *z_, int *count, uint atoms, double prmt_, uint *surfies_)
{
	int n = threadIdx.x + blockDim.x * blockIdx.x;

	extern __shared__ uint cache[];	
	uint temp = 0;
	

while (n < atoms)
{
	double xx_;
	double yy_;
	double zz_;
	double rij_;
	double x = x_[n];
	double y = y_[n];
	double z = z_[n];
	int kappa = 0; //save to thread register instead of reordering *count
	

		for (int neighbor = 0; neighbor < atoms; neighbor++)
		{
				xx_ =  x - x_[neighbor];
				yy_ = (y - y_[neighbor])*(y - y_[neighbor]);
				zz_ = (z - z_[neighbor])*(z - z_[neighbor]);
				rij_ = drsqrt(xx_*xx_ + yy_ + zz_);

			  if (1.0/rij_ < prmt_) 
			  {
				kappa++; //Count the atom itself to avoid branching 
			  }
			//if (1.0/rij_ < prmt_ && isinf(rij_) == 0)kappa++;//Check that one is not his own neighbor
		}

count[n] = kappa;

if(kappa < 13){
	temp++;
}  //Surface atoms will have less than the usual 12 neighbors of FCC (13 in this case because it is self-interacting)


n +=gridDim.x * blockDim.x;
}


cache[threadIdx.x] = temp;

__syncthreads();

//Perform sum reduction on temp values to obtain number of surface atoms

int u = blockDim.x/2;

while(u != 0)
	{
		if (threadIdx.x < u)
			{
			cache[threadIdx.x] += cache[threadIdx.x + u];
			}
	__syncthreads();
	u /= 2;
	}

if (threadIdx.x == 0) surfies_[blockIdx.x] = cache[0];



}



int main(void)
{

	int deviceCount;
	cudaGetDeviceCount (&deviceCount);
	if (deviceCount <1)
	{
		printf("CUDA supporting video card not detected. Go eat a sandwich or something.");
		return 0;
	}
	
	

	double *x, *y, *z, *dev_x, *dev_y, *dev_z;
	uint *surfies, *dev_surfies;
	int *count, *dev_count;

	int *lbl;
	double prmt;
	int natom;
	int i, j, p;
	int n_surface = 0;

	//Adjust the size of histogram accordingly to fit your maximum number of nearest neighbors
	int grp[14];

	//Change as suitable, this seems to work reasonably in general
	size_t block = 512;
	size_t thread = 512;

	size_t block_size = sizeof(uint)*block;

	


	FILE *finp, *fout, *fsupply;
	finp=fopen("coord_z.xyz","r");
	fsupply=fopen("num_edge.dat","w"); 
	fout=fopen("surface.xyz","w");
	
	fscanf(finp,"%d\n",&natom);
	fscanf(finp,"%lf\n",&prmt);
	

	cudaMallocHost((void**)&lbl, sizeof(int) * natom);
	cudaMallocHost((void**)&x,sizeof(double) * natom);
	cudaMallocHost((void**)&y,sizeof(double) * natom);
	cudaMallocHost((void**)&z,sizeof(double) * natom);
	cudaMallocHost((void**)&count,sizeof(int) * natom);
	cudaMallocHost((void**)&surfies, block_size);

	gpuErrchk( cudaMalloc((void**)&dev_surfies, block_size)	      );
	gpuErrchk( cudaMalloc((void**)&dev_x, sizeof(double) * natom) );
	gpuErrchk( cudaMalloc((void**)&dev_y, sizeof(double) * natom) );
	gpuErrchk( cudaMalloc((void**)&dev_z, sizeof(double) * natom) );
	gpuErrchk( cudaMalloc((void**)&dev_count, sizeof(int) * natom));

	j = 0;
	while (fscanf(finp,"%d %lf %lf %lf\n",&lbl[j], &x[j], &y[j], &z[j]) == 4)
	{

	j=j+1;

	}

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	cudaEventRecord(start,0);

	cudaMemcpyAsync(dev_x, x, sizeof(double)*natom, cudaMemcpyHostToDevice);
	cudaMemcpyAsync(dev_y, y, sizeof(double)*natom, cudaMemcpyHostToDevice);
	cudaMemcpyAsync(dev_z, z, sizeof(double)*natom, cudaMemcpyHostToDevice);
	
	for (p = 0; p < 14; p++)
	{
		grp[p] = 0;
	}

	prmt *= 0.8; //Lattice parameter/sqrt(2) ~ prmt * 0.7 is the nearest neighbor FCC distance, using 0.8 to make sure every neighbor is contained
		     //even if fast-math flags are used

	parameter_counter<<<thread, block, block_size>>>(dev_x,dev_y,dev_z,dev_count,natom, prmt, dev_surfies);
	
	cudaMemcpyAsync(count, dev_count, sizeof(int)*natom, cudaMemcpyDeviceToHost);
	cudaMemcpy(surfies, dev_surfies, block_size, cudaMemcpyDeviceToHost); //This is needed immediately by host


	for(p = 0; p < block; p++){
	n_surface+=surfies[p];
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop); 

	clock_t launch, finish;

	double cpu_time;
	launch = clock();


		fprintf(fout,"%d\n",n_surface);
		fprintf(fout,"%16.16lf\n",prmt/0.8);	

	for (i = 0; i < natom; i++) //Writing and wrap-up done on host
		{
			j = count[i];
			grp[j]++;
			
			if(count[i] < 13){

			fprintf(fout,"%d \t %16.15lf \t %16.15lf \t %16.15lf \n",lbl[i], x[i], y[i], z[i]);
			}
		}

	
	j = 0;
	for (i = 0; i < 14; i++)
		{
			fprintf(fsupply, "%d\n", grp[i] );
			j+=grp[i];
		}
		
			fprintf(fsupply,"%d\n",j);



	cudaFreeHost(x);
	cudaFreeHost(y);
	cudaFreeHost(z);
	cudaFreeHost(count);
	cudaFreeHost(surfies);
	cudaFreeHost(lbl);
	cudaFree(dev_x);
	cudaFree(dev_y);
	cudaFree(dev_z);
	cudaFree(dev_count);
	cudaFree(dev_surfies);

finish = clock();
cpu_time = ((double)(finish - launch)) / CLOCKS_PER_SEC;

fclose(finp);
fclose(fout);
fclose(fsupply);


printf("\n\nCPU process finished at %16.8lf seconds\n\n", cpu_time);
printf("\n\nGPU process finished at %.8f seconds\n\n", elapsedTime/1000);
printf("======================================================================~\n");


return 0;
}


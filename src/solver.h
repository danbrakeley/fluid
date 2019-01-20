/*
 * solver
 */


 /* macros */

#define nullptr 0
#define IX(i, j) ((i) + (N + 2) * (j))

/* FluidBox */

typedef struct {
	int InN, InDim;
	int OutN, OutDim;
	float dt;
	float Diffusion;
	float Viscosity;
} FluidBox;

FluidBox* CreateFluidBox(int n, float dt, float diffusion, float viscosity);
void DestroyFluidBox(FluidBox** box);
void StepDensity(FluidBox* box, float* x, float* x0, float* u, float* v);
void StepVelocity(FluidBox* box, float* u, float* v, float* u0, float* v0);

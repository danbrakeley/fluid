#include "solver.h"
#include <stdlib.h>

#define SWAP(x0, x)   \
	{                   \
		float* tmp = x0;  \
		x0         = x;   \
		x          = tmp; \
	}

FluidBox* CreateFluidBox(int n, float dt, float diffusion, float viscosity) {
	FluidBox* f = (FluidBox*)malloc(sizeof(FluidBox));

	f->InN       = n;
	f->InDim     = f->InN * f->InN;
	f->OutN      = n + 2;
	f->OutDim    = f->OutN * f->OutN;
	f->dt        = dt;
	f->Diffusion = diffusion;
	f->Viscosity = viscosity;

	return f;
}

void DestroyFluidBox(FluidBox** box) {
	if (box != nullptr) {
		free(*box);
		*box = nullptr;
	}
}

void addSource(FluidBox* box, float* x, float* s, float dt) {
	for (int i = 0; i < box->OutDim; i++) {
		x[i] += dt * s[i];
	}
}

void setBoundary(FluidBox* box, int b, float* x) {
	const int N = box->InN;

	// edges
	for (int i = 1; i <= N; i++) {
		x[IX(i, 0)]     = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	for (int j = 1; j <= N; j++) {
		x[IX(0, j)]     = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
		x[IX(N + 1, j)] = b == 1 ? -x[IX(N, j)] : x[IX(N, j)];
	}

	// corners
	x[IX(0, 0)]         = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)]     = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)]     = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void linearSolver(FluidBox* box, int b, float* x, float* x0, float a, float c) {
	const int N = box->InN;

	for (int k = 0; k < 20; k++) {
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
			}
		}
		setBoundary(box, b, x);
	}
}

void diffuse(FluidBox* box, int b, float* x, float* x0, float diff, float dt) {
	float a = dt * diff * box->InDim;
	linearSolver(box, b, x, x0, a, 1 + 4 * a);
}

void advect(FluidBox* box, int b, float* d, float* d0, float* u, float* v, float dt) {
	const int N = box->InN;
	float dt0   = dt * N;
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			float x = i - dt0 * u[IX(i, j)];
			float y = j - dt0 * v[IX(i, j)];

			if (x < 0.5f) {
				x = 0.5f;
			}
			if (x > N + 0.5f) {
				x = N + 0.5f;
			}

			int i0 = (int)x;
			int i1 = i0 + 1;

			if (y < 0.5f) {
				y = 0.5f;
			}
			if (y > N + 0.5f) {
				y = N + 0.5f;
			}

			int j0   = (int)y;
			int j1   = j0 + 1;
			float s1 = x - i0;
			float s0 = 1 - s1;
			float t1 = y - j0;
			float t0 = 1 - t1;

			d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		}
	}
	setBoundary(box, b, d);
}

void project(FluidBox* box, float* u, float* v, float* p, float* div) {
	const int N = box->InN;
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
			p[IX(i, j)]   = 0;
		}
	}
	setBoundary(box, 0, div);
	setBoundary(box, 0, p);

	linearSolver(box, 0, p, div, 1, 4);

	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
			v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
		}
	}
	setBoundary(box, 1, u);
	setBoundary(box, 2, v);
}

void StepDensity(FluidBox* box, float* x, float* x0, float* u, float* v) {
	addSource(box, x, x0, box->dt);
	SWAP(x0, x);
	diffuse(box, 0, x, x0, box->Diffusion, box->dt);
	SWAP(x0, x);
	advect(box, 0, x, x0, u, v, box->dt);
}

void StepVelocity(FluidBox* box, float* u, float* v, float* u0, float* v0) {
	addSource(box, u, u0, box->dt);
	addSource(box, v, v0, box->dt);
	SWAP(u0, u);
	diffuse(box, 1, u, u0, box->Viscosity, box->dt);
	SWAP(v0, v);
	diffuse(box, 2, v, v0, box->Viscosity, box->dt);
	project(box, u, v, u0, v0);
	SWAP(u0, u);
	SWAP(v0, v);
	advect(box, 1, u, u0, u0, v0, box->dt);
	advect(box, 2, v, v0, u0, v0, box->dt);
	project(box, u, v, u0, v0);
}

/*
	======================================================================
	 demo.c --- protoype to show off the simple solver
	----------------------------------------------------------------------
	 Author : Jos Stam (jstam@aw.sgi.com)
	 Creation Date : Jan 9 2003

	 Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

	=======================================================================
*/

#include "solver.h"
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>

/* global variables */

static FluidBox* box = nullptr;
static float force, source, diffusionStep, viscosityStep;
static int bDrawVelocity;

static float *u, *v, *u_prev, *v_prev;
static float *dens, *dens_prev;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;

static void printStatus() {
	printf("\n\n\n\n\n");
	printf("\n\n\n\n\n");
	printf("\n\n\n\n\n");
	printf("\n\n\n\n\n");
	printf(" right mouse  Add densities\n");
	printf("  left mouse  Add velocities\n");
	printf("           v  Draw velocity (currently %s)\n", bDrawVelocity ? "on" : "off");
	printf("           c  Clear density/velocity maps\n");
	printf("           q  Quit\n");
	printf("\n");
	printf("    r/f, t/g  Diffusion: %f (step %f)\n", box->Diffusion, diffusionStep);
	printf("    y/h, u/j  Viscosity: %f (step %f)\n", box->Viscosity, viscosityStep);
	printf("         o/l  Time Step: %f\n", box->dt);
	printf("\n");
	printf("  Velocity Force: %f\n", force);
	printf("  Density Source: %f\n", source);
}

/*
	----------------------------------------------------------------------
	 free/clear/allocate simulation data
	----------------------------------------------------------------------
*/

static void free_data() {
	if (u) {
		free(u);
		u = nullptr;
	}
	if (v) {
		free(v);
		v = nullptr;
	}
	if (u_prev) {
		free(u_prev);
		u_prev = nullptr;
	}
	if (v_prev) {
		free(v_prev);
		v_prev = nullptr;
	}
	if (dens) {
		free(dens);
		dens = nullptr;
	}
	if (dens_prev) {
		free(dens_prev);
		dens_prev = nullptr;
	}
}

static void clear_data(FluidBox* box) {
	for (int i = 0; i < box->OutDim; i++) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
}

static int allocate_data(FluidBox* box) {
	u         = (float*)malloc(box->OutDim * sizeof(float));
	v         = (float*)malloc(box->OutDim * sizeof(float));
	u_prev    = (float*)malloc(box->OutDim * sizeof(float));
	v_prev    = (float*)malloc(box->OutDim * sizeof(float));
	dens      = (float*)malloc(box->OutDim * sizeof(float));
	dens_prev = (float*)malloc(box->OutDim * sizeof(float));

	if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev) {
		fprintf(stderr, "cannot allocate data\n");
		return 0;
	}
	return 1;
}

/*
	----------------------------------------------------------------------
	 OpenGL specific drawing routines
	----------------------------------------------------------------------
*/

static void pre_display() {
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

static void post_display() {
	glutSwapBuffers();
}

static void draw_velocity() {
	const int N = box->InN;
	int i, j;
	float h, x, y;

	glLineWidth(1.0f);
	glBegin(GL_LINES);

	glColor3f(1.0f, 0.0f, 0.0f);
	h = 1.0f / (float)N;
	for (i = 1; i <= N; i++) {
		x = ((float)i - 0.5f) * h;
		for (j = 1; j <= N; j++) {
			y = ((float)j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + u[IX(i, j)], y + v[IX(i, j)]);
		}
	}

	glEnd();
}

static void draw_density() {
	const int N = box->InN;
	int i, j;
	float h, x, y, d00, d01, d10, d11;

	glBegin(GL_QUADS);

	h = 1.0f / N;
	for (i = 0; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 0; j <= N; j++) {
			y = (j - 0.5f) * h;

			d00 = dens[IX(i, j)];
			d01 = dens[IX(i, j + 1)];
			d10 = dens[IX(i + 1, j)];
			d11 = dens[IX(i + 1, j + 1)];

			glColor3f(d00, d00, d00);
			glVertex2f(x, y);
			glColor3f(d10, d10, d10);
			glVertex2f(x + h, y);
			glColor3f(d11, d11, d11);
			glVertex2f(x + h, y + h);
			glColor3f(d01, d01, d01);
			glVertex2f(x, y + h);
		}
	}

	glEnd();
}

/*
	----------------------------------------------------------------------
	 relates mouse movements to forces sources
	----------------------------------------------------------------------
*/

static void get_from_UI(FluidBox* box, float* d, float* u, float* v) {
	const int N = box->InN;
	int i, j;
	for (i = 0; i < box->OutDim; i++) {
		u[i] = v[i] = d[i] = 0.0f;
	}

	if (!mouse_down[0] && !mouse_down[2]) {
		return;
	}

	i = (int)((mx / (float)win_x) * N + 1);
	j = (int)(((win_y - my) / (float)win_y) * N + 1);

	if (i < 1 || i > N || j < 1 || j > N) {
		return;
	}

	if (mouse_down[0]) {
		u[IX(i, j)] = force * (mx - omx);
		v[IX(i, j)] = force * (omy - my);
	}

	if (mouse_down[2]) {
		d[IX(i, j)] = source;
	}

	omx = mx;
	omy = my;

	return;
}

/*
	----------------------------------------------------------------------
	 GLUT callback routines
	----------------------------------------------------------------------
*/

static void key_func(unsigned char key, int x, int y) {
	switch (key) {
	case 'c':
		clear_data(box);
		break;
	case 'q':
		free_data();
		exit(0);
		break;
	case 'v':
		bDrawVelocity = !bDrawVelocity;
		break;
	case 'r':
		box->Diffusion += diffusionStep;
		break;
	case 'f':
		box->Diffusion -= diffusionStep;
		break;
	case 't':
		diffusionStep *= 10.0f;
		break;
	case 'g':
		diffusionStep /= 10.0f;
		break;
	case 'y':
		box->Viscosity += viscosityStep;
		break;
	case 'h':
		box->Viscosity -= viscosityStep;
		break;
	case 'u':
		viscosityStep *= 10.0f;
		break;
	case 'j':
		viscosityStep /= 10.0f;
		break;
	case 'o':
		box->dt += 0.05f;
		break;
	case 'l':
		box->dt -= 0.05f;
		break;
	}

	printStatus();
}

static void mouse_func(int button, int state, int x, int y) {
	omx = mx = x;
	omx = my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func(int x, int y) {
	mx = x;
	my = y;
}

static void reshape_func(int width, int height) {
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}

static void idle_func(void) {
	get_from_UI(box, dens_prev, u_prev, v_prev);
	StepVelocity(box, u, v, u_prev, v_prev);
	StepDensity(box, dens, dens_prev, u, v);

	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func(void) {
	pre_display();

	draw_density();
	if (bDrawVelocity) {
		draw_velocity();
	}

	post_display();
}

/*
	----------------------------------------------------------------------
	 open_glut_window --- open a glut compatible window and set callbacks
	----------------------------------------------------------------------
*/

static void open_glut_window() {
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(100, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("Alias | wavefront");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	pre_display();

	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
}

/*
	----------------------------------------------------------------------
	 main --- main routine
	----------------------------------------------------------------------
*/

int main(int argc, char** argv) {
	// starting values
	int N         = 64;       // grid resolution
	float dt      = 0.1f;     // time step
	float diff    = 0.0f;     // diffusion rate of the density
	diffusionStep = 0.0001f;  // rate at which to raise/lower the diffusion
	float visc    = 0.0f;     // viscosity of the fluid
	viscosityStep = 0.0001f;  // rate at which to raise/lower the viscosity
	force         = 5.0f;     // scales the mouse movement that generate a force
	source        = 100.0f;   // amount of density that will be deposited
	bDrawVelocity = 1;

	box = CreateFluidBox(64, dt, diff, visc);

	glutInit(&argc, argv);

	printStatus();

	if (!allocate_data(box)) {
		exit(1);
	}
	clear_data(box);

	win_x = 1024;
	win_y = 1024;
	open_glut_window();

	glutMainLoop();

	free_data();
	DestroyFluidBox(&box);
	return 0;
}

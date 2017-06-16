// Mars lander simulator
// Version 1.9
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2016

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"
#include <math.h>

enum pilotMode { DESCENT = 0, LAUNCH = 1 };
pilotMode mode;

//orbit parameters
vector3d targetPlane;
double periapsis;
double apoapsis;

//variables to calculate current trajectory
double mechEnergy; //mechanical energy of the spacecraft
double coefB; //coefficient of the linear term in the quadratic equation whose solution is the min distance of the spacecraft from the centre of Mars
double constC; //constant term in the quadratic equation
double minRad; //min distance of the spacecraft trajectory from the centre of Mars
double maxRad; //max distance of the spacecraft trajectory from the centre of Mars

//control constants for orbital injection
const double vn0 = -5.0;
const double vtg0 = -0.1;
const double vr0 = 0.01;
const double k1 = 1E-5;
const double k2 = 1E-6;
const double k3 = 1E-3;
const double k4 = 9E-4;
//control variables
double vmax;
vector3d posIntegral;
vector3d throttleVect;
vector3d rUnit;

vector3d getGravity(vector3d pos, double mass, double centreMass);
vector3d getDrag(vector3d velocity, double area, double airDensity, double coeff);
double landerMass;

static vector3d prevPosition = position - velocity * delta_t;
static vector3d vertUnit;
static vector3d prevOrientation;	//to introduce angular velocity of the spacecraft itself
static vector3d nextOrientation;

static double relRot[16];
/*static double oMat[16];				//orientation matrix of the spacecraft
static double prevOMat[16];			//previous orientation matrix
static double invPrevOMat[16];		//inverse of the previous orientation matrix
//static double nextOMat[16];			//next orientation matrix
static double rotMat[16];			//rotation matrix
static double rotMatT[16];
static vector3d omega;*/

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{

	if (mode == DESCENT) {
		//stabilize attitude, opposite to velocity
		attitude_stabilization(-velocity);
		//version 1.
		float height	= position.abs() - MARS_RADIUS;
		float kh		= 0.045;
		float gain		= 1;
		float error		= 0;
		float delta		= 0.5;
		vector3d radUnit = position.norm();
		if (parachute_status == DEPLOYED) {
			delta = 0.2;
			kh = 0.1;
		}
		error = -(0.2 + kh*height + radUnit * velocity);
		if (height < EXOSPHERE/1.8) {
			throttle = delta + gain * error; //the simulation program will take care of the min/max thing
		} else {
			throttle = 0;
		}
		//throttle = delta + gain * error; //the simulation program will take care of the min/max thing

		//the parachute is useful after all
		if (parachute_status == NOT_DEPLOYED && safe_to_deploy_parachute() && height < EXOSPHERE/2) {
			//would-be acceleration:
			parachute_status = DEPLOYED;
		}

		//orbital reentry: slow the spacecraft does not intersect with the exosphere
		mechEnergy = pow(velocity.abs(), 2.0)/2 - GRAVITY * MARS_MASS/position.abs();
		coefB = GRAVITY * MARS_MASS / mechEnergy;
		constC = -(position ^ velocity).abs2()/(2*mechEnergy);
		minRad = coefB + pow(pow(coefB, 2) - 4*constC, 0.5);
		minRad = -minRad/2;
		if (minRad > MARS_RADIUS) { //it could be MARS_RADIUS + const*EXOSPHERE, but it would need negligibly less fuel and would take much more time
			throttle = 1;
		}
	} else if (mode == LAUNCH) {
		landerMass = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY;
		rUnit = ((targetPlane ^ position) ^ targetPlane); // unit radial vector
		if (rUnit.abs() < SMALL_NUM) {
			rUnit.x = -targetPlane.z; rUnit.y = 0.0; rUnit.z = targetPlane.x;
			if (rUnit.abs() < SMALL_NUM) {rUnit.x = 0.0; rUnit.y = targetPlane.z; rUnit.z = -targetPlane.y;}
		}
		rUnit = rUnit.norm();

		//calculating the max velocity from the orbit parameters
		vmax = pow(2.0*GRAVITY*MARS_MASS*apoapsis/(periapsis*(periapsis + apoapsis)), 0.5);

		/*//using PID control
		posIntegral += (position - rUnit*minRadius)*delta_t;
		if (posIntegral.abs() > 15) {
			posIntegral = posIntegral.norm()*15;
		}
		throttleVect = ((position - rUnit*minRadius)*k1 + posIntegral*k2 + velocity*k3 )*landerMass + getGravity(position, landerMass, MARS_MASS)*k4;
		printf("%f\n", position.abs()-minRadius);

		attitude_stabilization(-throttleVect);
		throttle = throttleVect.abs();*/

		//assumes that the launch is performed in the target orbital plane

		//calculates current trajectory
		mechEnergy = pow(velocity.abs(), 2.0)/2 - GRAVITY * MARS_MASS/position.abs();
		coefB = GRAVITY * MARS_MASS / mechEnergy;
		constC = -(position ^ velocity).abs2()/(2*mechEnergy);
		minRad = coefB + pow(pow(coefB, 2) - 4*constC, 0.5);
		minRad = -minRad/2;
		maxRad = coefB - pow(pow(coefB, 2) - 4*constC, 0.5);
		maxRad = -maxRad/2;

	}

}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{

	landerMass = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY;
	double landerArea = pow(LANDER_SIZE, 2.0) * M_PI;
	double atmDensity = atmospheric_density(position);
	vector3d totalFotrce = getGravity(position, landerMass, MARS_MASS) + getDrag(velocity, landerArea, atmDensity, DRAG_COEF_LANDER) + thrust_wrt_world();
	if (parachute_status == DEPLOYED) {
		double chuteArea = 20.0 * pow(LANDER_SIZE, 2.0);
		totalFotrce += getDrag(velocity, chuteArea, atmDensity, DRAG_COEF_CHUTE);
	}
	vector3d acceleration = totalFotrce / landerMass;

	vector3d nextPosition = 2.0 * position - prevPosition + acceleration * pow(delta_t, 2.0);

	prevPosition = position;
	position = nextPosition;

	velocity = velocity + acceleration * delta_t;

	// Here we can apply an autopilot to adjust the thrust, parachute and attitude
	if (autopilot_enabled) autopilot();

	// Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
	//if (stabilized_attitude) attitude_stabilization();
	if (stabilized_attitude) {
		xyz_euler_to_matrix(relativeAttitude, relRot);
		attitude_stabilization(matDotVect(position.norm(), relRot));
	}
	if (orientLocked) {
	}
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
	// The parameters to set are:
	// position - in Cartesian planetary coordinate system (m)
	// velocity - in Cartesian planetary coordinate system (m/s)
	// orientation - in lander coordinate system (xyz Euler angles, degrees)
	// delta_t - the simulation time step
	// boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
	// scenario_description - a descriptive string for the help screen

	scenario_description[0] = "circular orbit";
	scenario_description[1] = "descent from 10km";
	scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
	scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
	scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
	scenario_description[5] = "descent from 200km";
	scenario_description[6] = "areostationary orbit";
	scenario_description[7] = "";
	scenario_description[8] = "";
	scenario_description[9] = "";

	orientLocked = false;

	switch (scenario) {

	case 0:
		// a circular equatorial orbit
		/*position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		velocity = vector3d(0.0, 3247.087385863725, 0.0);*/
		velocity = vector3d(-3247.087385863725, 0.0, 0.0);
		position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
		orientation = vector3d(0.0, 0.0, 0.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	case 1:
		// a descent from rest at 10km altitude
		position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
		velocity = vector3d(0.0, 0.0, 0.0);
		orientation = vector3d(0.0, 0.0, 90.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = true;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	case 2:
		// an elliptical polar orbit
		position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
		velocity = vector3d(3500.0, 0.0, 0.0);
		orientation = vector3d(0.0, 0.0, 90.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	case 3:
		// polar surface launch at escape velocity (but drag prevents escape)
		/*position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
		velocity = vector3d(0.0, 0.0, 5027.0);
		orientation = vector3d(0.0, 0.0, 0.0);*/
		position = vector3d(MARS_RADIUS + LANDER_SIZE/2.0, 0.0, 0.0);
		velocity = vector3d(5027.0, 0.0, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	case 4:
		// an elliptical orbit that clips the atmosphere each time round, losing energy
		position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
		velocity = vector3d(4000.0, 0.0, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	case 5:
		// a descent from rest at the edge of the exosphere
		position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
		velocity = vector3d(0.0, 0.0, 0.0);
		orientation = vector3d(0.0, 0.0, 90.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = true;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	case 6:
		//areostationary orbit
		float orbitRad;
		orbitRad = pow(GRAVITY*MARS_MASS*pow(MARS_DAY, 2)/pow(2*M_PI, 2), (double)1/3);
		position = vector3d(orbitRad, 0.0, 0.0);
		velocity = vector3d(0.0, 2*M_PI*orbitRad/MARS_DAY, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	case 7:
		//injection into circular equatorial orbit
		position = vector3d(MARS_RADIUS + LANDER_SIZE/2.0, 0.0, 0.0);
		velocity = vector3d(1.0, 0.0, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		periapsis = 1.2*MARS_RADIUS;
		apoapsis = 1.2*MARS_RADIUS;
		targetPlane = vector3d(0.0, 0.0, 1.0);
		mode = LAUNCH;
		break;

	case 8:
		//injection into circular equatorial orbit
		position = vector3d(MARS_RADIUS + LANDER_SIZE/2.0, 0.0, 0.0);
		velocity = vector3d(0.0, 0.0, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		periapsis = 1.2*MARS_RADIUS;
		apoapsis = 1.2*MARS_RADIUS;
		targetPlane = vector3d(0.0, 0.0, 1.0);
		mode = LAUNCH;
		break;

	case 9:
		velocity = vector3d(0.0, -3247.087385863725, 0.0);
		position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
		orientation = vector3d(0.0, 0.0, 0.0);
		prevOrientation = orientation;
		delta_t = 0.1;
		prevPosition = position - velocity * delta_t;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = false;
		mode = DESCENT;
		break;

	}


	targetPlane = targetPlane.norm();
	throttleVect = vector3d(0.0, 0.0, 0.0);
}

vector3d getGravity(vector3d pos, double mass, double centreMass) {
	vector3d force;
	force = -pos.norm();
	force *= mass * centreMass * GRAVITY;
	force /= pos.abs2();
	return force;
}

vector3d getDrag(vector3d velocity, double area, double airDensity, double coeff) {
	vector3d force;
	force = -velocity.norm();
	force *= velocity.abs2();
	force *= area * airDensity * coeff;
	force /= 2.0;
	return force;
}

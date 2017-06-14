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
double minRadius;
double maxRadius;

//control constants for orbital injection
const double vn0 = -1.0;
const double vtg0 = 1;
const double vr0 = 1;
const double vtgSlope = 1;
const double vrSlope = 1;
const double kn = 1;
const double kr = 1;
const double ktg = 1;
//control variables
/*double vn;
double vtg;
double vr;
double vnTarget;
double vtgTarget;
double vrTarget;*/
double vmax;
vector3d thrustVect;
vector3d rPos;
vector3d nPos;

vector3d getGravity(vector3d pos, double mass, double centreMass);
vector3d getDrag(vector3d velocity, double area, double airDensity, double coeff);
double landerMass;

static vector3d prevPosition = position - velocity * delta_t;
static vector3d vertUnit;
static vector3d prevOrientation;	//to introduce angular velocity of the spacecraft itself
static vector3d nextOrientation;
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
		double orbitEnergy;
		double coefB; //coefficient of the linear term in the quadratic equation whose solution is the min distance of the spacecraft from the centre of Mars
		double constC; //constant term in the quadratic equation
		double minRad; //min distance of the spacecraft from the centre of Mars
		orbitEnergy = pow(velocity.abs(), 2.0)/2 - GRAVITY * MARS_MASS/position.abs();
		coefB = GRAVITY * MARS_MASS / orbitEnergy;
		constC = -(position ^ velocity).abs2()/(2*orbitEnergy);
		minRad = coefB + pow(pow(coefB, 2) - 4*constC, 0.5);
		minRad = -minRad/2;
		if (minRad > MARS_RADIUS) { //it could be MARS_RADIUS + const*EXOSPHERE, but it would need negligibly less fuel and would take much more time
			throttle = 1;
		}
	} else if (mode == LAUNCH) {
		landerMass = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY;
		rPos = (targetPlane ^ position) ^ targetPlane; // component of the position in the orbit plane
		nPos = (targetPlane * position) * targetPlane; // component of the position normal to the orbit plane (distance from orbit plane)

		//calculating the max velocity from the orbit parameters
		vmax = 2.0*GRAVITY*MARS_MASS*maxRadius/(minRadius*(minRadius + maxRadius));

		thrustVect = targetPlane*kn * ( nPos.abs()*vn0 - (targetPlane*velocity) ) //normal component of thrust
					+ kr * ((minRadius*rPos.norm() - rPos)*vr0 - (rPos.norm()*velocity)*rPos.norm()) //radial component of thrust
					+ ktg*(targetPlane ^ rPos).norm() * ( (vmax - velocity*(targetPlane ^ rPos).norm()) ); //tangential component of thrust
		thrustVect = thrustVect*landerMass;

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

	//the spacecraft has moment of inertia
	//TODO does not work
	/*xyz_euler_to_matrix(orientation, oMat);
	xyz_euler_to_matrix(prevOrientation, prevOMat);
	invert(prevOMat, invPrevOMat);
	dotMat(oMat, invPrevOMat, rotMat);
	dotMat(rotMat, oMat, nextOMat);
	prevOrientation = orientation;
	orientation = matrix_to_xyz_euler(nextOMat);*/
	/*xyz_euler_to_matrix(orientation, oMat);
	xyz_euler_to_matrix(prevOrientation, prevOMat);
	invert(prevOMat, invPrevOMat);
	dotMat(oMat, invPrevOMat, rotMat);
	double theta = acos((rotMat[0] + rotMat[5] + rotMat[10] - 1)/2);
	transpose(rotMat, rotMatT);
	for (int ii = 0; ii < 16; ii ++) {
		rotMat[ii] -= rotMatT[ii];
		rotMat[ii] *= theta;
		rotMat[ii] /= sin(theta)*2*delta_t;
	}
	omega = vector3d(rotMat[6], rotMat[1], rotMat[12]);

	if (omega.abs() > SMALL_NUM) {
		rotateOrientation(omega.abs()*delta_t, omega.norm());
	}*/

	// Here we can apply an autopilot to adjust the thrust, parachute and attitude
	if (autopilot_enabled) autopilot();

	// Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
	//if (stabilized_attitude) attitude_stabilization();
	if (stabilized_attitude) {
		attitude_stabilization();
		/*vector3d axis = -(position^prevPosition).norm();
		double dPhi = acos(position.norm()*prevPosition.norm());
		//printf("%f, %f, %f; %f\n", axis.x, axis.y, axis.z, dPhi);
		rotateOrientation(dPhi, axis);*/
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
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    prevOrientation = orientation;
    delta_t = 0.1;
    prevPosition = position - velocity * delta_t;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    mode = LAUNCH;
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

  case 8:

    break;

  case 9:
    break;

  }


	targetPlane = targetPlane.norm();
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

#include <dart/dart.hpp>
#include <dart/utils/utils.hpp>
#include <dart/gui/gui.hpp>
#include <dart/dart.hpp>
#include "BVHData.h"
#include <iostream>
#include <vector>
#include <map>

BVHData bvh;
FILE* fqpd=fopen("qpd.txt","r+");
FILE* fdqpd=fopen("dqpd.txt","r+");
FILE* fddqpd=fopen("ddqpd.txt","r+");
FILE* ftpd=fopen("tpd.txt","r+");
FILE* fmass=fopen("mass.txt","r+");
FILE* fcg=fopen("fcg.txt","r+");
FILE* fcons=fopen("fcons.txt","r+");
const double default_speed_increment = 0.5;

const int default_ik_iterations = 4500;

const double default_force =  50.0; // N
const int default_countdown = 100;  // Number of timesteps for applying force

using namespace dart::common;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::gui;
using namespace dart::gui::glut;
using namespace dart::utils;
using namespace dart::math;

class Controller
{
public:
    /// Constructor
    Controller(const SkeletonPtr& biped)
    : mBiped(biped),
    mSpeed(0.0)
    {
        int nDofs = mBiped->getNumDofs();
        mBiped->setTimeStep(0.01); //here
        
        mForces = Eigen::VectorXd::Zero(nDofs);
        
        mKp = Eigen::MatrixXd::Identity(nDofs, nDofs);
        mKd = Eigen::MatrixXd::Identity(nDofs, nDofs);
        
        for(std::size_t i = 0; i < 6; ++i)      
        {
            mKp(i, i) = 0.0;
            mKd(i, i) = 0.0;
        }
        
        for(std::size_t i = 6; i < mBiped->getNumDofs(); ++i)
        {
            mKp(i, i) = 100; //1000
            mKd(i, i) = 20; //50
        }
        setTargetPositions(mBiped->getPositions());
        //mBiped->setPositions(bvh.frame[0]);
    }
    
    /// Reset the desired dof position to the current position
    void setTargetPositions(const Eigen::VectorXd& pose)
    {
        mTargetPositions = pose;
    }
    
    /// Clear commanding forces
    void clearForces()
    {
        mForces.setZero();
    }
    
    /// Add commanding forces from PD controllers (Lesson 2 Answer)
    void addPDForces(int n)
    {
        Eigen::VectorXd q = mBiped->getPositions();
        for(int i=0;i<57;i++){
            std::cout<<q[i]<<", ";
        }
        std::cout<<std::endl;
        Eigen::VectorXd dq = mBiped->getVelocities();
        
        //Eigen::VectorXd p = -mKp * (q - mTargetPositions);
        Eigen::VectorXd p = -mKp * (q - bvh.frame[n]);
        Eigen::VectorXd d = -mKd * dq;
        
        mForces += p + d;
        mBiped->setForces(mForces);
    }
    
    /// Add commanind forces from Stable-PD controllers (Lesson 3 Answer)
    void addSPDForces(int n)
    {
        Eigen::VectorXd q = mBiped->getPositions();
        for(int i=0;i<57;i++){
            char s [100];
            sprintf(s,"%f\t",q.row(i)[0]);
            fputs(s,fqpd);
        }
        fputs("\n",fqpd);
        Eigen::VectorXd dq = mBiped->getVelocities();
        for(int i=0;i<57;i++){
            char s [100];
            sprintf(s,"%f\t",dq.row(i)[0]);
            fputs(s,fdqpd);
        }
        fputs("\n",fdqpd);
        //get massmatrix
        Eigen::MatrixXd mass= mBiped->getMassMatrix();
        for(int i=0;i<57;i++){
            for(int j=0;j<57;j++){
                char s [100];
                sprintf(s,"%f\t",mass.row(i)[j]);
                fputs(s,fmass);
            }
            fputs("\n",fmass);
        }
        
        Eigen::MatrixXd invM = (mBiped->getMassMatrix()
                                + mKd * mBiped->getTimeStep()).inverse();
        Eigen::VectorXd p =
        -mKp * (q + dq * mBiped->getTimeStep() - bvh.frame[n]);
        Eigen::VectorXd d = -mKd * dq;
        
        Eigen::VectorXd qddot =invM * (p + d);
        for(int i=0;i<57;i++){
            char s [100];
            sprintf(s,"%f\t",qddot.row(i)[0]);
            fputs(s,fddqpd);
        }
        fputs("\n",fddqpd);
        
        mForces += p + d - mKd * qddot * mBiped->getTimeStep();
        mBiped->setForces(mForces);
        for(int i=0;i<57;i++){
            char s [100];
            sprintf(s,"%f\t",mForces.row(i)[0]);
            fputs(s,ftpd);
        }
        fputs("\n",ftpd);
        
    }
    
protected:
    /// The biped Skeleton that we will be controlling
    SkeletonPtr mBiped;
    
    /// Joint forces for the biped (output of the Controller)
    Eigen::VectorXd mForces;
    
    /// Control gains for the proportional error terms in the PD controller
    Eigen::MatrixXd mKp;
    
    /// Control gains for the derivative error terms in the PD controller
    Eigen::MatrixXd mKd;
    
    /// Target positions for the PD controllers
    Eigen::VectorXd mTargetPositions;
    
    /// For velocity actuator: Current speed of the skateboard
    double mSpeed;
};


class MyWindow : public SimWindow
{
public:
    /// Constructor
    MyWindow(const WorldPtr& world)
    : mForceCountDown(0),
    mPositiveSign(true)
    {
        setWorld(world);
        
        mController = dart::common::make_unique<Controller>(
                                                            mWorld->getSkeleton("biped"));
    }
   
    int index = 0;
    /// Handle keyboard input
    void keyboard(unsigned char key, int x, int y) override
    {
        switch(key)
        {
            case ',':
                mForceCountDown = default_countdown;
                mPositiveSign = false;
                break;
            case '.':
                mForceCountDown = default_countdown;
                mPositiveSign = true;
                break;
	    case 'q':
		++index;
		bvh.setPositionAt(index);
		break;		
	    case 'e':
		--index;
		bvh.setPositionAt(index);
		break;		
            default:
                SimWindow::keyboard(key, x, y);
        }
    }
    int number=0;
    void timeStepping() override
    {
        mController->clearForces();
        
        // Lesson 3
        mController->addSPDForces(number);
        number++;
        
        // Step the simulation forward
        SimWindow::timeStepping();
    }
    
protected:
    std::unique_ptr<Controller> mController;
    
    /// Number of iterations before clearing a force entry
    int mForceCountDown;
    
    /// Whether a force should be applied in the positive or negative direction
    bool mPositiveSign;
    
};

// Load a biped model and enable joint limits and self-collision
// (Lesson 1 Answer)
SkeletonPtr loadBiped()
{
    // Create the world with a skeleton
    //WorldPtr world = SkelParser::readWorld("dart://sample/skel/biped.skel");
    //assert(world != nullptr);
    bvh.loadBVH("OptiTrack-IITSEC2007.bvh");
    SkeletonPtr biped = bvh.skeleton;
    //std::cout<<bvh.frame;
    
    // Set joint limits
    for(std::size_t i = 0; i < biped->getNumJoints(); ++i)
        biped->getJoint(i)->setPositionLimitEnforced(true);
    
    // Enable self collision check but ignore adjacent bodies
    //biped->enableSelfCollisionCheck();
    biped->disableAdjacentBodyCheck();
    
    return biped;
}


// Set the actuator type for four wheel joints to "VELOCITY" (Lesson 6 Answer)
void setVelocityAccuators(SkeletonPtr biped)
{
    Joint* wheel1 = biped->getJoint("joint_front_left");
    Joint* wheel2 = biped->getJoint("joint_front_right");
    Joint* wheel3 = biped->getJoint("joint_back_left");
    Joint* wheel4 = biped->getJoint("joint_back_right");
    wheel1->setActuatorType(Joint::VELOCITY);
    wheel2->setActuatorType(Joint::VELOCITY);
    wheel3->setActuatorType(Joint::VELOCITY);
    wheel4->setActuatorType(Joint::VELOCITY);
}

// Solve for a balanced pose using IK (Lesson 7 Answer)
Eigen::VectorXd solveIK(SkeletonPtr biped)
{
    // Modify the intial pose to one-foot stance before IK
    /*
    biped->setPosition(biped->getDof("j_shin_right")->
                       getIndexInSkeleton(), -1.4);
    biped->setPosition(biped->getDof("j_bicep_left_x")->
                       getIndexInSkeleton(), 0.8);
    biped->setPosition(biped->getDof("j_bicep_right_x")->
                       getIndexInSkeleton(), -0.8);
    */
    Eigen::VectorXd newPose = biped->getPositions();
    BodyNodePtr leftHeel = biped->getBodyNode("h_heel_left");
    BodyNodePtr leftToe = biped->getBodyNode("h_toe_left");
    double initialHeight = -0.8;
    
    for(std::size_t i = 0; i < default_ik_iterations; ++i)
    {
        Eigen::Vector3d deviation = biped->getCOM() - leftHeel->getCOM();
        Eigen::Vector3d localCOM = leftHeel->getCOM(leftHeel);
        LinearJacobian jacobian = biped->getCOMLinearJacobian() -
        biped->getLinearJacobian(leftHeel, localCOM);
        
        // Sagittal deviation
        double error = deviation[0];
        Eigen::VectorXd gradient = jacobian.row(0);
        Eigen::VectorXd newDirection = -0.2 * error * gradient;
        
        // Lateral deviation
        error = deviation[2];
        gradient = jacobian.row(2);
        newDirection += -0.2 * error * gradient;
        
        // Position constraint on four (approximated) corners of the left foot
        Eigen::Vector3d offset(0.0, -0.04, -0.03);
        error = (leftHeel->getTransform() * offset)[1] - initialHeight;
        gradient = biped->getLinearJacobian(leftHeel, offset).row(1);
        newDirection += -0.2 * error * gradient;
        
        offset[2] = 0.03;
        error = (leftHeel->getTransform() * offset)[1] - initialHeight;
        gradient = biped->getLinearJacobian(leftHeel, offset).row(1);
        newDirection += -0.2 * error * gradient;
        
        offset[0] = 0.04;
        error = (leftToe->getTransform() * offset)[1] - initialHeight;
        gradient = biped->getLinearJacobian(leftToe, offset).row(1);
        newDirection += -0.2 * error * gradient;
        
        offset[2] = -0.03;
        error = (leftToe->getTransform() * offset)[1] - initialHeight;
        gradient = biped->getLinearJacobian(leftToe, offset).row(1);
        newDirection += -0.2 * error * gradient;
        
        newPose += newDirection;
        biped->setPositions(newPose);
        biped->computeForwardKinematics(true, false, false);
    }
    return newPose;
}

SkeletonPtr createFloor()
{
    SkeletonPtr floor = Skeleton::create("floor");
    
    // Give the floor a body
    BodyNodePtr body =
    floor->createJointAndBodyNodePair<WeldJoint>(nullptr).second;
    
    // Give the body a shape
    double floor_width = 10.0;
    double floor_height = 0.01;
    std::shared_ptr<BoxShape> box(
                                  new BoxShape(Eigen::Vector3d(floor_width, floor_height, floor_width)));
    auto shapeNode
    = body->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(box);
    shapeNode->getVisualAspect()->setColor(dart::Color::Black());
    
    // Put the body into position
    Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
    tf.translation() = Eigen::Vector3d(0.0, -1, 0.0);
    body->getParentJoint()->setTransformFromParentBodyNode(tf);
    
    return floor;
}

int main(int argc, char* argv[])
{
    SkeletonPtr floor = createFloor();
    
    // Lesson 1
    SkeletonPtr biped = loadBiped();
    
    // Lesson 2
    //setInitialPose(biped);
    
    // Lesson 5
    //modifyBipedWithSkateboard(biped);
    
    // Lesson 6
    //setVelocityAccuators(biped);
    
    // Lesson 7
    //Eigen::VectorXd balancedPose = solveIK(biped);
    //biped->setPositions(balancedPose);
    
    WorldPtr world = std::make_shared<World>();
    world->setGravity(Eigen::Vector3d(0.0, 0.0, 0.0)); //-9.81
    
    if (dart::collision::CollisionDetector::getFactory()->canCreate("bullet"))
    {
        world->getConstraintSolver()->setCollisionDetector(
                                                           dart::collision::CollisionDetector::getFactory()->create("bullet"));
    }
    
    world->addSkeleton(floor);
    biped->setName("biped");
    
    //std::cout<<biped->frame;
    world->addSkeleton(biped);
    
    // Create a window for rendering the world and handling user input
    MyWindow window(world);
    
    // Print instructions
    std::cout << "'.': forward push" << std::endl;
    std::cout << "',': backward push" << std::endl;
    std::cout << "'s': increase skateboard forward speed" << std::endl;
    std::cout << "'a': increase skateboard backward speed" << std::endl;
    std::cout << "space bar: simulation on/off" << std::endl;
    std::cout << "'p': replay simulation" << std::endl;
    std::cout << "'v': Turn contact force visualization on/off" << std::endl;
    std::cout <<
    "'[' and ']': replay one frame backward and forward" << std::endl;
    
    // Initialize glut, initialize the window, and begin the glut event loop
    glutInit(&argc, argv);
    window.initWindow(640, 480, "character from .bvh file");
    glutMainLoop();
}



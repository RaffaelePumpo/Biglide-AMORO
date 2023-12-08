from lab_amoro.parallel_robot import *
from lab_amoro.plot_tools import *
from biglide_models import *  # Modify this to use the biglide
import sys
import numpy as np
import matplotlib.pyplot as plt
# write down a function to compute the trajectory in Cartesian space as a numpy array given some initial position, final position and duration. The function should output position, velocity and acceleration for each coordinate. The trajectory should be sampled with a 10 ms
#sampling time. Write down a second function to compute the trajectory (position, velocity,acceleration) in joint space from Cartesian space using the inverse geometric
# and kinematic models.

def compute_trajectory(initial_position, final_position, duration):
    t = np.arange(0,duration,0.01) 
    s= 10*(t/duration)**3 - 15*(t/duration)**4 + 6*(t/duration)**5
    position = initial_position+s*(final_position-initial_position)
    sd = (30/duration)*(t/duration)**2-(60/duration)*(t/duration)**3+(30/duration)*(t/duration)**4
    sdd = (60/duration**2)*(t/duration)-(180/duration**2)*(t/duration)**2+(120/duration**2)*(t/duration)**3
    velocity = sd*(final_position-initial_position)
    acceleration  = sdd*(final_position-initial_position)
    return position, velocity, acceleration

def compute_trajectory_joint(position, velocity, acceleration):
    num_samples = len(position[0])
    q11, q21 = np.zeros(num_samples), np.zeros(num_samples)
    q12, q22 = np.zeros(num_samples), np.zeros(num_samples)
    q11D, q21D = np.zeros(num_samples), np.zeros(num_samples)
    q12D, q22D = np.zeros(num_samples), np.zeros(num_samples)
    q11DD, q21DD = np.zeros(num_samples), np.zeros(num_samples)
    q12DD, q22DD = np.zeros(num_samples), np.zeros(num_samples)
    
    for t in range(num_samples):
        q11[t], q21[t] = igm(position[0, t], position[1, t], -1, -1)
        q12[t], q22[t] = dgm_passive(q11[t], q21[t], -1)
        q11D[t], q21D[t] = ikm(q11[t], q12[t], q21[t], q22[t], velocity[0, t], velocity[1, t])
        q12D[t], q22D[t] = dkm_passive(q11[t], q12[t], q21[t], q22[t], q11D[t], q21D[t], velocity[0,t], velocity[1,t])
        q11DD[t], q21DD[t] = ikm2(q11[t], q12[t], q21[t], q22[t], q11D[t], q12D[t], q21D[t], q22D[t], acceleration[0, t], acceleration[1, t])
        q12DD[t], q22DD[t] = dkm2_passive(q11[t], q12[t], q21[t], q22[t], q11D[t], q12D[t], q21D[t], q22D[t], q11DD[t], q21DD[t], acceleration[0,t], acceleration[1,t])
    
    active_joint_p = np.array([q11, q21])
    passive_joint_p = np.array([q12, q22])
    active_joint_v = np.array([q11D, q21D])
    passive_joint_v = np.array([q12D, q22D])
    active_joint_a = np.array([q11DD, q21DD])
    passive_joint_a = np.array([q12DD, q22DD])
    return active_joint_p, passive_joint_p, active_joint_v, passive_joint_v, active_joint_a, passive_joint_a


def main(args=None):
    # Initialize and start the ROS2 robot interface
    rclpy.init(args=args)
    robot = Robot("biglide")  # Modify this to use the biglide
    start_robot(robot)

    # Prepare plots
    app = QtGui.QApplication([])
    scope_joint1 = Scope("Joint 1", -0.5, 1.5)
    scope_joint2 = Scope("Joint 2", -1.5, 1.5)
    q11 = robot.active_left_joint.position
    q21 = robot.active_right_joint.position
    x,y = dgm(q11, q21, -1)
    initial_position = np.array([[x], [y]])
    #final_position = np.array([[x], [y+1]])   # For the trajectory along y
    final_position = np.array([[x+0.08], [y]])   # For the  trajectory moving joints in different way
    duration = 2.0
    # Create the trajectory as arrays in Cartesian space (initial position, final position, duration)
    positions, velocities, accelerations = compute_trajectory(initial_position, final_position, duration)
    # Create the trajectory as arrays in joint space using the inverse models (position, velocity, acceleration)

    active_joint_p, passive_joint_p, active_joint_v, passive_joint_v, active_joint_a, passive_joint_a = compute_trajectory_joint(positions, velocities, accelerations)
    index = 0
	

    #error_x = []
    #error_y = []
    # Controller gains
    # Define Kp and Kd for the active joints as diagonal matrices
    Kp = np.diag([30.0, 30.0])
    Kd = np.diag([15.0, 15.0])
    # Controller
    try:
        robot.apply_efforts(0.0, 0.0)  # Required to start the

        while True:
            if robot.data_updated():
                # Robot available data - This is the only data thet you can get from a real robot (joint encoders)
                q11 = robot.active_left_joint.position
                q21 = robot.active_right_joint.position
                q11D = robot.active_left_joint.velocity
                q21D = robot.active_right_joint.velocity

                # CTC controller
                q12,q22 = dgm_passive(q11, q21, 1)
                q12D,q22D = dkm_passive(q11, q12, q21, q22, q11D, q21D, 1, 1)
                M,c = dynamic_model(q11, q12, q21, q22, q11D, q12D, q21D, q22D)
               
                qa = np.array([q11,q21])
                qaD = np.array([q11D,q21D])
                alpha = active_joint_a[:,index]+Kd.dot(active_joint_v[:,index]-qaD)+Kp.dot(active_joint_p[:,index]-qa)
                tau = np.matmul(M,alpha)+c
                tau_left = tau[0]
                tau_right = tau[1]
                robot.apply_efforts(tau_left,tau_right)
                # Scope update
                time = robot.get_time()
                if time < 5.0:
                    scope_joint1.update(time,active_joint_p[0,index],q11)
                    scope_joint2.update(time, active_joint_p[1,index],q21)
                    #error_x.append(active_joint_p[0,index]-q11)
                    #error_y.append(active_joint_p[1,index]-q21)
                if index < len(positions[0])-1:
                    index += 1  # Next point in trajectory

    except KeyboardInterrupt:
        
        """plt.figure()
        plt.plot(error_x)
        plt.title("Error q11")
        plt.xlabel("Time")
        plt.ylabel("Error")
        plt.grid()
        plt.figure()
        plt.plot(error_y)
        plt.title("Error q21")
        print("mean = ",np.mean(error_y))
        print("max = ",np.max(error_y))
        plt.xlabel("Time [ms]")
        plt.ylabel("Error[m]")
        plt.grid()

        plt.show()"""
        pass


if __name__ == "__main__":
    main(sys.argv)


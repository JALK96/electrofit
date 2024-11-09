How to proceed: 
--------------


First run 1_pis.py
    Wait untill processing comletion (when the .acpype folder has been created in the directories: process/IP_XXXXXX/run_gau_create_gmx_in)

Secondly run 2_setup_sim_dir.py
    Then manually go to each directory (process/IP_XXXXXX/run_gmx_simulation) and specify the qcmXX machine u want to run the simulation on in the "gmx.sh" bash script. 

Thirdly run 3_start_sim.py
    This will run each gmx.sh bash script in the run_gmx_simulation directories and therby start a simulation for each config of IP6 on the specified machines (qcmXX).
    
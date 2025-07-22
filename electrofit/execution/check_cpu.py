#!/usr/bin/env python3
import paramiko
import threading

# List of remote machines
MACHINES = ["qcm04", "qcm05", "qcm06", "qcm07", "qcm08"]

# Command to check CPU utilization per core using mpstat.
CPU_USAGE_CMD = "mpstat -P ALL 1 1"

# Command to get process CPU usage by user.
PROCESS_USAGE_CMD = "ps -eo user,pcpu,comm --sort=-pcpu"

# SSH credentials
USERNAME = "johannal96"  # Replace with your SSH username
KEY_FILE = "/home/johannal96/.ssh/id_rsa"  # Replace with your SSH key path or set to None if using password.
PASSWORD = None  # Set a password if not using key-based authentication.

# Optional: set SSH connection timeout in seconds
SSH_TIMEOUT = 10

# Global dictionaries to store summary results.
# summary_results maps host -> (active_cores, total_cores)
summary_results = {}
# process_summary maps host -> { user: aggregated_cpu_usage }
process_summary = {}

# Locks for thread-safe updates
summary_lock = threading.Lock()
process_lock = threading.Lock()

def check_cpu_usage(host):
    """
    Connects to the given host via SSH, executes commands to get CPU core utilization
    and process-level CPU usage per user.
    """
    active_count = 0
    total_cores = 0
    user_cpu = {}  # Dictionary to aggregate CPU usage by user

    try:
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        # Connect using either key-based or password authentication
        if KEY_FILE:
            client.connect(hostname=host, username=USERNAME, key_filename=KEY_FILE, timeout=SSH_TIMEOUT)
        else:
            client.connect(hostname=host, username=USERNAME, password=PASSWORD, timeout=SSH_TIMEOUT)

        # Execute the mpstat command for per-core utilization
        stdin, stdout, stderr = client.exec_command(CPU_USAGE_CMD)
        mpstat_output = stdout.read().decode('utf-8')
        mpstat_error = stderr.read().decode('utf-8')

        if not mpstat_error:
            for line in mpstat_output.splitlines():
                tokens = line.split()
                if len(tokens) < 2:
                    continue
                if tokens[1].lower() in ["cpu", "all"]:
                    continue
                try:
                    idle = float(tokens[-1])
                    utilization = 100.0 - idle
                    total_cores += 1
                    if utilization > 1.0:
                        active_count += 1
                except ValueError:
                    continue

        # Execute the ps command to get process CPU usage per user
        stdin, stdout, stderr = client.exec_command(PROCESS_USAGE_CMD)
        ps_output = stdout.read().decode('utf-8')
        ps_error = stderr.read().decode('utf-8')

        if not ps_error:
            lines = ps_output.splitlines()
            if lines:
                for line in lines[1:]:  # Skip header line.
                    tokens = line.split(None, 2)  # Split into 3 parts: user, pcpu, command
                    if len(tokens) < 2:
                        continue
                    user = tokens[0]
                    try:
                        cpu_usage = float(tokens[1])
                    except ValueError:
                        continue
                    user_cpu[user] = user_cpu.get(user, 0.0) + cpu_usage

        # Update global summaries in a thread-safe manner
        with summary_lock:
            summary_results[host] = (active_count, total_cores)
        with process_lock:
            process_summary[host] = user_cpu

    except Exception as e:
        # In production, consider logging exceptions.
        with summary_lock:
            summary_results[host] = ("Error", "Error")
        with process_lock:
            process_summary[host] = {}
    finally:
        client.close()

def main():
    """
    Main function spawns a thread for each host, gathers data, and prints a summary
    including active cores and top 5 users by CPU usage.
    """
    threads = []
    for host in MACHINES:
        thread = threading.Thread(target=check_cpu_usage, args=(host,))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

    print("\nSummary of aggregated CPU usage per user (top 5) on each machine:")
    print("-----------------------------------------------------------------")
    for host in MACHINES:
        print(f"\n{host}:")
        user_cpu = process_summary.get(host, {})
        if user_cpu:
            top_users = sorted(user_cpu.items(), key=lambda item: item[1], reverse=True)[:5]
            for user, usage in top_users:
                print(f"  {user}: {usage:.2f}% CPU usage")
        else:
            print("  No process data available.")

    # Print final summaries
    print("\nSummary of active CPU cores per machine:")
    print("----------------------------------------\n")
    for host in MACHINES:
        active, total = summary_results.get(host, ("No data", "No data"))
        print(f"{host}: {active} active cores out of {total} total cores.")

if __name__ == "__main__":
    main()
import os
import sys
import subprocess
import shutil
import json
from typing import Dict, Any, List
import time
import glob
from ..util import strip_comment

class Simulation:
    """Class for handling NERDSS simulation configurations and running simulations.

    Attributes:
        work_dir (str): The working directory for the simulation.
    """
    
    def __init__(self, work_dir: str) -> None:
        """Initializes the Simulation class.
        
        Args:
            work_dir (str): The working directory for the simulation.
        
        Examples:
            sim = Simulation("~/my_simulation")
        """
        
        # Expand user home directory if path starts with ~
        if work_dir.startswith("~"):
            work_dir = os.path.expanduser(work_dir)
        
        # Convert to absolute path and store
        self.work_dir = os.path.abspath(work_dir)
        
        # Create working directory if it doesn't exist
        os.makedirs(self.work_dir, exist_ok=True)
        print(f"Working directory set to: {self.work_dir}")

        self.parmfile = 'parms.inp'
        self.coordinatefile = 'fixCoordinates.pdb'


    def modify_mol_file(self, mol_name: str, modifications: Dict[str, Any]) -> None:
        """Modifies the parameters of an existing .mol file.
        
        Args:
            mol_name (str): The name of the molecule to modify.
            modifications (Dict[str, Any]): A dictionary containing parameter modifications.
        
        Raises:
            FileNotFoundError: If the specified molecule file does not exist.
        """
        input_dir = os.path.join(self.work_dir, "nerdss_input")
        mol_file = os.path.join(input_dir, f"{mol_name}.mol")
        
        if not os.path.exists(mol_file):
            available_mols = [f.split(".mol")[0] for f in os.listdir(input_dir) if f.endswith(".mol")]
            raise FileNotFoundError(f"Molecule '{mol_name}' not found. Available molecules: {', '.join(available_mols)}")
        
        with open(mol_file, "r") as f:
            lines = f.readlines()
        
        with open(mol_file, "w") as f:
            for line in lines:
                key = line.split("=")[0].strip()
                if key in modifications:
                    f.write(f"{key} = {modifications[key]}\n")
                else:
                    f.write(line)

    def modify_inp_file(self, modifications: Dict[str, Any], filename: str = "parms.inp") -> None:
        """
        Modifies the parameters of the parms.inp file. If `isSphere` and `sphereR` are provided, 
        removes the `WaterBox` line and adds the new lines accordingly. If `WaterBox` is provided, 
        removes `isSphere` and `sphereR` if they exist.
        
        Args:
            modifications (Dict[str, Any]): A dictionary containing parameter modifications.
            filename (str): The name of the input file to modify. Defaults to "parms.inp".
        """
        if filename.strip() not in ["parms.inp", ""]:
            self.parmfile = filename
        inp_file = os.path.join(self.work_dir, "nerdss_input", self.parmfile)
        
        if not os.path.exists(inp_file):
            raise FileNotFoundError(f"{filename} file not found.")
        
        with open(inp_file, "r") as f:
            lines = f.readlines()
        
        modified_lines = []
        in_boundaries_section = False
        waterbox_removed = False
        sphere_removed = False
        in_molecules_section = False
        in_reactions_section = False
        current_reaction = None
        molecule_types = []

        for line in lines:
            stripped_line = line.strip()
            
            if stripped_line.startswith("start boundaries"):
                in_boundaries_section = True
                modified_lines.append(line)
                continue
            
            if stripped_line.startswith("end boundaries"):
                in_boundaries_section = False
                
                if "isSphere" in modifications and "sphereR" in modifications:
                    modified_lines.append(f"\tisSphere = {modifications['isSphere']}\n")
                    modified_lines.append(f"\tsphereR = {modifications['sphereR']}\n")
                elif "WaterBox" in modifications:
                    modified_lines.append(f"\tWaterBox = {modifications['WaterBox']}\n")
                
                modified_lines.append(line)
                continue
            
            if in_boundaries_section:
                if stripped_line.startswith("WaterBox") and "isSphere" in modifications and "sphereR" in modifications:
                    waterbox_removed = True
                    continue
                elif (stripped_line.startswith("isSphere") or stripped_line.startswith("sphereR")) and "WaterBox" in modifications:
                    sphere_removed = True
                    continue

            if stripped_line.startswith("start molecules"):
                in_molecules_section = True
                modified_lines.append(line)
                continue
            if stripped_line.startswith("end molecules"):
                in_molecules_section = False
                modified_lines.append(line)
                continue
            
            if stripped_line.startswith("start reactions"):
                in_reactions_section = True
                modified_lines.append(line)
                continue
            if stripped_line.startswith("end reactions"):
                in_reactions_section = False
                modified_lines.append(line)
                continue

            if in_molecules_section and ":" in stripped_line:
                mol_name, count = map(str.strip, stripped_line.split(":"))
                molecule_types.append(mol_name)
                if mol_name in modifications:
                    modified_lines.append(f"\t{mol_name} : {modifications[mol_name]}\n")
                else:
                    modified_lines.append(line)
                continue

            if in_reactions_section:
                if '=' not in stripped_line and stripped_line:
                    current_reaction = stripped_line.strip()
                    modified_lines.append(line)
                    continue
                if current_reaction and current_reaction in modifications:
                    param_name = stripped_line.split("=")[0].strip()
                    if param_name in modifications[current_reaction]:
                        modified_lines.append(f"\t\t{param_name} = {modifications[current_reaction][param_name]}\n")
                        continue
            
            key = stripped_line.split("=")[0].strip()
            if key in modifications:
                modified_lines.append(f"\t{key} = {modifications[key]}\n")
            else:
                modified_lines.append(line)
        
        with open(inp_file, "w") as f:
            f.writelines(modified_lines)

    def add_interface_state(self, mol_name: str, interface_name: str, states: List[str]) -> None:
        """Adds states to a specified interface of a molecule.
        
        Args:
            mol_name (str): The name of the molecule.
            interface_name (str): The name of the interface.
            states (List[str]): List of single-character state names.
        
        Raises:
            FileNotFoundError: If the molecule file does not exist.
            ValueError: If no valid states are provided.
        """
        if not states or any(len(state) != 1 for state in states):
            raise ValueError("States must be single-character values.")

        input_dir = os.path.join(self.work_dir, "nerdss_input")
        mol_file = os.path.join(input_dir, f"{mol_name}.mol")
        
        if not os.path.exists(mol_file):
            available_mols = [f.split(".mol")[0] for f in os.listdir(input_dir) if f.endswith(".mol")]
            raise FileNotFoundError(f"Molecule '{mol_name}' not found. Available molecules: {', '.join(available_mols)}")
        
        with open(mol_file, "a") as f:
            state_line = f"state = {interface_name}~" + "~".join(states) + "\n"
            f.write(state_line)

    def print_mol_parameters(self, mol_name: str) -> None:
        """Prints all parameters of a given .mol file.
        
        Args:
            mol_name (str): The name of the molecule to display.
        
        Raises:
            FileNotFoundError: If the specified molecule file does not exist.
        """
        input_dir = os.path.join(self.work_dir, "nerdss_input")
        mol_file = os.path.join(input_dir, f"{mol_name}.mol")
        
        if not os.path.exists(mol_file):
            available_mols = [f.split(".mol")[0] for f in os.listdir(input_dir) if f.endswith(".mol")]
            raise FileNotFoundError(f"Molecule '{mol_name}' not found. Available molecules: {', '.join(available_mols)}")
        
        with open(mol_file, "r") as f:
            print(f"Parameters for molecule '{mol_name}':")
            print(f.read())

    def print_inp_file(self, file_name: str = "parms.inp") -> None:
        """
        Prints the contents of the parms.inp file.
        
        Args:
            file_name (str): The name of the input file to print. Defaults to "parms.inp".
        """

        if file_name.strip() not in ["parms.inp", ""]:
            self.parmfile = file_name

        inp_file = os.path.join(self.work_dir, "nerdss_input", self.parmfile)
        
        if not os.path.exists(inp_file):
            print(f"{self.parmfile} file not found.")
            return
        
        with open(inp_file, "r") as f:
            print(f.read())

    def install_nerdss(self, nerdss_path: str = None) -> None:
        """Installs the NERDSS package.

        Args:
            nerdss_path (str): The path to install NERDSS. If None, uses the current directory.
        """
        if nerdss_path is None:
            nerdss_path = os.getcwd()

        if nerdss_path.startswith("~"):
            nerdss_path = os.path.expanduser(nerdss_path)
        nerdss_path = os.path.abspath(nerdss_path)
        
        nerdss_repo_path = os.path.join(nerdss_path, "NERDSS")
        
        # Ensure target directory exists
        os.makedirs(nerdss_path, exist_ok=True)
        print(f"Installing NERDSS to {nerdss_path}...")

        # Check if git and make are installed
        for cmd in ["git", "make"]:
            if shutil.which(cmd) is None:
                print(f"Error: {cmd} is not installed. Please install it and try again.")
                return

        # Clone the repository if it doesn't exist
        if not os.path.exists(nerdss_repo_path):
            result = subprocess.run(["git", "clone", "https://github.com/mjohn218/NERDSS.git", nerdss_repo_path], check=False)
            if result.returncode != 0:
                print("Error: Failed to clone the NERDSS repository.")
                return
        else:
            print("NERDSS repository already exists. Pulling latest updates...")
            subprocess.run(["git", "-C", nerdss_repo_path, "pull"], check=False)

        def detect_package_manager():
            """Detects the package manager for the current Linux distribution."""
            if os.path.exists("/etc/os-release"):
                with open("/etc/os-release", "r") as f:
                    os_release = f.read().lower()
                    if "ubuntu" in os_release or "debian" in os_release:
                        return "apt"
                    elif "fedora" in os_release:
                        return "dnf"
                    elif "centos" in os_release or "rhel" in os_release:
                        return "yum"
                    elif "opensuse" in os_release:
                        return "zypper"
            return None

        # Install dependencies based on the platform
        if sys.platform.startswith("linux"):
            package_manager = detect_package_manager()
            install_command = None

            if package_manager == "apt":
                install_command = ["sudo", "apt-get", "install", "-y", "build-essential", "libgsl-dev"]
            elif package_manager == "dnf":
                install_command = ["sudo", "dnf", "install", "-y", "gcc", "gcc-c++", "make", "gsl-devel"]
            elif package_manager == "yum":
                install_command = ["sudo", "yum", "install", "-y", "gcc", "gcc-c++", "make", "gsl-devel"]
            elif package_manager == "zypper":
                install_command = ["sudo", "zypper", "install", "-y", "gcc", "gcc-c++", "make", "gsl-devel"]
            else:
                print("Skipping system package installation. Ensure GSL is installed manually.")
                print("sudo apt-get install build-essential libgsl-dev  # For Debian/Ubuntu")
                print("sudo dnf install gcc gcc-c++ make gsl-devel     # For Fedora")
                print("sudo yum install gcc gcc-c++ make gsl-devel     # For CentOS/RHEL")
                print("sudo zypper install gcc gcc-c++ make gsl-devel  # For openSUSE")
            
            if install_command:
                result = subprocess.run(["sudo", package_manager, "update", "-y"], check=False)
                if result.returncode == 0:
                    subprocess.run(install_command, check=False)
                else:
                    print("Skipping system package installation. Ensure GSL is installed manually.")
                    print("sudo apt-get install build-essential libgsl-dev  # For Debian/Ubuntu")
                    print("sudo dnf install gcc gcc-c++ make gsl-devel     # For Fedora")
                    print("sudo yum install gcc gcc-c++ make gsl-devel     # For CentOS/RHEL")
                    print("sudo zypper install gcc gcc-c++ make gsl-devel  # For openSUSE")
        elif sys.platform == "darwin":
            subprocess.run(["brew", "install", "gsl"], check=False)

        # Compile NERDSS
        make_result = subprocess.run(["make", "serial"], cwd=nerdss_repo_path, check=False)
        if make_result.returncode == 0:
            print("NERDSS installation complete.")
        else:
            print("Error: Compilation failed. Please check the logs and dependencies.")

    def run_new_simulations(
            self, sim_indices: List[int] = None, sim_dir: str = None, nerdss_dir: str = None, parallel: bool = False,
            coordinate: bool = False, progress: bool = True, verbose=True,
        ) -> None:
        """Runs NERDSS simulations based on the given parameters.
        
        Args:
            sim_indices (List[int], optional): List of simulation indices to run. If None, runs one simulation with index = 1.
            sim_dir (str, optional): Directory where simulation results should be stored. Defaults to `self.work_dir/nerdss_output`.
            nerdss_dir (str, optional): Directory where NERDSS is installed. Defaults to `self.work_dir/NERDSS`.
            parallel (bool, optional): Whether to run simulations in parallel. Defaults to False.

        Notes:
            FIXME: Doesn't work on Fedora OS using Jupyter notebook. Doesn't test on other OS. Doesn't test using Python script.
        """
        if sim_dir is None:
            sim_dir = os.path.join(self.work_dir, "nerdss_output")
        elif sim_dir.startswith("~"):
            sim_dir = os.path.expanduser(sim_dir)
        sim_dir = os.path.abspath(sim_dir)
        os.makedirs(sim_dir, exist_ok=True)

        if nerdss_dir is None:
            nerdss_dir = os.path.join(self.work_dir, "NERDSS")
        elif nerdss_dir.startswith("~"):
            nerdss_dir = os.path.expanduser(nerdss_dir)
        nerdss_dir = os.path.abspath(nerdss_dir)
        # check whether nerdss executable exists
        nerdss_exec = os.path.join(nerdss_dir, "bin", "nerdss")
        if not os.path.exists(nerdss_exec):
            raise FileNotFoundError(f"NERDSS executable not found at {nerdss_exec}. Make sure it is installed and compiled.")
        
        input_dir = os.path.join(self.work_dir, "nerdss_input")
        parms_file = os.path.join(input_dir, self.parmfile)
        
        if not os.path.exists(parms_file):
            raise FileNotFoundError(f"NERDSS input file not found: {parms_file}")
        
        if sim_indices is None:
            sim_indices = [1]
        
        processes = []
        progress_bars = {}
        
        for index in sim_indices:
            sim_subdir = os.path.join(sim_dir, f"{index}")
            os.makedirs(sim_subdir, exist_ok=True)
            
            for file in os.listdir(input_dir):
                shutil.copy(os.path.join(input_dir, file), sim_subdir)
            shutil.copy(nerdss_exec, sim_subdir)
            
            output_log = os.path.join(sim_subdir, "output.log")
            with open(output_log, "w") as log_file:
                cmd = ["./nerdss", "-f", self.parmfile]
                if coordinate:
                    cmd.append("-c")
                    cmd.append(self.coordinatefile)
                
                if parallel:
                    process = subprocess.Popen(cmd, cwd=sim_subdir, stdout=log_file, stderr=log_file)
                    processes.append((index, process))
                else:
                    if verbose: print(f"Running simulation {index}...")
                    process = subprocess.Popen(cmd, cwd=sim_subdir, stdout=log_file, stderr=log_file)

                    if verbose:
                        if progress:
                            from tqdm import tqdm
                            progress_bars[index] = tqdm(total=100, desc=f"Simulation {index}")
                            while process.poll() is None:
                                progress = self.calculate_progress_percentage(sim_subdir)
                                progress_bars[index].n = progress
                                progress_bars[index].refresh()
                                time.sleep(2)
                            
                            progress_bars[index].close()
        
        if parallel and verbose:
            for index, process in processes:
                print(f"Waiting for simulation {index} to complete...")
                process.wait()
                print(f"Simulation {index} completed.")
        
        if verbose:
            if progress: 
                print("All simulations completed.")

    def calculate_progress_percentage(self, sim_subdir: str) -> int:
        """
        Calculates the progress percentage of a running simulation.
        
        Args:
            sim_subdir (str): The directory of the simulation.
        """
        current_time = 0.0
        copy_numbers_file = os.path.join(sim_subdir, "copy_numbers_time.dat")
        
        try:
            cmd = f"tail -n 2 {copy_numbers_file}"
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = process.communicate()
            output = output.decode("utf-8").strip()
            lines = output.split("\n")
            last_line = lines[-1] if len(lines) > 1 else lines[0]
            current_time = float(last_line.split(",")[0])
        except Exception:
            pass
        
        total_time = 0.0
        try:
            nItr = 0
            timeStep = 0.0
            parm_file = os.path.join(sim_subdir, self.parmfile)
            with open(parm_file, "r") as inp_file:
                for line in inp_file:
                    line = strip_comment(line.strip())
                    if line.startswith("nItr"):
                        nItr = int(line.split("=")[1])
                    if line.startswith("timeStep"):
                        timeStep = float(line.split("=")[1]) * 1e-6
            total_time = nItr * timeStep
        except Exception as e:
            raise e
        
        if total_time == 0.0:
            return 0
        else:
            return int(current_time / total_time * 100)

    def run_restart_simulations(self, sim_indices: List[int] = None, sim_dir: str = None, nerdss_dir: str = None, restart_from: str = "", restart_sim_name: str = "restart_sim", parallel: bool = False) -> None:
        """Runs NERDSS simulations from a restart file.
        
        Args:
            sim_indices (List[int], optional): List of simulation indices to restart. If None, restarts one simulation with index = 1.
            sim_dir (str, optional): Directory where restarted simulation results should be stored. Defaults to `self.work_dir/nerdss_output`.
            nerdss_dir (str, optional): Directory where NERDSS is installed. Defaults to `self.work_dir/NERDSS`.
            restart_from (str): Path to the directory containing the restart file.
            restart_sim_name (str): Name of the folder where restarted simulations will be stored.
            parallel (bool, optional): Whether to run simulations in parallel. Defaults to False.
        """
        if sim_dir is None:
            sim_dir = os.path.join(self.work_dir, "nerdss_output")
        os.makedirs(sim_dir, exist_ok=True)

        if nerdss_dir is None:
            nerdss_dir = os.path.join(self.work_dir, "NERDSS")
        
        nerdss_exec = os.path.join(nerdss_dir, "bin", "nerdss")
        if not os.path.exists(nerdss_exec):
            raise FileNotFoundError(f"NERDSS executable not found at {nerdss_exec}. Make sure it is installed and compiled.")
        
        if sim_indices is None:
            sim_indices = [1]
        
        for index in sim_indices:
            restart_subdir = os.path.join(sim_dir, f"{index}", f"{restart_sim_name}")
            os.makedirs(restart_subdir, exist_ok=True)
            
            if restart_from == "":
                restart_file = os.path.join(sim_dir, f"{index}", "DATA", "restart.dat")
            else:
                restart_file = os.path.join(sim_dir, f"{index}", restart_from, "DATA", "restart.dat")
            if not os.path.exists(restart_file):
                raise FileNotFoundError(f"Restart file not found at {restart_file}.")
            
            shutil.copy(restart_file, restart_subdir)
            shutil.copy(nerdss_exec, restart_subdir)
            
            output_log = os.path.join(restart_subdir, "output.log")
            with open(output_log, "w") as log_file:
                cmd = ["./nerdss", "-r", "restart.dat"]
                
                if parallel:
                    subprocess.Popen(cmd, cwd=restart_subdir, stdout=log_file, stderr=log_file)
                else:
                    print(f"Restarting simulation {index}...")
                    process = subprocess.Popen(cmd, cwd=restart_subdir, stdout=log_file, stderr=log_file)
                    process.wait()
                    print(f"Simulation {index} restarted successfully.")
        
        print("All restart simulations completed.")


    def _print_dict(self,dict):
        '''
        prints the output of pull reaction information, pull parameter information, and pull mol information for copying by the user. 


        This function iterates over each key-value pair in the input dictionary and prints it in the format:
        ['key'] = value. The user can than copy and paste their the key-value pairs they want to edit to their defined dictionary.
        Subsquently, that defined dictionary can be passed back through modify inp, or modify mol prior to simulation start. Useful for looping through conditions. 


        Parameters:
        -----------
        dict : dict
            The dictionary whose contents are to be printed.

        Returns:
        --------
        None
            This function prints to standard output and does not return any value.
        '''
        for key, value in dict.items():
            key_str = f"'{key}'" if isinstance(key, str) else str(key)

            if isinstance(value, list):
                value_str = f"[{', '.join(map(str, value))}]"
            else:
                value_str = str(value)
            print(f'[{key_str}] = {value_str}')
        print("\n")

    def pull_reaction_information(self,file: str):
        '''
        Extracts reaction information from a given input file and returns it as a dictionary.

        The function parses a file (e.g., "parms.inp") to extract details about reactions within a specific block 
        labeled by "start reactions" and "end reactions". For each reaction, the function captures relevant 
        information and organizes it into a dictionary, where each key corresponds to a reaction equation 
        (e.g., "A <-> B") and its associated parameters.

        Reaction details are stored as nested dictionaries, with the reaction equation as the outer key, 
        and each parameter as an inner key-value pair.

        Parameters:
        -----------
        file : str
            The path to the input file containing the reaction information.

        Returns:
        --------
        dict
            A dictionary where each key is a reaction equation (e.g., "A <-> B") and the value is another dictionary
            containing parameters and values associated with that reaction. The parameters may include things like 
            exclusion conditions and numerical values for reaction conditions.
        
        Example:
        --------
        Given an input file containing reaction data, the function will return a dictionary like:
        {
            "A <-> B": {
                "norm1": [1.0, 2.0, 3.0],
                "sigma": "1.20302012"
            },
            "C -> D": {
                "onRate": "0"
            }
        }
        
        Notes:
        ------
        - The function assumes the input file contains structured reaction information in blocks marked by "start reactions" and "end reactions".
        - Each reaction line may contain additional parameters, which are processed as key-value pairs.
        - Lines with "exclude" are treated specially, storing them in the dictionary under the respective reaction.
        '''

        rxn_dict = {}
        with open(file,"r") as f:
            lines = f.readlines()
            in_reactions: bool = False
            in_rxn_block: bool = False
            current_rxn: str = ""
            for line in lines:
                line = line.strip()
                if line.startswith("start reactions"):
                    in_reactions = True
                    continue
                if line.startswith("end reactions"):
                    in_reactions = False
                    continue
                if in_reactions == False:
                    continue
                if "<->" in line or "<-" in line or "->" in line:
                    in_rxn_block = True
                    current_rxn = line
                    rxn_dict[line] = {}
                    continue
                if "exclude" in line and in_rxn_block == True:
                    rxn_dict[current_rxn][line.split()[0]] = line.split()[-1]
                    in_rxn_block = False
                    continue
                if in_rxn_block == True and in_reactions == True:
                    #generally the .split() outputs [condition, =, vals] so anything from 2 onward is our vals
                    if len((line.split()[2:])) > 2:
                        rxn_dict[current_rxn][line.split()[0]] = [float(x) for x in line.replace("[","").replace(",", "").replace(']',"").split()[2:]]
                        continue
                    print(line.split())
                    if line.split() == []:
                        continue
                    else:
                        rxn_dict[current_rxn][line.split()[0]] = line.split()[2]
            print("The following lines can be used to access your reaction information. Copy and Paste the reactions into your code you wish to modify. Be sure to include the dictionary name.")
            self._print_dict(rxn_dict)
            return rxn_dict
            
    def pull_parameter_file_information(self,file: str):
        '''
        Parses a simulation input file and extracts parameter, boundary, and molecule information into a dictionary.

        Parameters:
        -----------
        file : str
            Path to the input file containing the simulation parameters and configuration data.

        Returns:
        --------
        dict
            A dictionary containing key-value pairs from the parameters, boundaries, and molecules blocks.
            - Keys are parameter names (e.g., "dt", "runtime").
            - Values are strings, floats, or lists of floats depending on the format of the line in the file.

        Example:
        --------
        Given an input file, the function might return:
        {
            "nItr": "10000",
            "timeStep": "0.1",
            'WaterBox' = [100, 100, 100],
        }
        '''
        
        with open(file,"r") as f:
            lines = f.readlines()
            param_dict: dict = {}
            in_params: bool = False
            in_bounds: bool = False
            in_mol: bool = False
        for line in lines:
            line = line.strip()
            if line.startswith("start parameters"):
                in_params = True
                continue
            if line.startswith("end parameters"):
                in_params = False
                continue
            if in_params == True:
                if "#iterations" in line:
                    param_dict[line.split()[0]] = line.split()[2]
                param_dict[line.split()[0]] = line.split()[2]
            line = line.strip()
            if line.startswith("start boundaries"):
                in_bounds = True
                continue
            if line.startswith("end boundaries"):
                in_bounds = False
                continue
            if line.startswith("start molecules"):
                in_mol = True
                continue
            if line.startswith("end molecules"):
                in_mol = False
            if in_bounds == True:
                if len((line.split()[2:])) > 2:
                        param_dict[line.split()[0]] = [float(x) for x in line.replace("[","").replace(",", "").replace(']',"").split()[2:]]
                        continue
                param_dict[line.split()[0]] = line.split()[2]
            if in_mol == True:
                param_dict[line.split()[0]] = line.split()[2]
        print("The following lines can be used to access your parameter information. Copy and Paste the parameters into your code you wish to modify. Be sure to include the dictionary name.")

        self._print_dict(param_dict)

        return param_dict

    def pull_mol_file_information(self,file: str):
        '''
        Extracts molecular configuration information from a .mol-style input file into a dictionary.

        Parameters:
        -----------
        file : str
            Path to the input .mol file containing molecular relationship and configuration information.

        Returns:
        --------
        dict
            A dictionary where keys are molecular attribute names (e.g., "mass", "COM", "D" for example) and 
            values are either strings or lists of floats depending on the format in the file.

        Example:
        --------
        Given a section of a .mol file, the output might look like:
        {
            "mass": "1.0",
            "COM": [0.0, 0.0, 1.0],
            "D" = [13.0, 13.0, 13.0]
        }
        '''
        with open(file,'r') as f:
            lines = f.readlines()
            mol: dict = {}
            in_rel = False
        for line in lines:
            line = line.strip()
            if line.startswith("Name"):
                in_rel = True
                continue
            if len(line.split()) == 0: #skip empty lines
                continue
            if line.startswith("COM"):
                in_rel = False
                continue
            if in_rel == True:
                if len((line.split()[2:])) > 2:
                        mol[line.split()[0]] = [float(x) for x in line.replace("[","").replace(",", "").replace(']',"").split()[2:]]
                        continue
                mol[line.split()[0]] = line.split()[2]
            
        print(("The following lines can be used to access your mol information. Copy and paste this output into your code to modify the .mol file"))
        self._print_dict(mol)
        return mol

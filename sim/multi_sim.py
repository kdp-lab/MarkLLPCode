import subprocess
import multiprocessing
import argparse
import os
from codes import functions

def process_input_file(input_file, output_directory, number_of_events, skipEvents):
    print("call extract_params")
    # Extract mass, width, and lifetime from the input file
    mass, width, lifetime = functions.extract_params(input_file)
    # Generate modified .tbl file
    print("call generate_modified_tbl")
    tbl_file = functions.generate_modified_tbl(mass, lifetime, width, input_file)
    # Run DDSIM simulation
    print("call run_ddsim")
    functions.run_ddsim(input_file, output_directory, tbl_file, number_of_events, skipEvents)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run multiple DDSIM simulations simultaneously.")
    parser.add_argument("input_files", nargs="+", help="List of input files for simulations.")
    parser.add_argument("-n", "--number_of_events", help="Number of events to simulate.", type=int, default=1000)
    parser.add_argument("-o", "--output_directory", help="Output directory for simulation results.", default="/home/larsonma/MarkLLPCode/sim")
    #parser.add_argument("-o", "--output_directory", help="Output directory for simulation results.", default="/scratch/larsonma/tutorial2024/LLPStudies/MarkLLPCode/sim")
    parser.add_argument("-j", "--ncpu", help="Number of CPU cores to use.", type=int, default=1)
    parser.add_argument("-s", "--skipEvents", help="Number of events to simulate.", type=int, default=0)
    args = parser.parse_args()
    
    # Prepend "/local/d1/mu+mu-/samples/' to each input file path
    input_files = [f"/scratch/larsonma/tutorial2024/LLPStudies/samples/{input_file}" for input_file in args.input_files]

    # Create a pool of processes for parallel execution
    pool = multiprocessing.Pool(args.ncpu)

    # Process each input file
    pool.starmap(process_input_file, [(input_file, args.output_directory, args.number_of_events, args.skipEvents) for input_file in input_files])

    # Close the pool of processes
    pool.close()
    pool.join()


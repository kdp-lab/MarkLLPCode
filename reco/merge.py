import os
import pyLCIO

# Directory where the split files are located
input_dir = "reco_osg_condor2"

# Output file to store the combined events
output_file = "4500_10_osg_reco.slcio"

# Create LCWriter to write combined events to the output file
writer = pyLCIO.IOIMPL.LCFactory.getInstance().createLCWriter()
writer.open(output_file)

# Loop over the 100 event files (starting from i=0)
for i in range(100):
    input_filename = os.path.join(input_dir, f"4500_10_reco{i}.slcio")
    
    if os.path.exists(input_filename):
        # Create LCReader to read the single event file
        reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
        reader.open(input_filename)
        
        # Read the event and write it to the combined output file
        for event in reader:
            writer.writeEvent(event)
        
        # Close the reader for the current event file
        reader.close()
        
        print(f"Event from {input_filename} added to {output_file}")
    else:
        print(f"File {input_filename} does not exist")

# Close the writer for the combined output file
writer.close()

print(f"Successfully combined events into {output_file}")

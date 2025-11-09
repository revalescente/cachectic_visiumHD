# the comands need to be run from the Home "/mnt/europa/valerio/" or paths must be specified

# Step 1: Create a new directory for the converted files
# mkdir -p Fluo_images/converted_tifs

# Step 2: Loop through all .ome.tif files in the source directory
for infile in Fluo_images/overlayed_ome_tif/*.ome.tif
do
  # Get the original filename without the extension
  base_name=$(basename "$infile" .ome.tif)
  
  # Set the output file path and name
  outfile="Fluo_images/warped_tif/${base_name}.tif"
  
  # Print a message to show progress
  echo "Converting ${infile}  ==>  ${outfile}"
  
  # Run the conversion command
  ./bftools/bfconvert -series 0 -bigtiff "$infile" "$outfile" # with the command '-series 0' I ask to convert only the first series of the ome.tif 
done

echo "Batch conversion complete!"


# single use command
# ./bftools/bfconvert -series 0 -bigtiff -overwrite Fluo_images/overlayed_ome_tif/blocco1.ome.tif Fluo_images/warped_tif/blocco1_try.tif
# try to no flatten the file 
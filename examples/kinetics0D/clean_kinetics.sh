#!/bin/sh

echo
echo "Cleaning logs and creating plot friendly data!"
echo "Usage: ./clean_logs.sh my_output_dir"
echo "Note: dont include trailing '/' on subdir or plotter will error!"
echo

# get current working directory for later
CWD=$(pwd)

if [ -z "$1" ]
then
    echo "Error: No argument supplied :("
    echo "Usage: provide log directory to this command"
    echo "IE:./clean_logs.sh my_kinetics_run_data_dir"
    exit 1
fi
FILEDIR=$1
mkdir -p $FILEDIR
echo "Cleaning logs from: $FILEDIR"

OUTFILE="./$FILEDIR/clean_data.csv"
echo "Writing to: $OUTFILE"

# remove old data if it exists
if [ -f $OUTFILE ]
then
   echo "Cleaning old data from $OUTFILE"
   rm $OUTFILE
else
   echo "No old data to clean from $OUTFILE ... continuing."
fi

#clean out old plots and data since it could have changed.
rm -f $FILEDIR/*.xda
rm -f $FILEDIR*.png

#will use number of dirs as the number of data points to collect (one point per proc)
NUMFILES=$( find . -name "*.xda" | wc -l )
echo
echo "Beginning processing of $NUMFILES total files."
echo

find .  -name "*.xda" -exec mv -t $FILEDIR {} +

cd $FILEDIR

# prepend header info to clean data csv file
HEADERSTR="Timestep, Temp, O, O2, O3"
echo $HEADERSTR >> $CWD/$OUTFILE

# use a counter to ensure each file gets hit.
proccessed_count=0

for file in $(find . -name '*.xda')
do
    # Make sure we have files to work on ...
    if [ -f $file ]
    then
        echo " Processing $file ..."
        PROCESSED_COUNT=$((PROCESSED_COUNT + 1 ))

        #first check file for unexpected errors
        ERR_DATA=$(grep -i "error" $file)
        #echo $ERR_DATA

        #clear string used to write data
        WRITESTR=""

        # first extract the time_step from filename
        # strip all (including _np) up to _npxx.log and then just trim to only get digits
        tmp=${file#*_np}
        NP_VAL=$(echo $tmp | tr -d -c 0-9)
        WRITESTR="$NP_VAL"

        # next extract data
        RAW_DATA=$(grep -i "Solution Vector" $file)

        # interpret grep result as string array to extract timing
        DATA_ARRAY=($RAW_DATA)

        #loop over requested array to extract data, use sed to seperate the csv into a bash useable array:
        for ITEM in ${DATA_ARRAY[@]}
        do
            if [ $ITEM != "#" ]
            then
                # append data to write string. if the str wasnt found we just write 0 (when we searchfor a fn not in logs)
                WRITESTR="$WRITESTR, $ITEM"
            else
              break
            fi
        done

        echo $WRITESTR >> $CWD/$OUTFILE

        #endif some log files found
    else
        echo "Error: No log files found in: $FILEDIR"
        exit 1
    fi
done #end file loop

# check to make sure all files were parsed
echo
if [ $PROCESSED_COUNT -eq $NUMFILES ]
then
    echo "Total processed files: $PROCESSED_COUNT"
    echo "All done! Exiting Cleaning Process..."
else
    echo "Error: Total processed files ($PROCESSED_COUNT) not equal to number of files ($NUMFILES) in logs directory. Some logs contain errors."
    exit 1
fi

# now run make_plots.py with the generated data
echo
echo "Running make_plots.py with generated clean_data.csv ...."
echo "Warning: this overwrites old *.png data in $FILEDIR"
echo
cd $CWD/$FILEDIR
rm -f *.png
echo "Raw data preview:"
cat clean_data.csv
echo

# pass in  timing type and filedir to use in plots
python $CWD/plot_species.py clean_data.csv $FILEDIR

echo "Plotting complete. Plots located at $CWD/$FILEDIR Exiting."
echo

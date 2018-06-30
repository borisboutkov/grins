#!/bin/sh

echo
echo "Cleaning logs and creating plot friendly data!"
echo "Usage: ./clean_kinetics.sh my_output_dir"
echo "Note: dont include trailing '/' on subdir or plotter will error!"
echo

# get current working directory for later
CWD=$(pwd)

if [ -z "$1" ]
then
    echo "Error: No argument supplied :("
    echo "Usage: provide log directory to this command"
    echo "IE:./clean_kinetics.sh my_kinetics_run_data_dir my_input.in"
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

if [ -z "$2" ]
then
    echo "Error: No input file argument supplied :("
    echo "Usage: provide log directory to this command"
    echo "IE:./clean_kinetics.sh my_kinetics_run_data_dir my_input.in"
    exit 1
fi
INPUTFILE=$2

# quite flexible :)
SPECIES_STR=$( grep -i "species = " $INPUTFILE | head -1 | cut -c 20- )
SPECIES_CLEAN=$( echo $SPECIES_STR | sed 's/.$//' )
SPECIES_HEADER="TS Temp $SPECIES_CLEAN"

echo "found header columns: $SPECIES_HEADER"


DELTA_T_STR=$( grep -i "delta_t = " $INPUTFILE)
#split on the single quote, take field 2
DELTA_T=$( echo $DELTA_T_STR | cut -d "'" -f2  )
echo "dt:" $DELTA_T

# prepend header info to clean data csv file
echo $SPECIES_HEADER >> $CWD/$OUTFILE

rm $FILEDIR".xda"

NUMFILES=$( find . -maxdepth 1 -name "*.xda" | wc -l )
echo
echo "Beginning processing of $NUMFILES total files."
echo
find . -maxdepth 1 -name "*.xda"  -exec mv -t $FILEDIR {} +

cd $FILEDIR

# prepend header info to clean data csv file
echo $HEADERSTR >> $CWD/$OUTFILE

# use a counter to ensure each file gets hit.
proccessed_count=0

#'alphabetize' by the timestep counter
for file in $(find . -name '*.xda' | ls -v *.xda)
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
        # strip all nonintegers from filename to only get digits
        TS_VAL=$(echo $file | cut -d "." -f2)
        if [[ -z "${TS_VAL// }" ]]
        then
           WRITESTR="0" #first file doesnt append a ts...
        else
            WRITESTR="$TS_VAL"
        fi

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
#echo "Raw data preview:"
#cat clean_data.csv
echo

# pass in timing type, filedir and delta_t to use in plots
python $CWD/plot_species.py clean_data.csv $FILEDIR $DELTA_T

echo "Plotting complete. Plots located at $CWD/$FILEDIR Exiting."
echo

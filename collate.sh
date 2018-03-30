# If the data files that contain the collated data do not already exists create them.
rm omegaConvergence.dat

touch omegaConvergence.dat

# Iterate through each set of results.
for results in "$@"*
do
    # For omega-convergence.
    # Get the values of omega from the input files and append it to collated data file.
    awk '/^(SOR-parameter:) /{printf $NF}' $results/input.txt >> omegaConvergence.dat

    # Put an empty column in the file.
    printf " " >> omegaConvergence.dat

    # Get the values of the order parameter from the results files and append it to collated data file.
    awk '/^(Number-of-iterations-until-convergence:) /{print $(NF)}' $results/results.txt >> omegaConvergence.dat

done

# Remove duplicate lines from files.
sort -u omegaConvergence.dat -o omegaConvergence.dat
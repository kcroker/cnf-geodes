#!/bin/bash

# Iterate over datasets
#
# WARNING: HiZ actually means 0.1Z_\odot, or low metallicity.  Hi(gh mass)Z.
# This is super awful, I apologize.
#
# Different variations of data set are commented out, so you can see how this works.
# Cuts, binning, etc can also be disabled.
# So once you cut and bin, you can just remake the graphs rapidly.
# (You have to manually comment the lines out though.)
#

# These sets are stable
for i in "NSNS"; do
    # Iterate over analyses
    #    for j in "control_timeuniform" "geodegeode" "geodens"; do
    for j in "geodegeode" "control"; do
    #for j in "control_timeuniform" "geodegeode_timeuniform" "geodens"; do
    
	# If the data set does not exist, do not try to make graphs for it
	# (not a combination we care about)
	if [ ! -f "results/A${i}_$j" ]; then
	    echo "Skipping set $i, analysis $j"
	    continue
	fi
	
	# # Redo cuts? (usually keep it commented unless you modified the cuts)
	cat "results/A${i}_$j" | ./horizoncut-stub.py > "results/A${i}_${j}_cut"

	# Do binning of cut results?
	top=0
	bins2d=0
	
	# Control populations have low statistics, so bin more coarsely so that something useful happens
	if [ "$j" == "control" ] || [ "$j" == "control_timeuniform" ]; then
	    bins2d=5
	else
	    bins2d=10
	fi

	# # Redo heatmap binnings?
	cat "results/A${i}_${j}_cut" | ./fast2dbin.py $bins2d > "results/A${i}_${j}_cut_2dbinned"

	# # Redo histogram binnings?
	# Eh, changed my mind.  Use the same scale for all graphs.
	if [ "$j" == "geodens" ]; then
	    top=150
	    cat "results/A${i}_${j}_cut" | ./fastbin-stub.py "4" "0 $top 5" > "results/A${i}_${j}_cut_1dbinned"
	else
	    top=150
	    cat "results/A${i}_${j}_cut" | ./fastbin-stub.py "4 5" "0 $top 5" > "results/A${i}_${j}_cut_1dbinned"
	fi
	
	cat heatmap.gpt | sed "s/__SET__/$i/g" | sed "s/__FIGURE__/$j/g" | sed "s/__BINS__/${bins2d}/g" | gnuplot
	pdflatex "chirpz$i-$j.tex"
	pdfcrop "chirpz$i-$j.pdf"

	# Only make the 1D histograms for non-control data, since control results
	# are included in these plots
	if [ ! "$j" == "control" ] && [ ! "$j" == "control_timeuniform" ]; then
	    
	    # This one assumes that the binned versions are alread there
	    cat prototype.gpt | sed "s/__FIGURE__/$j/g" | sed "s/__SET__/$i/g" | sed "s/__TOP__/$top/" | gnuplot
	    pdflatex "A$i-$j.tex"
	    pdfcrop "A$i-$j.pdf"
	fi
	
	# ######################################################################
	# # Do the CUTLESS graphs.  This was just for sanity checking things.
	# # Same code as above, with some hardcoded filename changes.
	# #
	# # NOTE: We keep the string "cut" inside the filenames
	# #       because that cut is hardcoded into the strings in the prototype plots
	# #
	# #######################################################################

	# # # Do predictions for 2x in redshift horizon cutoff (same shape)
	# #	cat "results/A${i}_$j" | ./horizoncut-double-stub.py > "results/A${i}_${j}_DOUBLE_cut"

       	# # # Redo heatmap binnings?
	# cat "results/A${i}_${j}" | ./fast2dbin.py $bins2d > "results/A${i}_${j}_RAW_cut_2dbinned"

	# # # Redo histogram binnings?
	# # Eh, changed my mind.  Use the same scale for all graphs.
	# if [ "$j" == "geodens" ]; then
	#     top=150
	#     cat "results/A${i}_${j}" | ./fastbin-stub.py "4" "0 $top 5" > "results/A${i}_${j}_RAW_cut_1dbinned"
	# else
	#     top=150
	#     cat "results/A${i}_${j}" | ./fastbin-stub.py "4 5" "0 $top 5" > "results/A${i}_${j}_RAW_cut_1dbinned"
	# fi
	
	# cat heatmap.gpt | sed "s/__SET__/$i/g" | sed "s/__FIGURE__/${j}_RAW/g" | sed "s/__BINS__/${bins2d}/g" | gnuplot
	# pdflatex "chirpz$i-${j}_RAW.tex"
	# pdfcrop "chirpz$i-${j}_RAW.pdf"

	# # Only make the 1D histograms for non-control data, since control results
	# # are included in these plots
	# if [ ! "$j" == "control" ] && [ ! "$j" == "control_timeuniform" ]; then
	    
	#     # This one assumes that the binned versions are alread there
	#     cat prototype.gpt | sed "s/__FIGURE__/${j}_RAW/g" | sed "s/__SET__/$i/g" | sed "s/__TOP__/$top/" | gnuplot
	#     pdflatex "A$i-${j}_RAW.tex"
	#     pdfcrop "A$i-${j}_RAW.pdf"
	# fi

    done
done

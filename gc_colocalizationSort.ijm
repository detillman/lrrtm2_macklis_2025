// ask user to select a folder
dir = getDir("Choose a Directory");

// make temporary directories
File.makeDirectory(dir + "deleteME");
File.makeDirectory(dir + "deleteME//first");
File.makeDirectory(dir + "deleteME//second");

// open rois
open(dir + "Analysis_RoiSet.zip");

// open results and get headings
open(dir + "Analysis_Results.csv");
IJ.renameResults("Analysis_Results.csv", "Results");
headings = split(String.getResultsHeadings);

// get column labels
first_label = replace(headings[5], "Bkgd_Corr_Intensity_", "");
second_label = replace(headings[10], "Bkgd_Corr_Intensity_", "");

// make empty columns
first_column = newArray(nResults);
second_column = newArray(nResults);

// get background corrected intensities
for (i=0; i<nResults; i++) {
	first_column[i] = getResult(headings[5], i);
	second_column[i] = getResult(headings[10], i);
}

// get rank indices
first_positions = Array.rankPositions(first_column);
second_positions = Array.rankPositions(second_column);

// get ranks
first_ranks = Array.rankPositions(first_positions);
second_ranks = Array.rankPositions(second_positions);
close("Results");

// open image and split
open(dir + "Analysis_Gallery.zip");
selectWindow("Analysis_Gallery.tif");
run("Split Channels");

// delete excluded rois
for (i = roiManager("Count") - 1; i >= 0; i--) {
	roiManager("Select", i);
	strokeColor = Roi.getStrokeColor;
	if(strokeColor == "red") {
		roiManager("delete");
	}
}

// rank particles
for (i = roiManager("Count") - 1; i >= 0; i--) {
	selectWindow("C1-Analysis_Gallery.tif");
    roiManager("Select", i);
	detection = Roi.getName;
	run("Duplicate...", "title=first_crop");
    run("Crop");
	first_name = IJ.pad(first_ranks[i], 6) + "_" + round(first_column[i]) + "_" + detection + "_" + first_label;
	saveAs("Tiff", dir + "deleteME//first//" + first_name);
	close();

	selectWindow("C2-Analysis_Gallery.tif");
    roiManager("Select", i);
	detection = Roi.getName;
	run("Duplicate...", "title=second_crop");
 	run("Crop");
	second_name = IJ.pad(second_ranks[i], 6) + "_" + round(second_column[i]) + "_" + detection + "_" + second_label;
	saveAs("Tiff", dir + "deleteME//second//" + second_name);
	close();
}

// close images
close("C1-Analysis_Gallery.tif");
close("C2-Analysis_Gallery.tif");

// get montage dimensions
nCol = floor(sqrt(roiManager("Count")))
nRow = round(roiManager("Count")/nCol+0.5)

// make first montage
File.openSequence(dir + "deleteME//first");
saveAs("zip", dir + first_label + "_stack");
run("Make Montage...", "columns="+nCol+" rows="+nRow+" scale=1");
rename("first_gallery");
saveAs("zip", dir + first_label + "_gallery");

// make second montage
File.openSequence(dir + "deleteME//second");
saveAs("zip", dir + second_label + "_stack");
run("Make Montage...", "columns="+nCol+" rows="+nRow+" scale=1");
rename("second_gallery");
saveAs("zip", dir + second_label + "_gallery");
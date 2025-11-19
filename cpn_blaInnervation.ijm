// rotate so that somata to be in top-left
// rotate so that midline is vertical and crop edges
// draw bla_contra and bla_ipsi first draft based on dapi
// delete dapi
// auto-adjust contrast
// update bla_contra and bla_ipsi rois
// draw ipsi_axons and ipsi_cortex rois
// save image as ".tif" and rois as "_ROI.zip"
// horizontally crop where cortex is widest
// vertically crop at midline
// save as "_contra.tif" and "_ipsi.tif"

//// 4 ROIs
// bla_contra = dapi-based contralateral bla, modified based on fluor
// bla_ipsi = dapi-based ipsialateral bla, modified based on fluor
// ipsi_axons = start at midline, trace ipsilateral axons to bla
// ipsi_cortex = start at midline, trace ipsilateral axons to bla, then do cortical surface

// select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);

// parameters
numRows = 1
numCols = 500
total_boxes = numRows * numCols;
run("Set Measurements...", "area mean min redirect=None decimal=4");

// run macro
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = dir + list[i];
	name = File.getName(path);
	strip_name = replace(name, "_ROI.zip", "");
	if (endsWith(path, "zip")) {
		// open rois in non-batch mode
		setBatchMode(false);
		open(dir + strip_name + "_ROI.zip");
		roiManager("Sort");
		
		// enter batch mode
		setBatchMode(true);

		// open coronal section
		open(strip_name + ".tif");
		rename("coronal_thresh");
		run("Duplicate...", "title=coronal");
		selectImage("coronal");
		run("32-bit");
		getPixelSize(unit, pixelWidth, pixelHeight);

		// make arrays
		name_results = newArray(total_boxes);
		box_results = newArray(total_boxes);
		pix_results = newArray(total_boxes);
		dist_results = newArray(total_boxes);
		mean_results = newArray(total_boxes);
		threshMean_results = newArray(total_boxes);
		threshFrac_results = newArray(total_boxes);
		threshArea_results = newArray(total_boxes);
		threshNumber_results = newArray(total_boxes);
		height_results = newArray(total_boxes);
		contraBLAarea_results = newArray(total_boxes);
		contraBLAmean_results = newArray(total_boxes);
		contraBLAtotal_results = newArray(total_boxes);
		contraBLAdensity_results = newArray(total_boxes);
		contraBLAmeanThresh_results = newArray(total_boxes);
		contraBLAtotalThresh_results = newArray(total_boxes);
		contraBLAdensityThresh_results = newArray(total_boxes);
		ipsiBLAarea_results = newArray(total_boxes);
		ipsiBLAmean_results = newArray(total_boxes);
		ipsiBLAtotal_results = newArray(total_boxes);
		ipsiBLAdensity_results = newArray(total_boxes);
		ipsiBLAmeanThresh_results = newArray(total_boxes);
		ipsiBLAtotalThresh_results = newArray(total_boxes);
		ipsiBLAdensityThresh_results = newArray(total_boxes);

		// quick process
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");

		// measure contralateral BLA
		roiManager("Select", 0); // bla_contra
		run("Measure");
		contraBLAarea = getResult("Area", 0);
		contraBLAmean = getResult("Mean", 0);
		close("Results");
		roiManager("Deselect");

		// measure ipsilateral BLA
		roiManager("Select", 1); // bla_ipsi
		run("Measure");
		ipsiBLAarea = getResult("Area", 0);
		ipsiBLAmean = getResult("Mean", 0);
		close("Results");
		roiManager("Deselect");

		// turn non-ROI into NaN
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 3); // ipsi_cortex
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// get coordinates for ipsi lines
		roiManager("Select", 2); // ipsi_axons
		getSelectionCoordinates(x_ipsiSmall, y_ipsiSmall);
		run("Scale... ", "x=1.50 y=1.50 centered");
		run("Fit Spline");
		getSelectionCoordinates(x_ipsiBig, y_ipsiBig);

		// remove contra values
		while(x_ipsiBig[0] > x_ipsiSmall[0]) {
			x_ipsiBig = Array.deleteIndex(x_ipsiBig, 0);
		}
		y_ipsiBig = Array.slice(y_ipsiBig, y_ipsiBig.length-x_ipsiBig.length);

		// make final ipsi line
		makeSelection(6, x_ipsiBig, y_ipsiBig);
		Roi.setStrokeWidth(4500);
		run("Straighten...");
		rename("ipsi_straight");

		// calculate box dimensions
		boxWidth = getWidth() / numCols;
		boxHeight = getHeight() / numRows;

		// loop through rows and columns to calculate box intensities
		for (row = 0; row < numRows; row++) {
  			for (col = 0; col < numCols; col++) {
    			
      			// calculate box coordinates
      			boxX = 0 + col * boxWidth;
      			boxY = 0 + row * boxHeight;

      			// define ROI for current box
      			makeRectangle(boxX, boxY, boxWidth, boxHeight);

      			// calculate box intensity
      			run("Measure");
      
			}			
		}
		close("ipsi_straight");
		close("coronal");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			box_results[box] = box;
			pix_results[box] = boxWidth * (box + 1);
			dist_results[box] = boxWidth * (box + 1) * pixelWidth;
			mean_results[box] = getResult("Mean", box);
			contraBLAarea_results[box] = contraBLAarea;
			contraBLAmean_results[box] = contraBLAmean;
			contraBLAtotal_results[box] = contraBLAmean * contraBLAarea;
			contraBLAdensity_results[box] = contraBLAmean / contraBLAarea;
			ipsiBLAarea_results[box] = ipsiBLAarea;
			ipsiBLAmean_results[box] = ipsiBLAmean;
			ipsiBLAtotal_results[box] = ipsiBLAmean * ipsiBLAarea;
			ipsiBLAdensity_results[box] = ipsiBLAmean / ipsiBLAarea;
		}
		close("Results");

		// get and set threshold
		selectImage("coronal_thresh");
		roiManager("Select", 3); // ipsi_cortex
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
		// measure contralateral BLA
		roiManager("Select", 0); // bla_contra
		run("Measure");
		contraBLAareaThresh = getResult("Area", 0);
		contraBLAmeanThresh = getResult("Mean", 0);
		close("Results");
		roiManager("Deselect");

		// measure ipsilateral BLA
		roiManager("Select", 1); // bla_ipsi
		run("Measure");
		ipsiBLAareaThresh = getResult("Area", 0);
		ipsiBLAmeanThresh = getResult("Mean", 0);
		close("Results");
		roiManager("Deselect");
		
		// set non-roi pixels to NaN
		run("32-bit");
		roiManager("Select", 3); // ipsi_cortex
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// make final ipsi line
		makeSelection(6, x_ipsiBig, y_ipsiBig);
		Roi.setStrokeWidth(4500);
		run("Straighten...");
		rename("ipsi_straight");

		// calculate box dimensions
		boxWidth = getWidth() / numCols;
		boxHeight = getHeight() / numRows;

		// loop through rows and columns to calculate thresholded intensities
		for (row = 0; row < numRows; row++) {
		  	for (col = 0; col < numCols; col++) {
    			
     		 	// calculate box coordinates
      			boxX = 0 + col * boxWidth;
      			boxY = 0 + row * boxHeight;

      			// define ROI for current box
      			makeRectangle(boxX, boxY, boxWidth, boxHeight);

      			// calculate thresholded intensity of box
      			run("Measure");
			}			
		}
		close("ipsi_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			threshMean_results[box] = getResult("Mean", box);
			threshFrac_results[box] = threshMean_results[box] / 255;
			threshArea_results[box] = getResult("Area", box);
			threshNumber_results[box] = threshArea_results[box] * threshFrac_results[box];
			height_results[box] = threshArea_results[box] / boxWidth;
			contraBLAmeanThresh_results[box] = contraBLAmeanThresh;
			contraBLAtotalThresh_results[box] = contraBLAmeanThresh * contraBLAarea;
			contraBLAdensityThresh_results[box] = contraBLAmeanThresh / contraBLAarea;
			ipsiBLAmeanThresh_results[box] = ipsiBLAmeanThresh;
			ipsiBLAtotalThresh_results[box] = ipsiBLAmeanThresh * ipsiBLAarea;
			ipsiBLAdensityThresh_results[box] = ipsiBLAmeanThresh / ipsiBLAarea;
		}
		close("Results");
		close("coronal_thresh");

		// store arrays in table and save
		Table.create("coronalTable");
		Table.setColumn("file", name_results);
		Table.setColumn("box", box_results);
		Table.setColumn("distance_pixels", pix_results);
		Table.setColumn("distance_microns", dist_results);
		Table.setColumn("height_microns", height_results);
		Table.setColumn("mean", mean_results);
		Table.setColumn("thresholdMean", threshMean_results);
		Table.setColumn("thresholdFraction", threshFrac_results);
		Table.setColumn("thresholdArea", threshArea_results);
		Table.setColumn("thresholdNumber", threshNumber_results);
		Table.setColumn("contraBLAarea", contraBLAarea_results);
		Table.setColumn("contraBLAmean", contraBLAmean_results);
		Table.setColumn("contraBLAtotal", contraBLAtotal_results);
		Table.setColumn("contraBLAdensity", contraBLAdensity_results);
		Table.setColumn("contraBLAmean_thresh", contraBLAmeanThresh_results);
		Table.setColumn("contraBLAtotal_thresh", contraBLAtotalThresh_results);
		Table.setColumn("contraBLAdensity_thresh", contraBLAdensityThresh_results);
		Table.setColumn("ipsiBLAarea", ipsiBLAarea_results);
		Table.setColumn("ipsiBLAmean", ipsiBLAmean_results);
		Table.setColumn("ipsiBLAtotal", ipsiBLAtotal_results);
		Table.setColumn("ipsiBLAdensity", ipsiBLAdensity_results);
		Table.setColumn("ipsiBLAmean_thresh", ipsiBLAmeanThresh_results);
		Table.setColumn("ipsiBLAtotal_thresh", ipsiBLAtotalThresh_results);
		Table.setColumn("ipsiBLAdensity_thresh", ipsiBLAdensityThresh_results);
		Table.save(dir + "/" + strip_name + ".csv");
		close("coronalTable");

		// close ROI manager
		close("ROI Manager");
	}
}

beep();
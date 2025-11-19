//// whole coronal tif with dapi
//// ipsi and contra tifs without dapi
//// ROIs
// contra_cortex = polygon of superficial contra layers
// contra_axons = segmented line from midline to distal axon curve
// ipsi_cortex = polygon of superficial ipsi layers
// ipsi_axons = segmented line from midline to proximal axon curve

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
		
		// open ipsi image
		open(strip_name + "_ipsi.tif");
		rename("ipsi");
		run("32-bit");
		getPixelSize(unit, pixelWidth, pixelHeight);

		// make ipsi arrays
		name_results = newArray(total_boxes);
		feature_results = newArray(total_boxes);
		box_results = newArray(total_boxes);
		pix_results = newArray(total_boxes);
		dist_results = newArray(total_boxes);
		mean_results = newArray(total_boxes);
		threshMean_results = newArray(total_boxes);
		threshFrac_results = newArray(total_boxes);
		threshArea_results = newArray(total_boxes);
		threshNumber_results = newArray(total_boxes);
		height_results = newArray(total_boxes);

		// quick process
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");

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

		// get and set threshold
		selectImage("ipsi");
		roiManager("Select", 3); // ipsi_cortex
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
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
			name_results[box] = strip_name;
			feature_results[box] = "ipsi";
			box_results[box] = box;
			pix_results[box] = boxWidth * (box + 1);
			dist_results[box] = boxWidth * (box + 1) * pixelWidth;
			mean_results[box] = getResult("Mean", box);
			threshMean_results[box] = getResult("Mean", box + total_boxes);
			threshFrac_results[box] = threshMean_results[box] / 255;
			threshArea_results[box] = getResult("Area", box + total_boxes);
			threshNumber_results[box] = threshArea_results[box] * threshFrac_results[box];
			height_results[box] = threshArea_results[box] / boxWidth;
		}
		close("Results");
		close("ipsi");

		// store arrays in table and save
		Table.create("ipsiTable");
		Table.setColumn("file", name_results);
		Table.setColumn("feature", feature_results);
		Table.setColumn("box", box_results);
		Table.setColumn("distance_pixels", pix_results);
		Table.setColumn("distance_microns", dist_results);
		Table.setColumn("height_microns", height_results);
		Table.setColumn("mean", mean_results);
		Table.setColumn("thresholdMean", threshMean_results);
		Table.setColumn("thresholdFraction", threshFrac_results);
		Table.setColumn("thresholdArea", threshArea_results);
		Table.setColumn("thresholdNumber", threshNumber_results);
		Table.save(dir + "/" + strip_name + "_ipsi.csv");
		close("ipsiTable");


		// open contra image
		open(strip_name + "_contra.tif");
		rename("contra");
		run("32-bit");
		getPixelSize(unit, pixelWidth, pixelHeight);

		// make contra arrays
		name_results = newArray(total_boxes);
		feature_results = newArray(total_boxes);
		box_results = newArray(total_boxes);
		pix_results = newArray(total_boxes);
		dist_results = newArray(total_boxes);
		mean_results = newArray(total_boxes);
		threshMean_results = newArray(total_boxes);
		threshFrac_results = newArray(total_boxes);
		threshArea_results = newArray(total_boxes);
		threshNumber_results = newArray(total_boxes);
		height_results = newArray(total_boxes);

		// quick process
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");

		// turn non-ROI into NaN
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 1); // contra_cortex
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// get coordinates for contra lines
		roiManager("Select", 0); // contra_axons
		getSelectionCoordinates(x_contraSmall, y_contraSmall);
		run("Scale... ", "x=1.50 y=1.50 centered");
		run("Fit Spline");
		getSelectionCoordinates(x_contraBig, y_contraBig);

		// remove ipsilateral values
		while(x_contraBig[0] < x_contraSmall[0]) {
			x_contraBig = Array.deleteIndex(x_contraBig, 0);
		}
		y_contraBig = Array.slice(y_contraBig, y_contraBig.length-x_contraBig.length);

		// make final contra line
		makeSelection(6, x_contraBig, y_contraBig);
		Roi.setStrokeWidth(4500);
		run("Straighten...");
		rename("contra_straight");

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
		close("contra_straight");

		// get and set threshold
		selectImage("contra");
		roiManager("Select", 1); // contra_cortex
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
		// set non-roi pixels to NaN
		run("32-bit");
		roiManager("Select", 1); // contra_cortex
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// make final contra line
		makeSelection(6, x_contraBig, y_contraBig);
		Roi.setStrokeWidth(4500);
		run("Straighten...");
		rename("contra_straight");

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
		close("contra_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "contra";
			box_results[box] = box;
			pix_results[box] = boxWidth * (box + 1);
			dist_results[box] = boxWidth * (box + 1) * pixelWidth;
			mean_results[box] = getResult("Mean", box);
			threshMean_results[box] = getResult("Mean", box + total_boxes);
			threshFrac_results[box] = threshMean_results[box] / 255;
			threshArea_results[box] = getResult("Area", box + total_boxes);
			threshNumber_results[box] = threshArea_results[box] * threshFrac_results[box];
			height_results[box] = threshArea_results[box] / boxWidth;
		}
		close("Results");
		close("contra");

		// store arrays in table and save
		Table.create("contraTable");
		Table.setColumn("file", name_results);
		Table.setColumn("feature", feature_results);
		Table.setColumn("box", box_results);
		Table.setColumn("distance_pixels", pix_results);
		Table.setColumn("distance_microns", dist_results);
		Table.setColumn("height_microns", height_results);
		Table.setColumn("mean", mean_results);
		Table.setColumn("thresholdMean", threshMean_results);
		Table.setColumn("thresholdFraction", threshFrac_results);
		Table.setColumn("thresholdArea", threshArea_results);
		Table.setColumn("thresholdNumber", threshNumber_results);
		Table.save(dir + "/" + strip_name + "_contra.csv");
		close("contraTable");
		
		// close ROI manager
		close("ROI Manager");
	}
}

beep();
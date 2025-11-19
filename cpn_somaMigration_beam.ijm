//// WORKFLOW
/// rotate image so that cpn somata are in top-left quadrant and midline is vertical
/// crop rectangle using cortical surface as boundary and delete dapi channel
/// auto-adjust contrast in red and green
/// save as tif
/// make composite
/// draw polygon "threshold" roi around merged fluoresence
/// draw segmented "straight" roi along merged fluoresence
/// save rois and close image

// select directory and get files
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);
run("Set Measurements...", "area mean min redirect=None decimal=4");

// set parameters
numRows = 1;
numCols = 500;
total_boxes = numRows * numCols;

// run macro
for (i=0; i<list.length; i++) {
	
	// running through each file
	showProgress(i+1, list.length);
	path = dir + list[i];
	name = File.getName(path);
	strip_name = replace(name, ".tif", "");
	if (endsWith(path, "tif")) {
		// open rois in non-batch mode
		setBatchMode(false);
		open(dir + strip_name + "_ROI.zip");
		roiManager("Sort");
		
		// enter batch mode and open image
		setBatchMode(true);
		open(strip_name + ".tif");
		getPixelSize(unit, pixelWidth, pixelHeight);
		rename("original");
		
		// process green channel
		run("Duplicate...", "title=green duplicate channels=1");
		selectImage("green");
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");
		
		// duplicate for threshold
		run("Duplicate...", "title=green_threshold");
		selectImage("green_threshold");
		
		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 1); // threshold
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// set and get threhsold
		selectImage("green_threshold");
		roiManager("Select", 1); // threshold
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		close("green_threshold");
		
		// analyze green
		selectImage("green");
		
		// make arrays
		name_results = newArray(total_boxes);
		color_results = newArray(total_boxes);
		box_results = newArray(total_boxes);
		pix_results = newArray(total_boxes);
		dist_results = newArray(total_boxes);
		mean_results = newArray(total_boxes);
		threshMean_results = newArray(total_boxes);
		threshFrac_results = newArray(total_boxes);
		threshArea_results = newArray(total_boxes);
		threshNumber_results = newArray(total_boxes);
		height_results = newArray(total_boxes);
		
		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 1); // threshold
		run("Make Inverse");
		changeValues(0, max_value, NaN);
		
		// straighten straight
		roiManager("Select", 0); // straight
		Roi.setStrokeWidth(1000);
		run("Straighten...");
		rename("green_straight");
		
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
		
		close("green_straight");
		
		// threshold green
		selectImage("green");
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
		// set non-roi pixels to NaN
		run("32-bit");
		roiManager("Select", 1); // threshold
		run("Make Inverse");
		changeValues(0, max_value, NaN);
		
		// straighten working
		roiManager("Select", 0); // straight
		Roi.setStrokeWidth(1000); ////
		run("Straighten...");
		rename("green_straight");
		
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
		
		close("green_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			color_results[box] = "green";
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
		close("green");
		
		// store arrays in table and save
		Table.create("greenTable");
		Table.setColumn("file", name_results);
		Table.setColumn("color", color_results);
		Table.setColumn("box", box_results);
		Table.setColumn("distance_pixels", pix_results);
		Table.setColumn("distance_microns", dist_results);
		Table.setColumn("height_microns", height_results);
		Table.setColumn("mean", mean_results);
		Table.setColumn("thresholdMean", threshMean_results);
		Table.setColumn("thresholdFraction", threshFrac_results);
		Table.setColumn("thresholdArea", threshArea_results);
		Table.setColumn("thresholdNumber", threshNumber_results);
		Table.save(dir + "/" + strip_name + "_greenDistribution.csv");
		close("greenTable");
		
		
		// process red channel
		run("Duplicate...", "title=red duplicate channels=2");
		selectImage("red");
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");
		
		// duplicate for threshold
		run("Duplicate...", "title=red_threshold");
		selectImage("red_threshold");
		
		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 1); // threshold
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// set and get threhsold
		selectImage("red_threshold");
		roiManager("Select", 1); // threshold
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		close("red_threshold");
		
		// analyze red
		selectImage("red");
		
		// make arrays
		name_results = newArray(total_boxes);
		color_results = newArray(total_boxes);
		box_results = newArray(total_boxes);
		pix_results = newArray(total_boxes);
		dist_results = newArray(total_boxes);
		mean_results = newArray(total_boxes);
		threshMean_results = newArray(total_boxes);
		threshFrac_results = newArray(total_boxes);
		threshArea_results = newArray(total_boxes);
		threshNumber_results = newArray(total_boxes);
		height_results = newArray(total_boxes);
		
		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 1); // threshold
		run("Make Inverse");
		changeValues(0, max_value, NaN);
		
		// straighten straight
		roiManager("Select", 0); // straight
		Roi.setStrokeWidth(1000);
		run("Straighten...");
		rename("red_straight");
		
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
		
		close("red_straight");
		
		// threshold red
		selectImage("red");
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
		// set non-roi pixels to NaN
		run("32-bit");
		roiManager("Select", 1); // threshold
		run("Make Inverse");
		changeValues(0, max_value, NaN);
		
		// straighten working
		roiManager("Select", 0); // straight
		Roi.setStrokeWidth(1000); ////
		run("Straighten...");
		rename("red_straight");
		
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
		
		close("red_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			color_results[box] = "red";
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
		close("red");
		
		// store arrays in table and save
		Table.create("redTable");
		Table.setColumn("file", name_results);
		Table.setColumn("color", color_results);
		Table.setColumn("box", box_results);
		Table.setColumn("distance_pixels", pix_results);
		Table.setColumn("distance_microns", dist_results);
		Table.setColumn("height_microns", height_results);
		Table.setColumn("mean", mean_results);
		Table.setColumn("thresholdMean", threshMean_results);
		Table.setColumn("thresholdFraction", threshFrac_results);
		Table.setColumn("thresholdArea", threshArea_results);
		Table.setColumn("thresholdNumber", threshNumber_results);
		Table.save(dir + "/" + strip_name + "_redDistribution.csv");
		close("redTable");
		
		close("ROI Manager");
		close("original");
	
	}
}

beep();
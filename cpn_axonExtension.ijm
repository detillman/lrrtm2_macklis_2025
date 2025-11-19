//// **DONT'T USE NON-SOMA THRESH**
//// whole coronal tif without dapi
//// ROIs
// somata = polygon of iue, all layers
// somata_line = ipsilateral cortical surface, start at dorsal
// proximal = polygon from midline to lateral edge of iue
// distal = polygon from midline to distal axon tip
// distal_line = segmented line from midline to distal axon tip
// proximal_line = segmented line from lateral edge of iue to midline

//select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);
run("Set Measurements...", "area mean min redirect=None decimal=4");

// run macro
for (i=0; i<list.length; i++) {
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
		open(path);
		
		// get parameters
		numRows = 1;
		numCols = 500;
		total_boxes = numRows * numCols;
		getPixelSize(unit, pixelWidth, pixelHeight);

		// quick process
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");
		rename("somata");

		// duplicates for threshold and axons
		run("Duplicate...", "title=threshold");
		run("Duplicate...", "title=proximal");
		run("Duplicate...", "title=distal");

		// find threshold
		selectImage("threshold");

		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 4); // somata		
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// set and get threhsold
		selectImage("threshold");
		roiManager("Select", 4); // somata
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		close("threshold");

		//// analyze somata
		selectImage("somata");

		// make arrays
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

		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 4); // somata
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// straighten somata
		roiManager("Select", 5); // somata line
		run("Scale... ", "x=0.75 y=0.75 centered");
		Roi.setStrokeWidth(2500);
		run("Straighten...");
		rename("somata_straight");

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
		close("somata_straight");

		// threshold image
		selectImage("somata");
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
		// set non-roi pixels to NaN
		run("32-bit");
		roiManager("Select", 4); // somata
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// straighten somata
		roiManager("Select", 5); // somata line
		Roi.setStrokeWidth(1);
		run("Scale... ", "x=0.75 y=0.75 centered");
		Roi.setStrokeWidth(2500);
		run("Straighten...");
		rename("somata_straight");

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
		close("somata_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "somata";
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
		close("somata");

		// store arrays in table and save
		Table.create("somataTable");
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
		Table.save(dir + "/" + strip_name + "_somata.csv");
		close("somataTable");


		//// analyze proximal
		selectImage("proximal");

		// make arrays
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

		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 2); // proximal
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// straighten proximal
		roiManager("Select", 3); // proximal line
		Roi.setStrokeWidth(1000);
		run("Straighten...");
		rename("proximal_straight");

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
		close("proximal_straight");

		// threshold image
		selectImage("proximal");
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
		// set non-roi pixels to NaN
		run("32-bit");
		roiManager("Select", 2); // proximal
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// straighten proximal
		roiManager("Select", 3); // proximal line
		Roi.setStrokeWidth(1000);
		run("Straighten...");
		rename("proximal_straight");

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
		close("proximal_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "proximal";
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
		close("proximal");

		// store arrays in table and save
		Table.create("proximalTable");
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
		Table.save(dir + "/" + strip_name + "_proximal.csv");
		close("proximalTable");


		//// analyze distal
		selectImage("distal");

		// make arrays
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

		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 0); // distal
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// straighten distal
		roiManager("Select", 1); // distal line
		Roi.setStrokeWidth(1000);
		run("Straighten...");
		rename("distal_straight");

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
		close("distal_straight");

		// threshold image
		selectImage("distal");
		setThreshold(lowerThresh, upperThresh, "raw");
		run("Convert to Mask");
		
		// set non-roi pixels to NaN
		run("32-bit");
		roiManager("Select", 0); // distal
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// straighten distal
		roiManager("Select", 1); // distal line
		Roi.setStrokeWidth(1000);
		run("Straighten...");
		rename("distal_straight");

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
		close("distal_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "distal";
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
		close("distal");

		// store arrays in table and save
		Table.create("distalTable");
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
		Table.save(dir + "/" + strip_name + "_distal.csv");
		close("distalTable");
		close("ROI Manager");
	}
}

beep();
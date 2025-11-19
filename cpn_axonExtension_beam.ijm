//// **DONT'T USE NON-SOMA THRESH**
/// flip, rotate, and crop in dapi, then delete dapi
/// auto-adjust green and red, then save as "_split.tif"
/// make composite, then draw and save six ROIs
// somata = polygon of iue, all layers
// somata_line = ipsilateral cortical surface segemented line, midline to dorsal-most point (not midline)
// proximal = polygon from midline to lateral edge of iue
// distal = polygon from midline to distal axon tip
// dsital_line = segmented line from midline to distal axon tip
// proximal_line = segmented line from lateral edge of iue to midline
/// save composite as "_merge.jpg"

//select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);
run("Set Measurements...", "area mean min redirect=None decimal=4");

// parameters
numRows = 1;
numCols = 500;
total_boxes = numRows * numCols;

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
		open(strip_name + ".tif");
		getPixelSize(unit, pixelWidth, pixelHeight);
		rename("original");
		
		// quickly process green channel
		run("Duplicate...", "title=green_somata duplicate channels=1");
		selectImage("green_somata");
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");

		// duplicates for threshold and axons
		run("Duplicate...", "title=green_threshold");
		run("Duplicate...", "title=green_proximal");
		run("Duplicate...", "title=green_distal");

		// find threshold
		selectImage("green_threshold");

		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 4); // somata		
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// set and get threhsold
		selectImage("green_threshold");
		roiManager("Select", 4); // somata
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		close("green_threshold");

		//// analyze somata
		selectImage("green_somata");
		
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
		Roi.setStrokeWidth(1500);
		run("Straighten...");
		rename("green_somata_straight");

		// calculate box dimensions
		boxWidth = getWidth() / numCols;
		boxHeight = getHeight() / numRows;

		// loop through rows and columns to calculate coritcal position box intensities
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
		
		// calculate box dimensions for cortical layer analysis
		boxWidth = getWidth() / numRows; // numRows = 1
		boxHeight = getHeight() / numCols; // numCols = 500
		
		// loop through columns and rows to calculate cortical layer box intensities
		for (col = 0; col < numRows; col++) {
  			for (row = 0; row < numCols; row++) {
    			
 		     	// calculate box coordinates
 		     	boxX = 0 + col * boxWidth;
  		    	boxY = 0 + row * boxHeight;

   		   		// define ROI for current box
   		   		makeRectangle(boxX, boxY, boxWidth, boxHeight);

    		  	// calculate box intensity
    		  	run("Measure");
			}			
		}		
		close("green_somata_straight");

		// threshold image
		selectImage("green_somata");
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
		Roi.setStrokeWidth(1500);
		run("Straighten...");
		rename("green_somata_straight");

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
		
		// calculate box dimensions for cortical layer analysis
		boxWidth_cl = getWidth() / numRows; // numRows = 1
		boxHeight_cl = getHeight() / numCols; // numCols = 500
		
		// loop through columns and rows to calculate cortical layer box intensities
		for (col = 0; col < numRows; col++) {
  			for (row = 0; row < numCols; row++) {
    			
 		     	// calculate box coordinates
 		     	boxX = 0 + col * boxWidth_cl;
  		    	boxY = 0 + row * boxHeight_cl;

   		   		// define ROI for current box
   		   		makeRectangle(boxX, boxY, boxWidth_cl, boxHeight_cl);

    		  	// calculate box intensity
    		  	run("Measure");
			}			
		}		
		close("green_somata_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "green_somata";
			box_results[box] = box;
			pix_results[box] = boxWidth * (box + 1);
			dist_results[box] = boxWidth * (box + 1) * pixelWidth;
			mean_results[box] = getResult("Mean", box);
			threshMean_results[box] = getResult("Mean", box + total_boxes*2);
			threshFrac_results[box] = threshMean_results[box] / 255;
			threshArea_results[box] = getResult("Area", box + total_boxes*2);
			threshNumber_results[box] = threshArea_results[box] * threshFrac_results[box];
			height_results[box] = threshArea_results[box] / boxWidth;
		}

		// store arrays in table and save
		Table.create("green_somataTable");
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
		Table.save(dir + "/" + strip_name + "_green_somata.csv");
		close("green_somataTable");
		
		// make cortical layer arrays
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
		
		// store coritcal layer values in cortical layer arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "green_somataLayer";
			box_results[box] = box;
			pix_results[box] = boxHeight_cl * (box + 1);
			dist_results[box] = boxHeight_cl * (box + 1) * pixelHeight;
			mean_results[box] = getResult("Mean", box + total_boxes);
			threshMean_results[box] = getResult("Mean", box + total_boxes*3);
			threshFrac_results[box] = threshMean_results[box] / 255;
			threshArea_results[box] = getResult("Area", box + total_boxes*3);
			threshNumber_results[box] = threshArea_results[box] * threshFrac_results[box];
			height_results[box] = threshArea_results[box] / boxHeight_cl;
		}
		close("Results");
		close("green_somata");

		// store arrays in table and save
		Table.create("green_somataLayerTable");
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
		Table.save(dir + "/" + strip_name + "_green_somataLayer.csv");
		close("green_somataLayerTable");
		
		
		//// analyze proximal
		selectImage("green_proximal");

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
		rename("green_proximal_straight");

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
		close("green_proximal_straight");

		// threshold image
		selectImage("green_proximal");
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
		rename("green_proximal_straight");

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
		close("green_proximal_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "green_proximal";
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
		close("green_proximal");

		// store arrays in table and save
		Table.create("green_proximalTable");
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
		Table.save(dir + "/" + strip_name + "_green_proximal.csv");
		close("green_proximalTable");


		//// analyze distal
		selectImage("green_distal");

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
		rename("green_distal_straight");

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
		close("green_distal_straight");

		// threshold image
		selectImage("green_distal");
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
		rename("green_distal_straight");

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
		close("green_distal_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "green_distal";
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
		close("green_distal");

		// store arrays in table and save
		Table.create("green_distalTable");
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
		Table.save(dir + "/" + strip_name + "_green_distal.csv");
		close("green_distalTable");
		
		// quickly process red channel
		selectImage("original");
		run("Duplicate...", "title=red_somata duplicate channels=2");
		selectImage("red_somata");
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");

		// duplicates for threshold and axons
		run("Duplicate...", "title=red_threshold");
		run("Duplicate...", "title=red_proximal");
		run("Duplicate...", "title=red_distal");

		// find threshold
		selectImage("red_threshold");

		// turn non-ROI into NaN
		run("32-bit");
		run("Measure");
		max_value = getResult("Max", 0);
		close("Results");
		roiManager("Select", 4); // somata		
		run("Make Inverse");
		changeValues(0, max_value, NaN);

		// set and get threhsold
		selectImage("red_threshold");
		roiManager("Select", 4); // somata
		roiManager("Deselect");
		setAutoThreshold("Default dark");
		getThreshold(lowerThresh, upperThresh);
		close("red_threshold");

		//// analyze somata
		selectImage("red_somata");
		
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
		Roi.setStrokeWidth(1500);
		run("Straighten...");
		rename("red_somata_straight");

		// calculate box dimensions
		boxWidth = getWidth() / numCols;
		boxHeight = getHeight() / numRows;

		// loop through rows and columns to calculate coritcal position box intensities
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
		
		// calculate box dimensions for cortical layer analysis
		boxWidth_cl = getWidth() / numRows; // numRows = 1
		boxHeight_cl = getHeight() / numCols; // numCols = 500
		
		// loop through columns and rows to calculate cortical layer box intensities
		for (col = 0; col < numRows; col++) {
  			for (row = 0; row < numCols; row++) {
    			
 		     	// calculate box coordinates
 		     	boxX = 0 + col * boxWidth_cl;
  		    	boxY = 0 + row * boxHeight_cl;

   		   		// define ROI for current box
   		   		makeRectangle(boxX, boxY, boxWidth_cl, boxHeight_cl);

    		  	// calculate box intensity
    		  	run("Measure");
			}			
		}		
		close("red_somata_straight");

		// threshold image
		selectImage("red_somata");
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
		Roi.setStrokeWidth(1500);
		run("Straighten...");
		rename("red_somata_straight");

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
		
		// calculate box dimensions for cortical layer analysis
		boxWidth_cl = getWidth() / numRows; // numRows = 1
		boxHeight_cl = getHeight() / numCols; // numCols = 500
		
		// loop through columns and rows to calculate cortical layer box intensities
		for (col = 0; col < numRows; col++) {
  			for (row = 0; row < numCols; row++) {
    			
 		     	// calculate box coordinates
 		     	boxX = 0 + col * boxWidth_cl;
  		    	boxY = 0 + row * boxHeight_cl;

   		   		// define ROI for current box
   		   		makeRectangle(boxX, boxY, boxWidth_cl, boxHeight_cl);

    		  	// calculate box intensity
    		  	run("Measure");
			}			
		}		
		close("red_somata_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "red_somata";
			box_results[box] = box;
			pix_results[box] = boxWidth * (box + 1);
			dist_results[box] = boxWidth * (box + 1) * pixelWidth;
			mean_results[box] = getResult("Mean", box);
			threshMean_results[box] = getResult("Mean", box + total_boxes*2);
			threshFrac_results[box] = threshMean_results[box] / 255;
			threshArea_results[box] = getResult("Area", box + total_boxes*2);
			threshNumber_results[box] = threshArea_results[box] * threshFrac_results[box];
			height_results[box] = threshArea_results[box] / boxWidth;
		}

		// store arrays in table and save
		Table.create("red_somataTable");
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
		Table.save(dir + "/" + strip_name + "_red_somata.csv");
		close("red_somataTable");
		
		// make cortical layer arrays
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
		
		// store coritcal layer values in cortical layer arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "red_somataLayer";
			box_results[box] = box;
			pix_results[box] = boxHeight_cl * (box + 1);
			dist_results[box] = boxHeight_cl * (box + 1) * pixelHeight;
			mean_results[box] = getResult("Mean", box + total_boxes);
			threshMean_results[box] = getResult("Mean", box + total_boxes*3);
			threshFrac_results[box] = threshMean_results[box] / 255;
			threshArea_results[box] = getResult("Area", box + total_boxes*3);
			threshNumber_results[box] = threshArea_results[box] * threshFrac_results[box];
			height_results[box] = threshArea_results[box] / boxHeight_cl;
		}
		close("Results");
		close("red_somata");

		// store arrays in table and save
		Table.create("red_somataLayerTable");
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
		Table.save(dir + "/" + strip_name + "_red_somataLayer.csv");
		close("red_somataLayerTable");


		//// analyze proximal
		selectImage("red_proximal");

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
		rename("red_proximal_straight");

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
		close("red_proximal_straight");

		// threshold image
		selectImage("red_proximal");
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
		rename("red_proximal_straight");

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
		close("red_proximal_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "red_proximal";
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
		close("red_proximal");

		// store arrays in table and save
		Table.create("red_proximalTable");
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
		Table.save(dir + "/" + strip_name + "_red_proximal.csv");
		close("red_proximalTable");


		//// analyze distal
		selectImage("red_distal");

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
		rename("red_distal_straight");

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
		close("red_distal_straight");

		// threshold image
		selectImage("red_distal");
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
		rename("red_distal_straight");

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
		close("red_distal_straight");
		
		// store values in arrays
		for(box=0; box<total_boxes; box++) {
			name_results[box] = strip_name;
			feature_results[box] = "red_distal";
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
		close("red_distal");

		// store arrays in table and save
		Table.create("red_distalTable");
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
		Table.save(dir + "/" + strip_name + "_red_distal.csv");
		close("red_distalTable");
		close("original");
		close("ROI Manager");
	}
}

beep();
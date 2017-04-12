//console.log('test.js: started.');

// read and parse the json file in a javascript object.
// /!\ This process is asynchronous !
jQuery.getJSON("resultsDB.json", {}, function(json){
  // console output: useful for debugging purposes (indicates this function was called)
  // /!\ No output generally means the json file was not reached !
  //console.log("getJSON: started")
  // addDynamicTable will populate the html page. The argument is the javascript object.
  // This ensures the whole json file has already entirely been parsed !
  addDynamicTable(json);
});

function addDynamicTable(resultsDBJson) {
  // We first populate the statistics and run environment sections
  var cfgContent = document.createElement('p');
  cfgContent.id = "configDetails";
  document.getElementById('configSection').appendChild(cfgContent);

  // "Quick and dirty" method (print everything in bulk)
  /*
  for( cfg in resultsDBJson.configuration )
  {
    var cfgLine = document.createElement('div');
    cfgLine.appendChild(document.createTextNode(cfg + " : " + resultsDBJson.configuration[cfg] + ""));
    cfgContent.appendChild(cfgLine);
  }
  */

  // Line by line method: we control exactly what is displayed and how

  // We add the run statistics

  var statsContent = document.createElement('p');
  statsContent.id = "statsDetails";
  document.getElementById('statsSection').appendChild(statsContent);

  cfgLine = document.createElement('div');
  cfgKey = "found files";
  cfgLine.appendChild(document.createTextNode("Found files" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  statsContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "processed files";
  cfgLine.appendChild(document.createTextNode("Processed files" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  statsContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "successful reconstructions";
  cfgLine.appendChild(document.createTextNode("Successful reconstructions" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  statsContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "incomplete reconstructions";
  cfgLine.appendChild(document.createTextNode("Incomplete reconstructions" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  statsContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "total failure";
  cfgLine.appendChild(document.createTextNode("Total failures" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  statsContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "total time";
  cfgLine.appendChild(document.createTextNode("Total execution time" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  statsContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "processes";
  cfgLine.appendChild(document.createTextNode("Number of processes" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  statsContent.appendChild(cfgLine);


  // We add the run environment

  var cfgLine = document.createElement('div');
  cfgKey = "date";
  cfgLine.appendChild(document.createTextNode("Date" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  cfgContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "database";
  cfgLine.appendChild(document.createTextNode("Database (root directory)" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  cfgContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "system";
  cfgLine.appendChild(document.createTextNode("System details" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  cfgContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "user";
  cfgLine.appendChild(document.createTextNode("User" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  cfgContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "angiotk version";
  cfgLine.appendChild(document.createTextNode("AngioTK version (commit)" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  cfgContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "python";
  cfgLine.appendChild(document.createTextNode("Python version" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  cfgContent.appendChild(cfgLine);

  cfgLine = document.createElement('div');
  cfgKey = "paraview";
  cfgLine.appendChild(document.createTextNode("Paraview version" + " : " + resultsDBJson.configuration[cfgKey] + ""));
  cfgContent.appendChild(cfgLine);



  // ------------------------------ Results section ------------------------------


  var resultsTable = document.getElementById('results table');
  var resultsTableHeader = document.getElementById('results table head');
  var resultsTableBody = document.getElementById('results table body');

  var imageFileNames = resultsDBJson.results;
  //resultsTableHeader

  nPipelineSteps = 0;
  // get pipeline step names to create and fill table header columns
  for( imageName in imageFileNames )
  {
    for( step in imageFileNames[imageName])
    {
      nPipelineSteps++;
      headerCol = document.createElement('th');
      headerCol.appendChild(document.createTextNode(step));
      resultsTableHeader.appendChild(headerCol);
    }
    break; // We actually only need on image file here to get the step titles
  }

  // for each image file
  for( imageName in imageFileNames )
  {
    // create a title (image name) row in the table
    imageNameRow = document.createElement('tr');
    imageNameRow.className = "accordionHead";
    imageNameRow.valign-top
    // add a first column (or cell) with the image file name
    titleCell = document.createElement('td');
    titleCell.className = "titleCell";
    titleCell.colSpan = nPipelineSteps+1;
    //titleCell.setAttribute('colspan')
    titleCell.appendChild(document.createTextNode(imageName));
    imageNameRow.appendChild(titleCell);

    // create a row to display details of each pipeline step execution
    imageDataRow = document.createElement('tr');
    imageDataRow.className = "accordionBody";

    nFailedSteps = 0;
    // add a cell for each pipeline step
    for( step in imageFileNames[imageName] )
    {
      stepCellTD = document.createElement('td');
      stepCellTD.className = "accordionTD";
      stepCellDivContainer = document.createElement('div');
      stepCellDivContainer.className = "container";
      stepCellTD.appendChild(stepCellDivContainer);
      if(imageFileNames[imageName][step]["Success"] == false)
      {
        //stepCellDivContainer.setAttribute("background-color", "red");
        stepCellDivContainer.style.backgroundColor = "red";
        stepCellTD.style.backgroundColor = "red";
        nFailedSteps++;
      }
      var timeText = document.createElement('div');
      timeText.appendChild(document.createTextNode("Time" + ": " + imageFileNames[imageName][step]["Time"] + "s"));
      stepCellDivContainer.appendChild(timeText);

      // output link
      var output = document.createElement('a');
      output.appendChild(document.createTextNode("Output"));
      output.href = imageFileNames[imageName][step]["Output"];
      output.title = "link to output";
      stepCellDivContainer.appendChild(output);

      // screenshot link
      var screenshot = document.createElement('a');
      img = document.createElement('img');
      img.className = "thumbnail";
      imgFile = imageFileNames[imageName][step]["Screenshot"]; // screenshots/img.ext
      var imgPath, imgName;
      var n = imgFile.lastIndexOf('/'); // we want to split after the last '/' to have dirname and basename ('screenshots/' and 'img.ext')
      img.src = imgFile.substr(0, n+1) + "tn_" + imgFile.substr(n+1); // 'screenshots/tn_img.png': we want a thumbnail image linking to full size image
      // if no thumbnail image is found, we use the full size image (slower method)
      img.setAttribute( "onError", "this.onerror=null;this.src='" + imgFile + "';" );
      img.alt = "No screenshot"; // if neither image is found, we simply use a hyperlink
      screenshot.appendChild(img);
      screenshot.href = imgFile;
      screenshot.title = "See full size screenshot."
      stepCellDivContainer.appendChild(screenshot);

      // comments
      var commentText = document.createElement('div');
      commentText.appendChild(document.createTextNode("Comments" + ": " + imageFileNames[imageName][step]["Comments"]));
      stepCellDivContainer.appendChild(commentText);


      imageDataRow.appendChild(stepCellTD);
    }
    if(nFailedSteps > 0) {
      titleCell.appendChild(document.createTextNode(" - " + nFailedSteps + " failed steps"));
      //imageNameRow.style.backgroundColor = "rgba(255,0,0,0.1)";
      imageNameRow.className += " failureRow";
    }
    else {
      //imageNameRow.style.backgroundColor = "rgba(0,255,0,0.1)";
      imageNameRow.className += " successRow";
    }
    resultsTableBody.appendChild(imageNameRow);
    resultsTableBody.appendChild(imageDataRow);
  }

  // initial folding and click event (accordion effect) for title cells
  var heads = document.querySelectorAll(".accordionHead");
  var headsLength = heads.length;
  var i;
  for (i = 0; i < headsLength; i++) {
    $(heads[i]).next('.accordionBody').children().children().hide();
    $(heads[i]).click(function(){
      $(this).next('.accordionBody').children().children().slideToggle();
    });
  $(heads[0]).next('.accordionBody').children().children().show();
  }
}

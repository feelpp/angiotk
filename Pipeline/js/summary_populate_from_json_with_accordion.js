console.log('test.js started.');

// read and parse the json file in a javascript object.
// /!\ This process is asynchronous !
jQuery.getJSON("summary.json", {}, function(json){
  // console output: useful for debugging purposes (indicates this function was called)
  // /!\ No output generally means the json file was not reached !
  console.log('Number of runs: ' + Object.keys(json.runs).length);
  // addDynamicTable will populate the html page. The argument is the javascript object.
  // This ensures the whole json file has already entirely been parsed !
  addDynamicTable(json);
});

console.log('getJSON passed');

function addDynamicTable(summaryJson) {
  // console output: useful for debugging purposes (indicates this function was called)
  console.log('addDynamicTable started');
  var pageContent = document.getElementById('content');

  var runCountLine = document.createElement('div');
  runCountLine.appendChild(document.createTextNode(Object.keys(summaryJson.runs).length + " runs:"));
  pageContent.appendChild(runCountLine);

  var accordion = document.createElement('div');
  accordion.className = "master_p";
  accordion.id = "accordion";
  //accord.appendChild(document.createTextNode("Paragraphe du tableau"));
  pageContent.appendChild(accordion);

  for( run in summaryJson.runs )
  {
    // Each run is put in a section. This is the section's header.
    var sectionHeader = document.createElement('h3');
    // header color is determined by errors (incomplete/failed reconstructions)
    processed = summaryJson.runs[run]["processed files"];
    successfuls = summaryJson.runs[run]["successful reconstructions"];
    incompletes = summaryJson.runs[run]["incomplete reconstructions"];
    failures = summaryJson.runs[run]["total failure"];
    hints = []
    if ( successfuls != processed )
    {
      if ( (processed > 0) && (failures == processed) )
      {
        sectionHeader.style.color = "red"
        hints.push("total failure ");
      }
      else
      {
        sectionHeader.style.color = "orange"
        if ( failures > 0 )
        {
          hints.push("failures ");
        }
        if ( incompletes > 0 )
        {
          hints.push("incomplete reconstructions ");
        }
	if ( (processed > successfuls + incompletes + failures) )
	{
	    hints.push("aborted ");
	}
      }
    }
    else
    {
      sectionHeader.style.color = "Limegreen"
    }

    // header title building.
    headerTitle = run;
    if (hints.length > 0) headerTitle += ' ( ' + hints.join(', ') + ')';
    sectionHeader.appendChild(document.createTextNode(headerTitle));
    accordion.appendChild(sectionHeader);

    var sectionDiv = document.createElement('div');
    accordion.appendChild(sectionDiv);

    // This paragraph contains all the details about the run
    var sectionDivParagraph = document.createElement('p');
    sectionDiv.appendChild(sectionDivParagraph);

    for( key in summaryJson.runs[run] )
    {
      if ( key != 'link' && key != 'command lines' )
      {
        sectionDivParagraph.appendChild(document.createTextNode(key + ": " + summaryJson.runs[run][key]));
        sectionDivParagraph.appendChild(document.createElement('br'));
      }
    }
    var a = document.createElement('a');
    a.href = summaryJson.runs[run]['link'];
    a.appendChild(document.createTextNode("link to full page"));
    sectionDivParagraph.appendChild(a);
  }
  
}

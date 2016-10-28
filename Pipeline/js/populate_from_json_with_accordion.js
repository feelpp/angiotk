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

  /*
  var parag = document.createElement('div');
  parag.className = "master_p";
  parag.id = "paragraphe contenant le tableau";
  parag.appendChild(document.createTextNode("Paragraphe du tableau"));
  pageContent.appendChild(parag);

  for( run in summaryJson.runs )
  {
    var dataParagraph = document.createElement('div');
    dataParagraph.className = "data paragraph";
    pageContent.appendChild(dataParagraph);
    // colonne gauche: nom du run
    var leftDiv = document.createElement('div');
    leftDiv.className = "left column";
    leftDiv.style = "float:left; width:50%";
    dataParagraph.appendChild(leftDiv);

      var a = document.createElement('a');
      a.href = summaryJson.runs[run]['link'];
      a.appendChild(document.createTextNode(run));
      leftDiv.appendChild(a);
      leftDiv.appendChild(document.createElement('br'));

    // colonne droite: détails du run
    var rightDiv = document.createElement('div');
    rightDiv.className = "right column";
    rightDiv.style = "float:right; width:50%; margin-left:10px";
    dataParagraph.appendChild(rightDiv);

    for( key in summaryJson.runs[run] )
    {
      if ( key != 'link' )
      {
        rightDiv.appendChild(document.createTextNode(key + ": " + summaryJson.runs[run][key]));
        rightDiv.appendChild(document.createElement('br'));
      }
    }
  }
  */
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  /*
  // on créé le tableau
  var table = document.createElement('table');
  table.id = "Tableau créé dynamiquement.";
  table.className = "tableblock frame-all grid-all spread";
  pageContent.appendChild(table);
  // on créé un <caption> pour le tableau
  var caption = document.createElement('caption');
  caption.className = "title";
  caption.appendChild(document.createTextNode('Tableau dynamique:'));
  table.appendChild(caption);
  //on créé un <colgroup> pour le tableau
  var colzGroup = document.createElement('colgroup');
  // on lui intègre 2 colonnes
  for (var i = 0; i < 2; i++) {
    var column = document.createElement('col');
    column.style = "width: 50%;";
    colzGroup.appendChild(column);
  }
  table.appendChild(colzGroup);
  // on créé un <thead> pour le tableau,
  // c'est la ligne des titres de colonnes.
  var tHead = document.createElement('thead');
  table.appendChild(tHead);
  var tRow = document.createElement('tr');
  tHead.appendChild(tRow);
  // on créé un titre par colonne.
  for (var i = 0; i < 2; i++) {
    var tH = document.createElement('th');
    tH.className = "tableblock halign-left valign-top";
    if (i==0) tH.appendChild(document.createTextNode('Run ID'));
    else tH.appendChild(document.createTextNode('Details'));
    tRow.appendChild(tH);
  }

  // on créé un <tbody> pour le tableau
  var tBody = document.createElement('tbody');
  table.appendChild(tBody);
  // on créé les autre lignes du tableau,
  // celles qui contiennent les informations
  // sur les runs.
  for( run in summaryJson.runs )
  {
    var tRow = document.createElement('tr');
    tBody.appendChild(tRow);
    // colonne gauche: nom du run
    var tD = document.createElement('td');
    tD.className = "tableblock halign-left valign-top";
    var paragraph = document.createElement('p');
    tD.appendChild(paragraph);

      var a = document.createElement('a');
      a.href = summaryJson.runs[run]['link'];
      a.appendChild(document.createTextNode(run));
      paragraph.appendChild(a);
      paragraph.appendChild(document.createElement('br'));

    //paragraph.appendChild(document.createTextNode(run));
    tRow.appendChild(tD);
    // colonne droite: détails du run
    var tD = document.createElement('td');
    tD.className = "tableblock halign-left valign-top";
    var paragraph = document.createElement('p');
    paragraph.className = "testRow";
    tD.appendChild(paragraph);
    for( key in summaryJson.runs[run] )
    {
      if ( key != 'link' )
      {
        paragraph.appendChild(document.createTextNode(key + ": " + summaryJson.runs[run][key]));
        paragraph.appendChild(document.createElement('br'));
      }
    }
    tRow.appendChild(tD);
  }
  */

  /*
  $('.master_p').readmore(
    {
      collapsedHeight: 10,
      speed: 100,
      moreLink: '<a href="#">Read more.</a>'
    });
    */
}

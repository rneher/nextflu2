var tip_labels = true, branch_labels=false;
function treePlotHeight(width) {return 400 + 0.30*width;}
var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
var treeWidth = containerWidth;
var treeHeight = treePlotHeight(treeWidth);
var treeplot = d3.select("#treeplot")
    .attr("width", treeWidth)
    .attr("height", treeHeight);

var legendCanvas = d3.select("#legend")
    .attr("width", 280)
    .attr("height", 140);

var mutType='aa'; //mutations displayed in tooltip

var serumSymbol = '\uf0fe';
var focusSymbol = '\uf05b'
var focusNode;
var decorations=true;

function tipDecoration(d) {
    if (typeof(d.serum) === 'undefined'){
        return '';
    }else if (d.serum){
        if (d==focusNode) 
            {return focusSymbol;} 
        else 
            {return serumSymbol;}
    }else{
        return '';
    }
}


treeplot.left_margin = 10;
treeplot.bottom_margin = 16;
treeplot.top_margin = 32;
if (branch_labels) {treeplot.top_margin +=15;}
treeplot.right_margin = 10;

var regions = ["Africa", "SouthAmerica", "WestAsia", "Oceania", "Europe", "JapanKorea",
               "NorthAmerica", "SoutheastAsia", "SouthAsia", "China"]
var LBItau=0.0005;

function load_tree(){
    var myTree;
    var myTreeSearch;
    var myLegend;
    var cladeToSeq;
    var myDivChart;
    var myFreqChart;
    var oneYear = 365.25*24*60*60*1000; // days*hours*minutes*seconds*milliseconds
    var tw = 10.0;
    var colorScale;

    function legend_mouseover(legend_element){
        var lb = myLegend.lowerBound[legend_element];
        var ub = myLegend.upperBound[legend_element];
        if (typeof lb == "number"){
            myTree.tips.forEach(function (d){d._highlight = (lb<=d.coloring && ub>d.coloring);});
        }else{
            myTree.tips.forEach(function (d){d._highlight = legend_element==d.coloring;});
        }
        myTree.updateStyle(colorScale, branchColor="same");
        console.log("Mousover, updated radius " + legend_element +" "+ lb+" "+ub);
    }
    function legend_mouseout(){
        myTree.tips.forEach(function (d){d._highlight = false;});
        myTree.updateStyle(colorScale, branchColor="same");
    };

    function legend_click(legend_element){
        var lb = myLegend.lowerBound[legend_element];
        var ub = myLegend.upperBound[legend_element];
        if (typeof lb == "number"){
            myTree.tips.forEach(function (d){d._selected = d._selected||(lb<=d.coloring && ub>d.coloring);});
        }else{
            myTree.tips.forEach(function (d){d._selected = d._selected||(legend_element==d.coloring);});
        }
        myTree.updateStyle(colorScale, branchColor="same");
        console.log("Click, updated style " + legend_element +" "+ lb+" "+ub);
    }

/********************************************************************************************
TITER STUFF
*********************************************************************************************/
    function addTipDecorations(){
        console.log('adding tip decorations');
        treeplot.selectAll(".tipDeco")
            .data(myTree.nodes.filter(function(d) {return d.serum;}))
            .enter()
            .append("text")
            .style('font-family', 'FontAwesome')
            .attr('text-anchor', 'middle')
            .attr('dominant-baseline', 'central')
            .attr("class","tipDeco")
            .style("cursor", "pointer")
            .text(tipDecoration)
            .style("font-size", function (d) {if (d==focusNode) {return "30px";} else {return "12px";}})
            .on('click', changeTiterFocus)
            .on('mouseover', function(d) {
                virusTooltip.show(d, this);
            })
            .on('mouseout', virusTooltip.hide);
        myTree.updateGeometry();
    }

    function updateTipDecoration(){
        console.log('updating tip decorations');
        treeplot.selectAll(".tipDeco")
            .text(tipDecoration)
            .style("font-size", function (d) {if (d==focusNode) {return "30px";} else {return "12px";}});
    }

    function changeTiterFocus(d) {
        focusNode=d;
        updateTipDecoration();
        assignRawTiter(myTree.tips, focusNode);
        colorByTiter();
    }

    function colorByTiter(){
        myTree.nodes.forEach(function (d){d.coloring = d.HI_dist; d._valid=(d.HI_dist!='NaN');});
        var activeTips = myTree.tips.filter(function (d){return d.current;});
        var activeVals = activeTips.map(function(d){return (d.HI_dist=='NaN')?0:d.HI_dist;});
        d3.select("#legend-title").html(function(d){return "HI titers";});
        colorScale = myTree.currentColorScale;
        var curr_min = d3.min(activeVals);
        var curr_max = d3.max(activeVals);
        colorScale.domain(myTree.zero_one.map(function (d){return curr_min + d*(curr_max-curr_min);}));

        myTree.updateStyle(colorScale, branchColor="#DDD");
        if (typeof myLegend != "undefined"){
            myLegend.remove();
        }
        var legend_stack_height = d3.max([10,colorScale.domain().length])/2;
        legendCanvas.attr('height',20*(legend_stack_height+1));
        var label_fmt = function(d) {return d.toFixed(2).replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');}
        myLegend = new legend(legendCanvas, colorScale, label_fmt,
                              legend_mouseover, legend_mouseout, legend_click,
                              legend_stack_height);
    }


/********************************************************************************************
GENERIC COLORING
*********************************************************************************************/
    function updateColor(){
        focusNode='undefined';
        updateTipDecoration();
        var choice = document.getElementById("coloring").value;
        console.log("coloring by", choice);
        myTree.nodes.forEach(function (d){d._valid=true;});
        if (choice=="date"||choice=="ep"||choice=="ne"||choice=="rb"||choice=="cTiter"||choice=="LBI"){
            if (choice=="LBI")
                {calcLBI(myTree.rootNode, myTree.nodes);}
            
            if (choice=="date") {myTree.nodes.forEach(function (d){d.coloring = d._numDate;});}
            else {myTree.nodes.forEach(function (d){d.coloring = d[choice];});}
            var activeTips = myTree.tips.filter(function (d){return d.current;});
            var activeVals = activeTips.map(function(d){return (choice=="date")?d._numDate:d[choice];});

            d3.select("#legend-title").html(function(d){return "Collection date";});
            var label_fmt = function(d) {return d.toFixed(2).replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');}
            colorScale = myTree.currentColorScale;
            var curr_min = d3.min(activeVals);
            var curr_max = d3.max(activeVals);
            colorScale.domain(myTree.zero_one.map(function (d){return curr_min + d*(curr_max-curr_min);}));
        }else{
            myTree.nodes.forEach(function (d){d.coloring = d[choice];});
            var cats = [];
            myTree.nodes.forEach(function (d){cats.push(d[choice]);});
            var tmp_categories = d3.set(cats).values();
            d3.select("#legend-title").html(function(d){return (choice=="region")?"Region":"Country";});
            var label_fmt = function(d) {return d;}
            var colorScale = d3.scale.ordinal()
                .domain(tmp_categories)
                .range(genotypeColors);
        }

        myTree.updateStyle(colorScale, branchColor="same");
        if (typeof myLegend != "undefined"){
            myLegend.remove();
        }
        var legend_stack_height = d3.max([10,colorScale.domain().length])/2;
        legendCanvas.attr('height',20*(legend_stack_height+1));
        myLegend = new legend(legendCanvas, colorScale, label_fmt,
                              legend_mouseover, legend_mouseout, legend_click,
                              legend_stack_height);

    }

    d3.json(file_prefix + "tree.json", function (error, root){
        if (error) return console.warn(error);

        document.getElementById("timetree").checked=false;
        myTree = new PhyloTree(root, treeplot, d3.select('.treeplot-container'));
        addTipDecorations();
        d3.select("#reset").on("click", myTree.resetLayout);
        d3.select("#timetree").on("change", function(){
            var tmp_timetree = document.getElementById("timetree").checked;
            myTree.toggleTimeTree(tmp_timetree);
        });

        var mc = autocomplete(document.getElementById('search'))
            .keys(myTree.tips)
            .dataField("strain")
            .placeHolder("search strains...")
            .width(800)
            .height(500)
            .onSelected(onSelect)
            .render();

        var searchEvent;
        function onSelect(tip) {
            var strainName = (tip.strain).replace(/\//g, "");
            d3.select("#"+strainName)
                .call(function(d) {
                    virusTooltip.show(tip, d[0][0]);
                })
                .style("fill", function (d) {
                  d._selected = true;
                  searchEvent = setTimeout(function (){
                  virusTooltip.hide(d);}, 5000, d);
                  return d3.rgb(tipFillColor(d)).brighter();
                });
            myTree.updateStyle(colorScale, branchColor="same");
        }


        function draggedMinFun(start, end){
            myTree.setNodeState(start, end);
            myTree.updateVisibility();
        }

        function draggedEndFunc(start, end){
            myTree.setNodeState(start, end);
            updateColor();
            myTree.updateAnnotations();
        }

        function draggedFunc(start, end){
            myTree.setNodeState(start, end);
            myTree.updateVisibility();
        }
        myDateSlider = new dateSlider(draggedFunc, draggedMinFun, draggedEndFunc);
        myDateSlider.date_init(myTree.earliestDate, myTree.latestDate, tw);
        var label_fmt = function(d) {return d.toFixed(2).replace(/([10-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');}
        myLegend =  new legend(legendCanvas, myTree.currentColorScale, label_fmt,
                               legend_mouseover, legend_mouseout);
    });


    function stateAtPosition(clade, gene, pos){
        if (typeof cladeToSeq[clade][gene][pos] == "undefined"){
            return cladeToSeq["root"][gene][pos];
        }else{
            return cladeToSeq[clade][gene][pos];
        }
    }

    d3.json(file_prefix + "sequences.json", function(error, json) {
        if (error) return console.warn(error);
        cladeToSeq=json;
        var searchClades = myTree.tips.map(function (d){return d.clade;});
        myTreeSearch = new TreeSearch(stateAtPosition, myTree.tips.map(function(d){return d.clade;}),
                                        cladeToSeq['root']['nuc'], seqSearchResult);
    });


    d3.select("#coloring").on("change", updateColor);

    var genotypeColoringEvent;
    d3.select("#gt-color")
        .on("keyup", function(){
            if (typeof genotypeColoringEvent != "undefined"){clearTimeout(genotypeColoringEvent);}
            genotypeColoringEvent = setTimeout(colorByGenotype, 200);
        });

    var genotypeColors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"];
    function colorByGenotype() {
        document.getElementById("coloring").value = 'none';
        var user_input = document.getElementById("gt-color").value.split(':')
        if (user_input.length==1){
            var pos = parseInt(user_input[0])-1;
            var gene=default_gene;
        }else{
            var pos = parseInt(user_input[1])-1;
            var gene=user_input[0];
        }
        colorByPosition(pos, gene);
    };

    function colorByPosition(pos, gene){
        var gts = myTree.nodes.map(function (d) {
            d.coloring = stateAtPosition(d.clade, gene, pos);
            return d.coloring;});

        var tmp_categories = d3.set(gts).values();
        var tmp_range = [], tmp_freqs = [];
        for (var ii=0; ii<tmp_categories.length; ii++)
            {tmp_range.push(genotypeColors[ii%genotypeColors.length]);
             tmp_freqs.push('global'+'_'+gene+':'+(pos+1)+tmp_categories[ii]);}
        colorScale = d3.scale.ordinal()
            .domain(tmp_categories)
            .range(tmp_range);
        myTree.updateStyle(colorScale, branchColor="same");
        myLegend.remove();
        var label_fmt = function(d) {return d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');}
        d3.select("#legend-title").html(function(d){return "State at " + gene+' '+(pos+1);});
        myLegend =  new legend(legendCanvas, colorScale, label_fmt,
                                legend_mouseover, legend_mouseout, legend_click);

        console.log(tmp_range);
        myFreqChart.addTrajectories(tmp_freqs, tmp_range);
    }
    this.colByPos = colorByPosition;

    // callback to highlight the result of a search by strain name
    var searchEvent;
    function highlightStrainSearch(tip) {
        var strainName = (tip.strain).replace(/\//g, "");
        d3.select("#"+strainName)
            .call(function(d) {
                markInTreeStrainSearch(tip);
                virusTooltip.show(tip, d[0][0]);
            });
    }

    var strainSearchEvent;
    d3.select('#seqinput').on('keyup', function(){
        if (typeof strainSearchEvent != "undefined"){clearTimeout(strainSearchEvent);}
        var lines = document.getElementById('seqinput').value.split('\n');
        strainSearchEvent = setTimeout(function (d) {myTreeSearch.parseSequences(lines);}, 500);
    });

    function seqSearchResult(clades){
        var nodesToHighlight = myTree.nodes.filter(function (d){
        var tmp=0;
        for (var clade in clades){tmp+= (d.clade==clade);}
            return tmp>0;});
        console.log(clades, nodesToHighlight);
        markSearchResult(nodesToHighlight);
    }

    // highlight clades in tree
    function markSearchResult(nodesToHighlight){
        treeplot.selectAll(".searchResult").data(nodesToHighlight)
            .enter()
            .append('text')
            .attr("class", "searchResult")
            .attr('text-anchor', 'middle')
            .attr('dominant-baseline', 'central')
            .style("font-size", "24px")
            .text(function(d) {return '\uf069'; })
            .on('mouseover', function(d) {virusTooltip.show(d, this);})
            .on('mouseout', virusTooltip.hide);
        myTree.updateStyle(colorScale, branchColor="same");
        myTree.updateGeometry(0.0);
    }

    d3.select('#searchinputclear').on('click', function (){
        myTree.tips.forEach(function(d) {d._selected=false;});
        //treeplot.selectAll('.searchResult').data([]).exit().remove();
        document.getElementById('seqinput').value = "";
        virusTooltip.hide();
        myTree.updateStyle(colorScale, branchColor='same');
    });

    function diversityCallback(d){
        console.log(d);
        colorByPosition(d.index, d.name);
    }

    d3.json(file_prefix+"entropy.json", function(error, S){
        var chart_data = {};
        for (var ii=0; ii<proteins.length; ii++){
            chart_data[proteins[ii]] = S[proteins[ii]]
        }
        myDivChart = new diversityChart(d3.select('.entropy-container'), '#entropy',
                                        chart_data, diversityCallback)
    });

    d3.json(file_prefix+"frequencies.json", function(error, freqs){
        myFreqChart = new frequencyChart(d3.select('.frequency-container'), '#frequencies',
                                         freqs['pivots'],freqs);
        myFreqChart.addTrajectories(['global_HA1:159Y','global_HA1:159F', 'global_HA1:159S'], genotypeColors);
    });

}

auspice_tree = load_tree();


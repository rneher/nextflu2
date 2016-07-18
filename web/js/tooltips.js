var virusTooltip = d3.tip()
	.direction(function(d){return (d.x<600)?'e':'w';})
	.attr('class', 'd3-tip')
	.offset(function (d){ return (d.x<600)?[0,12]:[0,-12];})
	.html(function(d) {

		string = "";

		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";

		string += "<div class=\"smallnote\">";

		// check if vaccine strain
		if (vaccineStrains.indexOf(d.strain) != -1) {
			string += "Vaccine strain<br>";
			var vaccine_date = new Date(vaccineChoice[d.strain]);

			string += "First chosen " + vaccine_date.toLocaleString("en-us", { month: "short" }) + " " + vaccine_date.getFullYear() + "<br>";
			string += "<div class=\"smallspacer\"></div>";
		}

		if (typeof d.country != "undefined") {
			string += d.country.replace(/([A-Z])/g, ' $1');
		}
		if (typeof d.date != "undefined") {
			string += ", " + d.date;
		}
		if (typeof d.isolate_id != "undefined") {
			string += "<br>GISAID ID: " + d.isolate_id;
		}
		if (typeof d.lab != "undefined") {
			if (d.lab != "") {
				string += "<br>Source: " + d.lab.substring(0,25);
				if (d.lab.length>25) string += '...';
			}
		}
		string += "</div>";
		string += "<div class=\"smallspacer\"></div>";
		// following may or may not be present
		string += "<div class=\"smallnote\">";
		if (typeof d.cTiter != "undefined") {
			string += "Antigenic adv: " + d.cTiter.toFixed(1) + "<br>";
		}
		if (typeof d.ep != "undefined") {
			string += "Epitope distance: " + d.ep + "<br>";
		}
		if (typeof d.rb != "undefined") {
			string += "Receptor binding distance: " + d.rb + "<br>";
		}
		if (typeof d.LBI != "undefined") {
			string += "Local branching index: " + d.LBI.toFixed(3) + "<br>";
		}
		if (typeof d.dfreq != "undefined") {
			string += "Freq. change: " + d.dfreq.toFixed(3) + "<br>";
		}
		if (typeof d.fitness != "undefined") {
			string += "Fitness: " + d.fitness.toFixed(3) + "<br>";
		}
		if (typeof d.pred_distance != "undefined") {
			string += "Predicted distance: " + d.pred_distance.toFixed(3) + "<br>";
		}
		string += "</div>";

		return string;
	});


var linkTooltip = d3.tip()
	.direction(function(d){return (d.x<600)?'e':'w';})
	.attr('class', 'd3-tip')
	.offset(function (d){ return (d.x<600)?[0,12]:[0,-12];})
	.html(function(d) {
		string = ""
		if (typeof d.frequency != "undefined") {
			string += "Frequency: " + (100 * d.frequency).toFixed(1) + "%"
		}
		string += "<div class=\"smallspacer\"></div>";
		string += "<div class=\"smallnote\">";
		if ((typeof d.aa_mut_str !="undefined")&&(mutType=='aa')){
			var ncount = 0;
			for (tmp_gene in d.aa_mut_str) {ncount+=d.aa_mut_str[tmp_gene].length;}
			if (ncount) {string += "<b>Mutations:</b><ul>";}
			for (tmp_gene in d.aa_mut_str){
				if (d.aa_mut_str[tmp_gene].length){
					string+="<li>"+tmp_gene+":</b> "+d.aa_mut_str[tmp_gene].replace(/,/g, ', ') + "</li>";
				}
			}
		}
		if ((typeof d.mut_str !="undefined")&&(d.mut_str.length)){
			var tmp_muts = d.mut_str.split(',');
			var nmuts = tmp_muts.length;
			tmp_muts = tmp_muts.slice(0,Math.min(10, nmuts))
			string += "<li>"+tmp_muts.join(', ');
			if (nmuts>10) {string+=' + '+ (nmuts-10) + ' more';}
			string += "</li>";
		}
		string += "</ul>";
		if (typeof d.region != "undefined") {
			string += "Region: "+d.region.replace(/([A-Z])/g, ' $1') + "<br>";
		}
		if (typeof d.dTiter != "undefined") {
			string += "Titer drop: "+d.dTiter.toFixed(1) + "<br>";
		}
		if (typeof d.cTiter != "undefined") {
			string += "cum. titer: "+d.cTiter.toFixed(1) + "<br>";
		}
		string += "click to zoom into clade"
		string += "</div>";
		return string;
	});

var raw_titers;
var tree_model;
var substitution_model;

d3.json(file_prefix + "titers.json", function (error, json){
    if (error) return console.warn(error);
    raw_titers = json;
});

d3.json(file_prefix + "titer_tree_model.json", function (error, json){
    if (error) return console.warn(error);
    tree_model = json;
});

d3.json(file_prefix + "titer_subs_model.json", function (error, json){
    if (error) return console.warn(error);
    substitution_model = json;
    console.log('Subsmodel:',substitution_model);
});


function assignRawTiter(tips, focusNode){
    var tmp_titers = raw_titers[focusNode.clade];
    for (var i=0; i<tips.length; i+=1){
        d = tips[i];
        if (typeof(tmp_titers[d.clade])!="undefined"){
            var tmp_HI=0;
            var serum_count=0;
            for (var tmp_serum in tmp_titers[d.clade]){
                tmp_HI += tmp_titers[d.clade][tmp_serum];
                serum_count+=1;
            }
            if (serum_count){
                d.HI_dist = tmp_HI/serum_count
            }else{
                d.HI_dist = 'NaN';             
            }
        }else{
            d.HI_dist = 'NaN';
        }
        console.log(d.HI_dist);
    }
};

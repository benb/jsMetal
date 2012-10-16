// GETHOMOLOGYSETS
var SSP = 0;
var SIM = 1;
var POS = 2;
var EVO = 3;
onmessage = function(e){
        var o = JSON.parse(e.data);
        getHomologySets(o.aln,o.tree,o.doEvo,o.seqNum);
}

function getHomologySets(aln,tree,doEvo,seqNum){	
	
	var aln=resortNonOverlapping(aln,seqNum);
	labeller(aln,tree,doEvo,seqNum);
	var gapsHere=[];
	var homologySets = [];
		
	for(var hom=0;hom<=POS+doEvo;hom++){
		
		homologySets[hom]=[];
		for (var i=0;i<aln.length;i++){
			homologySets[hom][i]=[];
			jNoGap=0;
			
			for ( var j=0;j<aln[i].content.length;j++){
				
				if(aln[i].content[j]  != "-"){
					
					homologySets[hom][i][jNoGap]=[];
					for(var k=0;k<aln.length;k++){
						if(aln[k].content[j]  == "-"){
							gapsHere[jNoGap]=true;
						}
						
						if(k!=i && aln[k].content[j]){
							
							homologySets[hom][i][jNoGap].push(aln[k].labeledContent[hom][j]);
							
						}
					
					
					}
					jNoGap++;
				}
					
					
			}
		}
	}
		
	postMessage(JSON.stringify([homologySets,gapsHere]));
}



// GLOBALBUBBLESORT

function resortNonOverlapping(alignment,seqNum) {
		var START=new Date();
	
	var repeat;
	var temp;
	do{
		repeat = false;
		
		for (var j = 0; j < (alignment[0].content.length - 1); j++){
			
			if( areNonOverlapping( j,alignment,seqNum) ){
				
				if( shouldBeFlipped( j,alignment,seqNum) ){
					
					for(var i=0;i<seqNum;i++){
						var newSeq= alignment[i].content.substr(0,j).concat(alignment[i].content[j+1]).concat(alignment[i].content[j]).concat(alignment[i].content.substr(j+2));
						alignment[i].content = newSeq;
						
						}
					repeat=true;
					
				}
			}
		}
	
	}while(repeat);
	var  END=new Date();
	
	return alignment;
}

// ARENONOVERLAPPING
// Tests if two columns are non-overlapping. We assume they are until
// we find a pair of non-gap characters.
function areNonOverlapping(j,aln,seqNum){
	
	var nonOverlapping = true;
	
	for (var i =0; i<seqNum;i++){
		
		//The moment we find a pair where neither character is gap, the columns are non-overlapping
		if( aln[i][j] != "-" && aln[i][j+1] != "-"){
			nonOverlapping = false;
			break;
		}
		
	}
	
	return nonOverlapping;
	
}

// SHOULDBEFLIPPED
// Tests if two *confirmed non-overlapping* columns should be flipped. Will do weird things
// if passed overlapping columns.
function shouldBeFlipped(j,aln,seqNum){
	
	var flipThem = false;
	
	for(var i=0; i<seqNum;i++){
		
		// Are both characters gaps? Let's continue
		if(aln[i][j] == "-" && aln[i][j+1]=="-"){
			continue;
		}
		// No? Then only one of them must be a gap (or something has gone horribly wrong).
		// If it's the one on the right, we must flip the columns. Whatever the case, we break the loop.
		else{
			flipThem = (aln[i][j+1]=="-")
			break;
		}		
	}
	return flipThem;
}

function labeller(alignment,tree,doEvo,seqNum){	
	var index;
	var nextLabel;
	var gapsHere=[];
	for(var i =0; i<seqNum;i++){
		alignment[i].labeledContent=[];
	}
	
	if(doEvo){
		evoLabeller(alignment,tree);
	
	}
	
	for(var i =0; i<seqNum;i++){
		alignment[i].labeledContent[SSP] = [];
		alignment[i].labeledContent[SIM] = [];
		alignment[i].labeledContent[POS] = [];
		
		
		index=0;
			
		
		for(var j=0;j<alignment[i].content.length;j++){
			
			if(alignment[i].content[j] != "-"){
				// Label character and increase index. Using pre-increment on index to start at 1 and thus allow gaps
				// that appear before any character to be labelled as 0.
				nextLabel = i + "X" + ++index;
				
				alignment[i].labeledContent[SSP].push(nextLabel);
				alignment[i].labeledContent[SIM].push(nextLabel);
				alignment[i].labeledContent[POS].push(nextLabel);
				
				if(doEvo){
				
					alignment[i].labeledContent[EVO][j]=nextLabel;
				}
				
				
			
			}
		
			else{
				//gapsHere[j]=true;
				// Do not label gaps
				alignment[i].labeledContent[SSP].push("-");
				// Label gaps by sequence
				alignment[i].labeledContent[SIM].push( i + "-");
				// Label gaps by position
				alignment[i].labeledContent[POS].push( i + "-" + index);
				// Add position information to evo-labelled gaps.
				if(doEvo){
					alignment[i].labeledContent[EVO][j]=alignment[i].labeledContent[EVO][j].concat(i + "-" + index);
					}
				
			}
		}
	}
	//return gapsHere;
}


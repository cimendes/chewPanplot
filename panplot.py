import csv, argparse, time, sys, os
import pandas
import mpld3
from mpld3 import utils, plugins 
import matplotlib.pyplot as plt
plt.switch_backend('agg') 


class ClickInfo2(plugins.PluginBase):
    """Plugin for getting info on click"""
    
    JAVASCRIPT = """
    mpld3.register_plugin("clickinfo2", ClickInfo2);
    ClickInfo2.prototype = Object.create(mpld3.Plugin.prototype);
    ClickInfo2.prototype.constructor = ClickInfo2;
    ClickInfo2.prototype.requiredProps = ["id"];
    ClickInfo2.prototype.defaultProps = {labels:null}   function ClickInfo2(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    ClickInfo2.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        labels = this.props.labels;
        obj.elements().on("mousedown",
                          function(d, i){ 
                            window.open(labels[i], '_blank')});
    }
    """
    def __init__(self, points, labels):
        self.points = points
        self.labels = labels
        """if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None"""
        self.dict_ = {"type": "clickinfo2","id": utils.get_id(points, "pts"),"labels": labels}

def distributionGraph(filename, outputname):
	#function to obtain the interactive gene frequency graph for the pangenome, in html
	genes=[]
	coreSize=0
	isolates=[]
	with open(filename, 'r') as fh:

		reader = csv.reader(fh, delimiter='\t')
		isolates=reader.next()[1:] 

		for row in reader:
			geneName=row[0]
			numbers = [ int(x) for x in row[1:]]
			freq=sum(numbers)
			genes.append((geneName, freq))
			if freq == len(isolates):
				coreSize+=1
	genes.sort(key=lambda tup: tup[1], reverse=True) 
	x=[]
	y=[]
	xLabels=[]

	for item in genes:
		y.append(item[1])
		x.append(genes.index(item))
		xLabels.append(item[0])

	fig,ax = plt.subplots(figsize=(12, 9))  

	line=ax.plot(x, y, '.', color="#3F5D7D")
	plt.title("Pan-Genome Frequency")
	plt.ylabel('Frequency')
	plt.xlabel('Gene')
	plt.text(1.1, 0.9,  s="%s Isolates | Pan-Genome: %s genes | Core Genome: %s genes | Acessory Genome: %s genes" % (str(len(isolates)),str(len(genes)), str(coreSize), str(len(genes)-coreSize)),fontsize=10, horizontalalignment='right',verticalalignment='center', transform = ax.transAxes) 
	ax.yaxis.labelpad = 40

	mpld3.plugins.connect(fig, plugins.PointLabelTooltip(line[0],labels=xLabels))

	mpld3.save_html(fig,outputname)



def main():

	parser = argparse.ArgumentParser(description='Generation pan-genomefrequency plot for a presence and absence profile file.', epilog='by C I Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-i', '--input', help='Inpu presence and absence profile file.')
	parser.add_argument('-o', '--output', help='Output file name.')

	args = parser.parse_args()

	print 'transposing profile...'
	df=pandas.read_csv(args.input, sep='\t', header=0, index_col=0)
	transpose=df.transpose()
		
	with open('clean_gene_presence_absence.tsv', 'w') as cleanPA:
		cleanPA.write('GENE')
		transpose.to_csv(cleanPA, sep='\t')

	print 'generating pan-genome frequency plot..'

	distributionGraph('clean_gene_presence_absence.tsv', args.output +'.html')

if __name__ == "__main__":
    main()
